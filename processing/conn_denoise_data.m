function conn_denoise_data(s,TR,task,derivDir,mask)
%% -------------------------------------------------------------------- %%
% This script denoises MNI-space functional time series by running CONN's
% confound regression pipeline. Also uses SPM and cosmoMVPA tools
%
% conn saves the denoised voxel level 4D data in:
% /results/preprocessing/niftiDATA_Subject*_Condition*.nii, which is
% moved and renamed to ~data/denoised/. Saves for a specific task,
% (condition), with runs concatenated
%
% After denoising, the script loads back in the .niis to convert zeros
% (outside of brain mask) to NaNs -- this is needed to maintaining the
% masking if the data is subsequently smoothed.
%
% Also, although you can give conn zipped files (.nii.gz), it will just go
% ahead and unzip them in their directory before running anything, so the
% unzipped version can be removed at the end (section 8)
%
% Note that my spm_defaults (saved as *spm_my_defaults*) implicit threshold
% is set to -inf - so we always need an explicit mask specified.
%
% Inputs:
% s = index of subject number (NOT string of subject ID) -- this function
%     is designed to run as a job array (calling the function per subject
%     in parallel)
% TR = tr of functional data (seconds) % NOTE - only one TR can be
%      specified. If TR differs between tasks, adapt this script to specify
%      which task to run.
% task = define which task runs to denoise (e.g. "rest" for task-rest runs)
% derivDir   = parent directory (e.g. ~/data/derivs) of fmriprep folder
% mask = file for explicit (whole brain) mask - this could be a template or
%        subject-specific file
%
%
% Rose Cooper
% last updated: March 2020
%
% see conn_batch for further detals of functional connectivity analysis
% options - e.g. Analysis.weight and Analysis.modulation and
% condition.param for specifying gPPI and modulators.
%------------------------------------------------------------------------%%


%where is preprocessed data
fmriprepDir = fullfile(derivDir, 'fmriprep');
%where are my denoising regressors
covDir   = fullfile(derivDir, 'confounds');
%where to put denoised data:
analysisDir  = fullfile(derivDir, 'denoised');
if ~exist(analysisDir,'dir'), mkdir(analysisDir); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear batch

%% 1. SETUP SUBJECTS and Parallel Processing

batch.Setup.isnew     = 1;    %0 if want to update existing project
batch.Setup.RT        = TR;
batch.Setup.done      = 1;    %0 to specify files files but not run step
batch.Setup.overwrite = 1;    %0 if updating existing project and want to keep previous files

fmrisubs = struct2cell(dir(fmriprepDir));
fmrisubs = fmrisubs(:,contains(fmrisubs(1,:),'sub'));
fmrisubs = fmrisubs(1,cell2mat(fmrisubs(5,:)) == 1); %directories

curSubj = fmrisubs(s);

batch.Setup.nsubjects = length(curSubj);
batch.filename=fullfile(analysisDir,curSubj{1},['task-' task '.mat']); % New conn_*.mat project name


% ---------------------------------------------------------------------- %
% %%% Alternative: parallel processing based on n subjects to put in a
% single group project
% batch.parallel.N = batch.Setup.nsubjects;
% batch.parallel.profile = 'PBS/Torque computer cluster'; %(check for your server)
% % configure this as you would a .pbs file
% batch.parallel.cmd_submitoptions = '-l nodes=1:ppn=2,walltime=0:30:00,mem=12g -m abe -M cooperrn@bc.edu';
% ---------------------------------------------------------------------- %


%% 2. Get task run numbers

files = dir(fullfile(fmriprepDir, curSubj{1}, 'func'));
files = struct2cell(files);
files = files(1,contains(files(1,:),'.tsv'));
task_files = files(contains(files,task));

%get run numbers by task:
run_IDs = {};
for nses = 1:length(task_files)
    parts = strsplit(task_files{nses},'_');
    idx = contains(parts,'run-');
    run_IDs{nses} = parts{idx};
end
nsessions = length(run_IDs);


%% 3. GET FUNCTIONALS by task and run

fprintf('\nGetting functional scans...\n');

batch.Setup.functionals = repmat({{}},[1,1]); % initialize main functional volumes for each subject/session

% grab all niis in preprocessed functional folder, as only valid runs were copied there
% a) main 4D functionals :
for nses = 1:length(run_IDs)
    for nsub = 1:length(curSubj)
        func_regexp = ['^' curSubj{nsub} '_task-' task '_' run_IDs{nses} '.*MNI.*preproc.*.nii.gz']; %fmriprep output (unsmoothed)
        scanRuns = cellstr(spm_select('FPList', fullfile(fmriprepDir, curSubj{nsub}, 'func'), func_regexp));
        if isempty(scanRuns{1})
            error('no functional files found.');
        end
        batch.Setup.functionals{nsub}{nses} = scanRuns{1};
    end
end


%% 4. GET STRUCTURALS

fprintf('\nGetting structural scans...\n');

batch.Setup.structurals = repmat({{}},[1,1]); %inirialize variable for structural scan names

% find MNI space T1 preprocessed image from fmriprep
for nsub = 1:length(curSubj)
    % T1
    strc_regexp = ['^' curSubj{nsub} '.*MNI.*preproc_T1w.nii.gz'];
    strcFile = cellstr(spm_select('FPList', fullfile(fmriprepDir, curSubj{nsub}, 'anat'), strc_regexp));
    if isempty(strcFile{1})
        error('no structural file found.');
    end
    batch.Setup.structurals{nsub} = strcFile{1};
    
    
    % GM/WM/CSF files
    classFiles = {};
    for m = 1:3
        if m == 1, class = 'CSF'; elseif m == 2, class = 'GM'; elseif m == 3; class = 'WM'; end
        mask_regexp = ['^' curSubj{nsub} '.*MNI.*' class '_probseg.nii.gz'];
        maskFile   = cellstr(spm_select('FPList', fullfile(fmriprepDir, curSubj{nsub}, 'anat'), mask_regexp));
        if isempty(maskFile{1})
            error('no file found for tissue type.');
        end
        classFiles{m} = maskFile;
    end
    
    % add to batch
    batch.Setup.masks.CSF.files{nsub}   = classFiles{1};
    batch.Setup.masks.Grey.files{nsub}  = classFiles{2};
    batch.Setup.masks.White.files{nsub} = classFiles{3};
    % specify that we only want 1 dimension
    batch.Setup.masks.CSF.dimensions   = 1;
    batch.Setup.masks.Grey.dimensions  = 1;
    batch.Setup.masks.White.dimensions = 1;
end


%% 5. GET CONFOUND REGRESSORS AND CONDITIONS

fprintf('\nAdding regressors...\n\n');

batch.Setup.conditions.names = {task}; %for denoising, we're just specifying the functional task as one 'condition'
covar_names = {'motion','aCompCor','spikes'}; %see fmriprep_confound_regressors

% add onsets and durations per session for target events, as well as
% all covariates
for nses = 1:nsessions
    for nsub = 1:length(curSubj)
        % mark whole session as this task (condition)
        batch.Setup.conditions.onsets{1}{nsub}{nses} = 0;
        batch.Setup.conditions.durations{1}{nsub}{nses} = inf;
        
        for c = 1:length(covar_names)
            c_name = covar_names{c};
            % grab covariates associated with this session for this subject:
            fileexp = ['.*' task '_' run_IDs{nses} '_' c_name '.txt'];
            sessFiles = cellstr(spm_select('FPList', fullfile(covDir, curSubj{nsub}), fileexp));
            %add to batch
            if ~isempty(sessFiles{1}) %not all subjects will have outlier spikes
                batch.Setup.covariates.names{c} = c_name;
                batch.Setup.covariates.files{c}{nsub}{nses} = sessFiles{1};
            end
        end
    end
end

% add confound covariates to denoising pipelines
batch.Denoising.confounds.names = batch.Setup.covariates.names; %all (and only) my listed nuisance covariates
for c = 1:length(batch.Denoising.confounds.names)
    batch.Denoising.confounds.deriv{c}=0; %do not add derviatives (note my motion parameters files already have temporal derivs)
end


%% 6. DEFINE STEPS AND RUN

% CONN Denoising --> detrending and regress out specified regressors from functional data
batch.Denoising.filter     = [0.008,inf];   %band-pass filter (Hz) - high pass and low pass
batch.Denoising.detrending = 1;             %1: linear detrending
batch.Denoising.regbp      = 1;             %1= reg then bp filter - best if you have spikes. 2= Simultaneous regression and BP filtering

% Additional setup steps
batch.Setup.analyses        = 3; %just vox-to-vox, to visualize denoising
batch.Setup.voxelmask       = 1; %1.Explicit mask -- note by default this would be SPM MNI brainmask.nii
% but you could specify the fmriprep template used from templateflow:
% e.g. mask = ['/path/to/tpl-MNI152NLin6Asym_res-02_desc-brain_mask.nii.gz/']
batch.Setup.voxelmaskfile   = mask;
batch.Setup.voxelresolution = 3; %same as functional volumes (not this is a code for conn, not the voxel size)
batch.Setup.analysisunits   = 2; %2 = raw, 1 = PSC
batch.Setup.outputfiles     = [0,1]; %(2) = create nifti for confound-corrected time series
batch.Setup.rois.names      = {'BrainMask'}; %note that this is just a filler given that we aren't running analyses
batch.Setup.rois.files      = {batch.Setup.voxelmaskfile};


batch.Denoising.done      = 1;
batch.Denoising.overwrite = 1;
batch.Analysis.done       = 0; % no analyses
batch.vvAnalysis.done     = 0; % no analyses
batch.dynAnalysis.done    = 0; % no analyses

%%%%%%%%%%%%%%%%%%
conn_batch(batch);
%%%%%%%%%%%%%%%%%%


%% 7. Now move and rename the doinoised .nii files for ease of use
% Also converts zeros outside of brain mask to NaNs

for nsub = 1:batch.Setup.nsubjects
    dataFile = fullfile(analysisDir,curSubj{nsub},['task-' task],'results','preprocessing',['niftiDATA_Subject' num2str(nsub,'%03.f') '_Condition000.nii']);
    newFile  = fullfile(analysisDir,curSubj{nsub},[curSubj{nsub} '_task-' task '_MNI_preproc_denoised.nii']);
    
    % move file:
    unix(['mv ' dataFile ' ' newFile]);
    
    %convert file, zeros to NaNs (using cosmo - convenient for 4D data):
    data_new = cosmo_fmri_dataset(newFile);
    data_new.samples(data_new.samples == 0) = NaN;
    cosmo_map2fmri(data_new, newFile); %overwrite
end


%% 8. Option to remove unzipped niftis from fmriprep

% fmriprep only saves .nii as zipped (.gz) for space.
% conn can accept .gz files, but it will automatically unzip them to work
% with. Here, we look for all of these unzipped niftis in the fmriprep sub-
% directories and delete them to save space:

for nsub = 1:batch.Setup.nsubjects
    unix(['rm -rf ' fullfile(fmriprepDir, curSubj) '/*/*.nii']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end