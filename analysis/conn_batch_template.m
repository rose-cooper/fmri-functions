%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Template script to run first level connectivity analyses in CONN, from
% fmriprep preprocessed data

% (https://sites.google.com/view/conn/)

% *** 'help conn_batch' for CONN batch code options ***

% Rose Cooper - last updated July 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up file paths

clearvars; clc;

TR  = 1.5; %TR for all functional scans

analysis_name = ['my-analysis-name'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create output directory for connectivity data
outDir = ['/path/to/my/output/'];
if ~exist(outDir,'dir')
    mkdir(outDir);
end

% SPM & CONN toolboxes
spmDir = ['/path/to/toolboxes/'];
addpath(genpath(spmDir));

% structural/functional data?
fmriDir  = ['/data/derivs/fmriprep/'];
% ROI files (seeds)?
ROIdir  = ['/path/to/masks/'];
% first level covariates (e.g. confounds)
covDir  = ['/path/to/txtfiles/'];
% first level condition onset/duration files (if any)
taskDir = ['/path/to/task.matfiles/'];


% define ROI names
myROIs = {'ROIa','ROIb','ROIc'};

% name of task conditions


clear batch
%% 1. SETUP SUBJECTS

batch.filename=fullfile(outDir,[analysis_name '.mat']); % New conn_*.mat project name

%get all subjects who have first level covariate files (text files) for CONN:
subjects = struct2cell(dir(covDir));
subjs    = subjects(1,contains(subjects(1,:),'sub'));
NSUBS    = size(subjs,1);

batch.Setup.isnew     = 1; %0 to update existing project
batch.Setup.nsubjects = NSUBS;
batch.Setup.RT        = TR; %TR (seconds)



%% 2. GET FUNCTIONALS

fprintf('\nAdding functional scans...\n');

batch.Setup.functionals = repmat({{}},[NSUBS,1]); %main functional volumes for each subject/session
%batch.Setup.roifunctionals.roiextract_functionals = repmat({{}},[NSUBS,1]); %additional (e.g. unsmoothed) functional volumes for ROI analyses

for nsub=1:NSUBS
    fprintf('...%s\n',subjs{nsub});
    
    % 4D functionals :
    func_regexp = ['\<' subjs{nsub} '.*task-MyTaskName.*MNI.*_preproc.nii'];
    clear funcFiles
    % find nii file for every scan run
    funcFiles = cellstr(spm_select('FPList', [fmriDir subjs{nsub} '/func/'], func_regexp)); 
    if isempty(funcFiles)
        error('no functional files found!');
    end
    nsessions = length(funcFiles);
    
    for nses=1:nsessions
        batch.Setup.functionals{nsub}{nses} = funcFiles{nses};
    end
    
    batch.Setup.roifunctionals.roiextract = 1; % 1 to extract roi data from main functional files % 4 indicate user defined data files
end



%% 3. GET STRUCTURALS

fprintf('\nAdding structural scans...\n');

batch.Setup.structurals = repmat({{}},[NSUBS,1]);

for nsub=1:NSUBS
    fprintf('...%s\n',subjs{nsub});
    
    % find MNI space preprocessed image from fmriprep
    strc_regexp = ['\<' subjs{nsub} '.*MNI.*_preproc.nii'];
    strcFile    = cellstr(spm_select('FPList', [fmriDir subjs{nsub} '/anat/'], strc_regexp));
    if isempty(strcFile)
        error('no structural file found!');
    end
    
    batch.Setup.structurals{nsub} = strcFile{1};
    
    %get fmriprep grey/WM/csf masks:
    targetFiles = {'CSF','GM','WM'};
    for m = 1:length(targetFiles)
        mask_regexp = ['\<' subjs{nsub} '.*MNI.*' targetFiles{m} '_probtissue.nii'];
        maskFile   = cellstr(spm_select('FPList', [fmriDir subjs{nsub} '/anat/'], mask_regexp));
        if isempty(maskFile)
            error(['no ' targetFiles{m} ' file found!']);
        end
        
        % add to batch
        if strcmp(targetFiles{m},'CSF')
             batch.Setup.masks.CSF{nsub} = maskFile{1};
        elseif strcmp(targetFiles{m},'GM')
             batch.Setup.masks.Grey{nsub} = maskFile{1};
        elseif strcmp(targetFiles{m},'WM')
             batch.Setup.masks.White{nsub} = maskFile{1};
        end
    end
end
% specify that we only want 1 dimension
batch.Setup.masks.CSF.dimensions   = 1;
batch.Setup.masks.Grey.dimensions  = 1;
batch.Setup.masks.White.dimensions = 1;



%% 4. GET ROIs

fprintf('\nAdding ROI files...\n');

% note. subject-specific ROIs can be specified (e.g. in native space) by
% using batch.Setup.rois.files{nroi}{nsub} = ['insert subject-specific ROI file name']
for nroi = 1:length(myROIs)
    batch.Setup.rois.names{nroi} = [myROIs{nroi}];
    batch.Setup.rois.files{nroi} = [ROIdir myROIs{nroi} '_ROI.nii'];
    batch.Setup.rois.dimensions{nroi} = 1; %extract single component (mean across voxels)
    batch.Setup.rois.roiextract(nroi) = 0; %1 indicates to extract from additionally specified functional data, or same as main functionals (0).
end



%% 5. GET COVARIATES *by run*

%%% this section adds in nuisance confounds and task conditions. You can
%%% also add in parametric modulators as first level covariates (same way
%%% as confounds) to modulate task condition effects. 
%%% Any parametric modulators would need to be convolved with the HRF to be
%%% a vector with length nTRs

fprintf('\nAdding covariates...\n');

for nsub=1:NSUBS
    fprintf('...%s\n',subjs{nsub});
    
    %get number of scan runs:
    nsessions = length(batch.Setup.functionals{nsub});
    
    for nses = 1:nsessions
        
        % grab all covariates associated with this run for this subject (might just be 1 text file with all confounds):
        sessFiles = cellstr(spm_select('FPList', [covDir subjs{nsub} '/'], ['.*_run0' num2str(nses) '.txt']));
        
        % add names and covariate files to conn_batch
        regNames = {}; %store all variables to include for denoising later
        for ncov = 1:length(sessFiles)
            %get regressor name
            curReg = strsplit(sessFiles{ncov},'/');
            curReg = curReg{end}; %grab final part, which is the regressor name
            curReg = strsplit(curReg,'_'); %now remove the runID from name
            curReg = curReg{1};
            regNames{ncov} = curReg;
            
            %add to batch
            batch.Setup.covariates.names{ncov} = curReg; % e.g. 'confounds'
            batch.Setup.covariates.files{ncov}{nsub}{nses} = sessFiles{ncov};
        end
        
        
        % now add any task conditions that you have (if any):
        % grab all onset files for this subject and run:
        sessFiles = cellstr(spm_select('FPList', [taskDir subjs{nsub} '/'], ['.*_run0' num2str(nses) '.mat']));
        
        for ntask = 1:length(sessFiles)
            %get condition name
            curReg = strsplit(sessFiles{ntask},'/');
            curReg = curReg{end}; %grab final part, which is the regressor name
            curReg = strsplit(curReg,'_'); %now remove the runID from name
            
            %get task onsets and durations from .mat file
            task_info = load(sessFiles{ntask});
            
            batch.Setup.conditions.names{ntask} = curReg{1};
            batch.Setup.conditions.onsets{ntask}{nsub}{nses} = task_info.onsets;       %can also be 0 if it's the whole run (e.g. rest)
            batch.Setup.conditions.durations{ntask}{nsub}{nses} = task_info.durations; %use inf for the whole functional run
            batch.Setup.conditions.param{ntask}{nsub}{nses} = 0; %or index number of first level covariate for temporal parametric modulation
            % note, if using parametric modulation, you should attribute it
            % to it's own "condition" with ons=0 and duration=inf
        end
    end %end of loop through sessions
end



%% 6. DEFINE ANALYSIS steps

% CONN Setup
batch.Setup.analyses=[1,2]; %ROI-to-ROI(1) and seed-to-voxel(2)
batch.Setup.analysisunits= 2; %2 = raw, 1 = PSC
batch.Setup.voxelmask=1; % 1. Explicit mask (brainmask.nii - default), 2. implicit subject specific
batch.Setup.voxelmaskfile = ['/templates/my-MNI-template.nii'];
batch.Setup.voxelresolution=3; %same as functional volumes (see options)
batch.Setup.outputfiles = [0,0,1,0,1,0]; %saves seed-to-voxel r maps and fdr-p maps, see options
    
% CONN Denoising --> detrending and regress out specified regressors from functional data
batch.Denoising.filter=[0.01,0.1]; %band-pass filter (Hz) - resting state default, or [0.008 inf] to match spm high pass filtering for tasks
batch.Denoising.detrending=1;       %1: linear detrending
batch.Denoising.confounds.names=regNames; %include all of my manually-added confounds (and not CONN's WM/GM/CSF)
for r = 1:length(regNames)
    batch.Denoising.confounds.deriv{r}=0; %do not add derviatives
end

% CONN First level
batch.Analysis.analysis_number='task_gppi'; %numerical index or custom string for analysis ID
batch.Analysis.measure=3;     % 1=bivariate correlation, 3=bivariate regression, 4=multivariate regression
batch.Analysis.weight=2;      % 1.no within-condition weighting, 2. HRF (if estimating connectivity for event-related design)
batch.Analysis.modulation=1;  % 0 = standard weighted, 1 = gPPI of condition-specific temporal modulation factor.
batch.Analysis.type=3; %1='ROI-to-ROI', 2=Seed-to-voxel, 3=all
batch.Analysis.sources=batch.Setup.rois.names; % define seeds (all source rois)



%% Run first level analyses

fprintf('\n ----- Running Setup, Denoising, and First level analyses -----\n');
batch.Setup.done          = 1; % 0 to specify files files but not run step
batch.Setup.overwrite     = 1; % 0 if updating existing project and want to keep previous files
batch.Denoising.done      = 1;
batch.Denoising.overwrite = 1;
batch.Analysis.done       = 1;
batch.Analysis.overwrite  = 1;

conn_batch(batch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%