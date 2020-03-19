%% bids_format_data %% --------------------------------------------------%%
%
% This example script converts dicom images into BIDS formatted data:
% https://bids.neuroimaging.io/
%
% It also 'de-faces' a T1 image.
%
%
% Rose Cooper
% last updated: March 2020
%------------------------------------------------------------------------%%

clearvars; clc;

%% Set up directories

% tools for dicom conversion and de-facing anatomical
% download these tools at:
% https://github.com/rordenlab/dcm2niix
% https://surfer.nmr.mgh.harvard.edu/fswiki/mri_deface
addpath(genpath('/path/to/dcm2niix'));
addpath(genpath('/path/to/mri_deface'));


% Directory information
b.dicomDir        = '/path/to/data/dicoms/';     %dicoms
b.sourcedataDir   = '/path/to/data/sourcedata/'; %to put BIDS niftis


% Grab DICOM folder names (per subject):
dicom_info  = struct2cell(dir(dicom_folder));
% **customize to select by common string in folder names**
sub_folders = dicom_info(1,contains(dicom_info(1,:),'SUB'));

% Extract subject ID number folder name (customize), e.g. '001','020' and
% create name for BIDS folders using the 'sub-' prefix
ID = {}; subjects = {};
for s = 1:length(sub_folders)
    info  = strsplit(sub_folders{s},'_');
    ID(s) = info(3);
    subjects{s} = ['sub-' ID{s}];
end


%% Create dataset_description.json and participants.tsv

fid = fopen(fullfile(b.sourcedataDir, 'dataset_description.json'),'w');
jsontext = strcat('{','\n\t"Name": "My Project Name",',...
    '\n\t"BIDSVersion": "1.3.12",',...
    '\n\t"Authors": ["Rose Cooper"]',...
    '\n}');
fprintf(fid,jsontext); fclose(fid);

fid = fopen(fullfile(b.sourcedataDir, 'participants.tsv'),'w');
tsvtext = [{'participant_id'};subjects];
fprintf(fid,'%s\t\n',tsvtext{:});
fclose(fid);


%% Define scan run numbers for participants

% default dicom file scan numbers
b.taskScans = [8 11 14 20 23 26];
b.fmapScans = [17 18];
% Single-band reference scans come before every functional scan
b.taskScans_sbref = b.taskScans - 1;

% name of task for functional runs
task = 'Memory';


%% Format files by subject

for s = 1:length(subjects)
    
    fprintf('\n\nWorking on subject %s...\n',ID{s});
    
    % Current subject
    subject_dicom  = fullfile(b.dicomDir, sub_folders{s});
    
    % Move to subject data folder
    subject_BIDS   = fullfile(b.sourcedataDir, subjects{s});
    if ~exist(subject_BIDS,'dir'),mkdir(subject_BIDS); end
    
    % ALL DICOM TO (zipped) NIFTI
    fprintf('\nRunning dicom conversion with dcm2niix\n');
    fprintf(' - Files from: %s\n',subject_dicom);
    fprintf(' - Output to: %s\n',subject_BIDS);
    % -b y = generates a BIDS data file (.json)
    % -z y = compress to .nii.gz
    command = ['/path/to/dcm2niix/dcm2niix -b y -z y -o ' subject_BIDS ...
        ' -f "' subjects{s} '_%p_%s" ' subject_dicom];
    status = system(command);
    if status
        error('Error on dicom conversion');
    end
    
    
    %% Move functionals
    
    funcDir = fullfile(subject_BIDS, 'func');
    if ~exist(funcDir), mkdir(funcDir); end
    
    fprintf('\nMoving functional files to: %s\n',funcDir);
    for i=1:length(b.taskScans)
        
        runID     = sprintf('%02d',b.taskScans(i));
        runID_new = sprintf('%02d',i);
        
        fname_base = [subjects{s} '_BOLD_Minn_HCP_2mm_' runID];
        F = spm_select('List',subject_BIDS,fname_base);
        if isempty(F)
            error('No functional file found to match run ID');
        end
        
        for ff = 1:size(F,1) % get all files labeled by fname_base
            fname = F(ff,:);
            fname_new = [subjects{s} '_task-' task '_run-' runID_new '_bold'];
            Fnew = strrep(fname,fname_base,fname_new);
            unix(['mv ' subject_BIDS filesep fname ' ' funcDir filesep Fnew]);
        end
    end
    
    
    %% Move single-band references
    
    for i=1:length(b.taskScans_sbref)
        
        runID     = sprintf('%02d',b.taskScans_sbref(i));
        runID_new = sprintf('%02d',i);
        
        fname_base = [subjects{s} '_BOLD_Minn_HCP_2mm_SBREF_' runID];
        F = spm_select('List',subject_BIDS,fname_base);
        if isempty(F)
            error('No single band reference file found to match run ID');
        end
        
        for ff = 1:size(F,1) % get all files labeled by fname_base
            fname = F(ff,:);
            fname_new = [subjects{s} '_task-' task '_run-' runID_new '_sbref'];
            Fnew = strrep(fname,fname_base,fname_new);
            unix(['mv ' subject_BIDS filesep fname ' ' funcDir filesep Fnew]);
        end
    end
    
    
    %% Move anatomicals
    
    anatDir = fullfile(subject_BIDS, 'anat');
    if ~exist(anatDir), mkdir(anatDir); end
    
    fprintf('\nMoving T1 files to: %s\n',anatDir);
    F = spm_select('List',subject_BIDS,'T1_MEMPRAGE');
    if isempty(F)
        error('No anatomical scan found');
    end
    
    for ff = 1:size(F,1) % get all files labeled by fname_base
        fname = F(ff,:);
        [~,~,ext] = fileparts(fname);
        Fnew =  [subjects{s} '_T1w' ext];
        unix(['mv ' subject_BIDS filesep fname ' ' anatDir filesep Fnew]);
    end
    
    
    %% Remove any residual files (e.g., AA, localizer)
    
    unix(['rm -f ' subject_BIDS filesep subjects{s} '_*']);
    
    
    %% De-face T1 - to anonymize data
    
    P = spm_select('FPList',anatDir,'T1w.nii');
    if size(P,1) > 1
        error('Found too many T1 files. There should be only one per subject.');
    elseif isempty(P)
        error('No T1 file found');
    end
    fprintf('\nDefacing T1 image using mri_deface:\n - %s\n',P)
    
    command = ['/path/to/mri_deface/mri_deface ' P ' talairach_mixed_with_skull.gca face.gca ' P '_defaced.nii'];
    status = system(command);
    if status
        error('Error on mri deface');
    end
    
    % Replace original T1 image with defaced image:
    system(['mv -f ' P '_defaced.nii ' P]);
    % Delete copy of defaced image:
    system(['rm -f ' P '_defaced.nii']);
    
end

%% END
fprintf('\n\nBIDS formatting complete!\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%