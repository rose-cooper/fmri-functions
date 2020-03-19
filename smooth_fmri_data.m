function smooth_fmri_data(s,task,fwhm,baseDir)
%% -------------------------------------------------------------------- %%
% Smooths functional .nii file in a subject directory (this could be either
% straight from fmriprep or a denoised data file) - prefix = 'sm_'
%
% This function is set up to run subjects in parallel as a job array
%
% Inputs:
% s = subject index, not actual 'sub-XXX' ID
% task = suffix to 'task-'
% fwhm = smoothing parameter, e.g. 6 (for [6 6 6])
% baseDir = base directory to grab functional files (e.g.
% ~/data/derivs/fmriprep or ~/data/derivs/denoised)
%
% Rose Cooper
% last updated: March 2020
%------------------------------------------------------------------------%%

% subjects (to get sub ID)
fmrisubs = struct2cell(dir(baseDir));
fmrisubs = fmrisubs(:,contains(fmrisubs(1,:),'sub'));
fmrisubs = fmrisubs(1,cell2mat(fmrisubs(5,:)) == 1); %directories


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get files to unsmooth

curSubj = fmrisubs{s};
fprintf('Smoothing data for %s...\n',curSubj);

% find all functional runs matching subject and task:
zip = 0;
func_regexp = ['^' curSubj '_task-' task '.*MNI.*preproc.*.nii$'];
scanRuns = cellstr(spm_select('FPListRec', fullfile(baseDir, curSubj), func_regexp));
if isempty(scanRuns{1})
   func_regexp = ['^' curSubj '_task-' task '.*MNI.*preproc.*.nii.gz$'];
   scanRuns = cellstr(spm_select('FPListRec', fullfile(baseDir, curSubj), func_regexp));
   zip = 1; gunzip(scanRuns);
   scanRuns = cellfun(@(x) x(1:(length(x)-3)),scanRuns,'UniformOutput',0);
end


%% Run smoothing

clear matlabbatch
matlabbatch{1}.spm.spatial.smooth.data = scanRuns;
matlabbatch{1}.spm.spatial.smooth.fwhm = repmat(fwhm,1,3);
matlabbatch{1}.spm.spatial.smooth.dtype = 0; %output data type is same as input data type
matlabbatch{1}.spm.spatial.smooth.im = 1;    %preserve any masking of input image (NaNs)
matlabbatch{1}.spm.spatial.smooth.prefix = 'sm_';
spm_jobman('run', matlabbatch);


%% clear unzipped files, if applicable

if zip == 1
    for r = 1:length(scanRuns)
        unix(['rm -rf ' scanRuns{r}]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end