function fmriprep_data_check(nAnat,nFunc,dataDir)
%% -------------------------------------------------------------------- %%
% This function checks the number of files in each subject's fmriprep folder,
% as well as the number of subjects run (relative to those in source data).
% It will also report the subject IDs with error log folders.
%
% Inputs:
%   nAnat = number of anat files that should be generated from
%       fmriprep (depending in your preprocessing options)
%   nFunc = number of func files that should be generated from
%       fmriprep (depending in your preprocessing options)  
%   dataDir = '/path/to/data/' (contains sourcedata and derivs/fmriprep)
%
%
% Rose Cooper
% last updated: March 2020
%------------------------------------------------------------------------%%

sourceDir   = fullfile(dataDir, 'sourcedata');
fmriprepDir = fullfile(dataDir, 'derivs', 'fmriprep');


% 1. number of source data subjects?
sourcesubs = struct2cell(dir(sourceDir));
sourcesubs = sourcesubs(1,contains(sourcesubs(1,:),'sub'));

fprintf('\nNumber of sourcedata subjects = %d\n',length(sourcesubs));


% 2. number of preprocessed subjects?
fmrisubs = struct2cell(dir(fmriprepDir));
fmrisubs = fmrisubs(:,contains(fmrisubs(1,:),'sub'));
fmrisubs = fmrisubs(1,cell2mat(fmrisubs(5,:)) == 1); %directories

fprintf('\nNumber of preprocessed subjects = %d\n',length(fmrisubs));

subs_left = setdiff(sourcesubs,fmrisubs);
fprintf('\tSubjects left to preprocess:\n');
fprintf(1,'\t\t%s\n',subs_left{:});


% 3. search through preprocessed subjects to check that they all have the
% correct number of files:

% number of anatomical files ok?
anat_files = cellfun(@(x) filenumber(x,'anat',fmriprepDir), fmrisubs);
fprintf('\nNumber of subjects with incorrect anat files: %d\n',sum(anat_files ~= nAnat));
subs_incorrect_anat = fmrisubs(anat_files ~= nAnat);
if ~isempty(subs_incorrect_anat)
    fprintf('\tSubjects with incorrect anat files:\n');
    fprintf(1,'\t\t%s\n',subs_incorrect_anat{:});
end


% number of functional files ok?
func_files = cellfun(@(x) filenumber(x,'func',fmriprepDir), fmrisubs);
fprintf('\nNumber of subjects with incorrect func files: %d\n',sum(func_files ~= nFunc));
subs_incorrect_func = fmrisubs(func_files ~= nFunc);
if ~isempty(subs_incorrect_func)
    fprintf('\tSubjects with incorrect func files:\n');
    fprintf(1,'\t\t%s\n',subs_incorrect_func{:});
end


% any error log folders in subject directories?
log_files = cellfun(@(x) logerror(x,fmriprepDir), fmrisubs);
log_files = logical(log_files);
fprintf('\nNumber of subjects with an error log folder: %d\n',sum(log_files));
subs_log = fmrisubs(log_files == 1);
if ~isempty(subs_log)
    fprintf('\tSubjects with error logs:\n');
    fprintf(1,'\t\t%s\n',subs_log{:});
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l = filenumber(x,folder, mydir)
files = struct2cell(dir(fullfile(mydir, x, folder)));
files = files(1,cell2mat(files(5,:)) == 0);
l = length(files);
end

function f = logerror(x,mydir)
f = exist(fullfile(mydir, x, 'log'),'dir');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%