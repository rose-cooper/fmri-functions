function fmriprep_confound_regressors(fd_spike,dvars_spike,derivDir)
%% -------------------------------------------------------------------- %%
% Generate confound regressors from fmriprep output
%
% Loads tsv files geneated from fmriprep to generate nuisance
% regressors per scan run/subject for data denoising
%
% saves 3 text files, with one column per confound:
% motion    --> FD, 6 realignment params and temporal derivatives (13)
% aCompCor  --> 6 PCs of a combined WM and CSF mask
% spikes    --> one regressors per outlying TR to effectively censor.
%               Criteria = FD > 0.6mm and.or STD DVARS > 2.
%
% Inputs:
%   fd_spike = threshold (mm) for marking a time point as a spike
%       (e.g. 0.5mm)
%   dvars_spike = threshold (dvars) for marking a time point as a spike
%       (e.g. 2)
%   derivDir = fmriprep parent 'derivs' directory
%
% Rose Cooper
% last updated: March 2020
%------------------------------------------------------------------------%%

%where is regressor information (.tsv files from fmriprep)
fmriprepDir = fullfile(derivDir, 'fmriprep');
%where to save confound text files
saveDir     = fullfile(derivDir, 'confounds');
if ~exist(saveDir,'dir'), mkdir(saveDir); end


fmrisubs = struct2cell(dir(fmriprepDir));
fmrisubs = fmrisubs(:,contains(fmrisubs(1,:),'sub'));
fmrisubs = fmrisubs(1,cell2mat(fmrisubs(5,:)) == 1); %directories

fprintf('\nNumber of preprocessed subjects = %d\n',length(fmrisubs));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EDIT these variable names as needed -----------------------------------
target_motion = {'framewise_displacement',...
    'trans_x','trans_x_derivative1',...
    'trans_y','trans_y_derivative1',...
    'trans_z','trans_z_derivative1',...
    'rot_x','rot_x_derivative1',...
    'rot_y','rot_y_derivative1',...
    'rot_z','rot_z_derivative1',...
    };
target_aComp = {'a_comp_cor_00','a_comp_cor_01','a_comp_cor_02',...
    'a_comp_cor_03','a_comp_cor_04','a_comp_cor_05',...
    };
% -----------------------------------------------------------------------

fprintf('\nCreating confound regressors...\n');
fprintf('\tFD spike threshold = %s mm\n',num2str(fd_spike));
fprintf('\tDVARS spike threshold = %d\n\n',dvars_spike);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:length(fmrisubs)
    
    subID = fmrisubs{s};
    fprintf('\tSubject: %s',subID);
    
    %make subject folder:
    subjDir   = fullfile(saveDir, subID);
    if ~exist(subjDir,'dir'), mkdir(subjDir); end
    
    %get .tsv files:
    files = struct2cell(dir(fullfile(fmriprepDir, subID, 'func')));
    files = files(1,contains(files(1,:),'.tsv'));
    
    for t = 1:length(files)
        
        %get task name and run ID (if applicable)
        taskname = strsplit(files{t},'_');
        curTask = taskname{contains(taskname,'task')};
        if ~isempty(contains(taskname,'run'))
            runID = taskname{contains(taskname,'run')};
            curTask = strcat(curTask,'_',runID);
        end
        
        fprintf('\n\t\t\t%s',curTask);
        confounds = tdfread(fullfile(fmriprepDir, subID, 'func', files{t}));
        confound_names = fieldnames(confounds);
        
        % convert to matrix, replacing first n/a with zero
        idx = structfun(@ischar, confounds);
        for c = 1:length(idx)
            if idx(c)
                temp = cellstr(confounds.(confound_names{c}));
                if strcmp(temp{1},'n/a')
                    temp{1} = '0';
                end
                confounds.(confound_names{c}) = cellfun(@str2num, temp);
            end
        end
        confounds = struct2array(confounds);
        
        
        %find target regressor columns for motion: ------------------------ %
        motion_Cols = cellfun(@(x) find(endsWith(confound_names,x)),target_motion);
        %get target confound data:
        R = confounds(:,motion_Cols);
        fileName = fullfile(subjDir, [subID '_' curTask '_motion.txt']);
        writetable(array2table(R),fileName,'Delimiter',' ','WriteVariableNames',false);
        
        
        %find target regressor columns for aCompCor: ---------------------- %
        aComp_Cols = cellfun(@(x) find(endsWith(confound_names,x)),target_aComp);
        %get target confound data:
        A = confounds(:,aComp_Cols);
        fileName = fullfile(subjDir, [subID '_' curTask '_aCompCor.txt']);
        writetable(array2table(A),fileName,'Delimiter',' ','WriteVariableNames',false);
        
        
        % now create spike (scrubbing) regressors --------------------
        spike_idx = [];
        
        % check to see if we need to remove non-steady-state outliers
        nss = find(contains(confound_names,'non_steady_state'))';
        if ~isempty(nss)
            for n = nss
                spike_idx = [spike_idx find(confounds(:,n) == 1)];
            end
        end
        
        % now FD and dvars
        fd_col    = contains(confound_names,'framewise_displacement');
        dvars_col = contains(confound_names,'std_dvars');
        idx = find(confounds(:,fd_col) > fd_spike | confounds(:,dvars_col) > dvars_spike)';
        if ~isempty(idx)
            spike_idx = [spike_idx idx];
        end
        spike_idx = unique(spike_idx); % just in case double up on nss and spikes
        
        
        % save scrubbing file is any to exclude:
        if ~isempty(spike_idx)
            spike_regs = zeros(size(R,1),length(spike_idx));
            for s = 1:length(spike_idx)
                spike_regs(spike_idx(s),s) = 1;
            end
            fprintf('\t%d spikes found',size(spike_regs,2));
            
            fileName = fullfile(subjDir, [subID '_' curTask '_spikes.txt']);
            writetable(array2table(spike_regs),fileName,'Delimiter',' ','WriteVariableNames',false);
        end
        
    end % end of loop through tasks/runs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n');
end %end of loop through subjects ---------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%