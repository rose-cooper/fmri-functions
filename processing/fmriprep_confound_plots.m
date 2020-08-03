function fmriprep_confound_plots(fd_th,fd_prop,fd_spike,dvars_spike,derivDir)
%% -------------------------------------------------------------------- %%
% Fmriprep's carpet plot does not show the 6 motion parameters, and it is
% also generated only in html reports with the MNI2009c template.
% This function uses the fmriprep-generated .tsv confounds to generate a 
% multi-panel confound plot per subject/run.
%
% Customize thresholds for motion, based on run inclusion criteria and
% spikes.
%
% Inputs:
%   fd_th   = threshold (mm) for judging high motion timepoints 
%       (e.g. 0.2mm)
%   fd_prop = proportion threshold for the number of high-motion time
%       points that are permitted per run before exclusion 
%       (e.g. 0.2 (20%))
%   fd_spike = threshold (mm) for marking a time point as a spike 
%       (e.g. 0.5mm) 
%   dvars_spike = threshold (dvars) for marking a time point as a spike 
%       (e.g. 2)
%   derivDir = fmriprep parent 'derivs' directory
%
% Saves in a single multi-panel plot:
% - FD, with marker at fd_th and fd_spike, and fd_prop annotation
% - std dvars, with a marker at dvars_spike
% - xyz translations
% - xyz rotations
% - global signal, acompcor, wm, csf --> plots all zeros if found NAs
% (which indicates a problem generating a mask during preprocessing)
%
%
% Rose Cooper
% last updated: March 2020
%------------------------------------------------------------------------%%


fmriprepDir = fullfile(derivDir, 'fmriprep');
imageDir    = fullfile(derivDir, 'fmriprep_plots');
if ~exist(imageDir,'dir'), mkdir(imageDir); end


fmrisubs = struct2cell(dir(fmriprepDir));
fmrisubs = fmrisubs(:,contains(fmrisubs(1,:),'sub'));
fmrisubs = fmrisubs(1,cell2mat(fmrisubs(5,:)) == 1); %directories

fprintf('\nNumber of preprocessed subjects = %d\n',length(fmrisubs));
fprintf('\nGenerating plots...\n');


for s = 1:length(fmrisubs)
    
    subID = fmrisubs{s};
    fprintf('\tSubject: %s\n',subID);
    
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
        
        if ~exist(fullfile(imageDir, [subID '_' curTask '.png']),'file')
            
            fprintf('\t\t%s\n',curTask);
            confounds = tdfread(fullfile(fmriprepDir, subID, 'func', files{t}));
            
            % FD
            fd = cellstr(confounds.framewise_displacement);
            fd{1} = '0';
            c.fd = cellfun(@str2num, fd);
            meanfd = round(mean(c.fd(2:end)),3);
            percfd = round((sum(c.fd(2:end) > fd_th)/length(c.fd(2:end)))*100,2);
            outliersfd = sum(c.fd(2:end) > fd_spike);
            
            % STD DVARS
            sdvars = cellstr(confounds.std_dvars);
            sdvars{1} = '0';
            c.dvars = cellfun(@str2num, sdvars);
            outliersdvars = sum(c.dvars(2:end) > dvars_spike);
            
            % xyx
            c.xyz = [confounds.trans_x,...
                confounds.trans_y,...
                confounds.trans_z];
            
            % rotation
            c.rot = [confounds.rot_x,...
                confounds.rot_y,...
                confounds.rot_z];
            
            % global signal (average within brainmask)
            c.global = confounds.global_signal;
            if ischar(c.global)
                c.global = zeros(1,length(c.fd));
            end
            
            % csf
            c.csf = confounds.csf;
            if ischar(c.csf)
                c.csf = zeros(1,length(c.fd));
            end
            
            % WM
            c.wm = confounds.white_matter;
            if ischar(c.wm)
                c.wm = zeros(1,length(c.fd));
            end
            
            % acompcor (first principal component)
            c.acomp = confounds.a_comp_cor_00;
            
            
            %%%%%%%%%%%%%%%% make plot %%%%%%%%%%%%%%%%%%
            x = (1:length(c.fd))';
            xmax = round(length(c.fd)*1.02);
            h = figure('Visible','off');
            
            % FD ----------------%
            subplot(4,2,1);
            plot(x,c.fd,'LineWidth',1.5, 'color',[0.2 0.6 1])
            hold on;
            plot(x,(fd_th*ones(1,length(c.fd)))',...
                'LineWidth',1, 'color',[0.4 0.4 0.4], 'LineStyle','--');
            plot(x,(fd_spike*ones(1,length(c.fd)))',...
                'LineWidth',1, 'color',[0 0 0], 'LineStyle','--');
            title('Framewise Displacement');
            xlim([0 xmax]);
            texty = max(0.5,max(c.fd)*(1-fd_prop));
            message = text(5,texty,['Mean = ',num2str(meanfd),...
                ', % > ',num2str(fd_th),' = ',num2str(percfd),...
                ', n > ',num2str(fd_spike),' = ',num2str(outliersfd)]);
            set(message,'FontSize',9);
            
            % DVARS --------------%
            subplot(4,2,3);
            plot(x,c.dvars,'LineWidth',1.5, 'color',[1 0.6 0.2]);
            hold on;
            plot(x,(dvars_spike*ones(1,length(c.fd)))',...
                'LineWidth',1, 'color',[0 0 0], 'LineStyle','--');
            title('Standardized DVARS');
            xlim([0 xmax]);
            message = text(5,1.6,['n > ',num2str(dvars_spike),' = ',num2str(outliersdvars)]);
            set(message,'FontSize',9);
            
            % Motion Translations -----%
            subplot(4,2,5);
            plot(x,c.xyz(:,1),'LineWidth',1.5, 'color',[0.1 0.1 0.1]);
            hold on;
            plot(x,c.xyz(:,2),'LineWidth',1.5, 'color',[0.4 0.4 0.4]);
            plot(x,c.xyz(:,3),'LineWidth',1.5, 'color',[0.7 0.7 0.7]);
            title('Translations');
            xlim([0 xmax]);
            
            % Motion Rotations --------%
            subplot(4,2,7);
            plot(x,c.rot(:,1),'LineWidth',1.5, 'color',[0.1 0.1 0.1]);
            hold on;
            plot(x,c.rot(:,2),'LineWidth',1.5, 'color',[0.4 0.4 0.4]);
            plot(x,c.rot(:,3),'LineWidth',1.5, 'color',[0.7 0.7 0.7]);
            title('Rotations');
            legend(["x","y","z"],'Location','best');
            xlim([0 xmax]);
            
            % Global signal -----------%
            subplot(4,2,2);
            plot(x,c.global,'LineWidth',1.5, 'color',[0.8 0 0.8]);
            title('Global Signal (brain mask)');
            xlim([0 xmax]);
            
            % a comp cor PC ----------%
            subplot(4,2,4);
            plot(x,c.acomp,'LineWidth',1.5, 'color',[0 0.8 0.6]);
            title('aCompCor PC 00');
            xlim([0 xmax]);
            
            % white matter signal ----------%
            subplot(4,2,6);
            plot(x,c.wm,'LineWidth',1.5, 'color',[0.2 0.2 0.2]);
            title('White Matter Signal');
            xlim([0 xmax]);
            
            % CSF signal ----------%
            subplot(4,2,8);
            plot(x,c.csf,'LineWidth',1.5, 'color',[0.6 0.6 0.6]);
            title('CSF Signal');
            xlim([0 xmax]);
            
            
            %%%%%%%%%%% save plot %%%%%%%%%%%%
            set(h,'position',[0 0 1250 750]);
            print([imageDir subID '_' curTask], '-dpng', '-r400');
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%