function [h_matrices] = roi_homogeneity(func, rois)
%% -------------------------------------------------------------------- %%
% This function uses spm and cosmomvpa tools to measure roi homogeniety -
% defined as within-roi connectivity (similarity in voxel time series)
%
% Inputs:
% func = cell array of 4D functional file names (e.g. one per subject/per run)
% rois = cell array of roi file names to analyse
%
% Outputs:
% h_matrices = {roi}{func}[voxel x voxel] array of within-roi connectivity
% values (fisher z transformed)
%
% Rose Cooper, March 2020
%------------------------------------------------------------------------%%

cosmo_warning('off');

fprintf('Calculating ROI homogeneity ...\n\n');

h_matrices = {};
for f = 1:length(func)
    func_file = func{f};
    fprintf('\tFunctional file: %s\n',func_file);
    
    % use cosmo to grab data: -------------------------------------
    data = cosmo_fmri_dataset(func_file);
    % in case missing values are 0 not NaN
    data.samples(data.samples == 0) = NaN;
    fprintf('\t\tdimensions: %d time points x %d voxels\n\n',size(data.samples,1),size(data.samples,2));
    
    for r = 1:length(rois)
        roi_file = rois{r};
        fprintf('\t\tROI file: %s\n',roi_file);
        
        roi_mask = logical(spm_read_vols(spm_vol(roi_file)));
        
        %check dimensions:
        if length(roi_mask(:)) ~= size(data.samples,2)
            error('ROI and functional files are not in the same space');
        end
        fprintf('\t\tcontains %d voxels\n',sum(roi_mask(:)));
        
        % mask func with roi
        roi_data = data.samples(:,roi_mask);
        if sum(isnan(roi_data(:))) > 0 %any missing voxels in mask?
            error('NaNs present in ROI data');
        end
        
        % correlate voxels (columns)
        vox_cor = atanh(corr(roi_data, 'type', 'Spearman'));
        %mask out diagonal
        diag_mask = logical(diag(ones(size(vox_cor,1),1)));
        vox_cor(diag_mask) = NaN;
        
        fprintf('\t\thomogeneity = %s\n\n', num2str(round(mean(vox_cor(:),'omitnan'),3)));
        
        % add to output:
        h_matrices{r}{f} = vox_cor;
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%