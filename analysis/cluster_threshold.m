% Threshold a whole-brain nifti file by voxel cluster extent.

% Rose Cooper Nov 2019

% inputs:
% data = 3D data to threshold
% vox = voxel cluster extent (minimum)
% val = minimum statistical value to keep (e.g. z or t) if not thresholded
%       already
% type = type of connected components classified as a 'cluster':
%       6  = adjoining faces
%       18 = adjoining edges
%       26 = adjoining corners
%
% outputs a 3D matrix (new_data) with data cluster-thresholded
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_data] = cluster_threshold(data, vox, type, val)

% 1. value-based threshold
data(data < val) = 0;

% 2. extent-based threshold
% logical mask to get cluster sizes
clusters = bwconncomp(logical(data),type);
extent = cellfun('length', clusters.PixelIdxList);

% mask by cluster threshold
extent(extent < vox) = 0;
clusters.PixelIdxList = clusters.PixelIdxList(logical(extent));

% create new data file with only clusters > X voxels
new_data = zeros(size(data));
idx = vertcat(clusters.PixelIdxList{:})';
new_data(idx) = data(idx);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
