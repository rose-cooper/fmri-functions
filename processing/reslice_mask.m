function reslice_mask(ref,source)
%% -------------------------------------------------------------------- %%
% This function uses spm's reslice tool to match a template to functional
% data or a target file on dimensions/vox size/orientation
%
% NOTE. the reference and to-be-resliced images MUST already be normalized
% or coregistered to the same template space. 
%
% Uses nearest neighbor interpolation as appropriate for mask/labeled
% images
%
% Inputs:
% ref = string of reference image name (single file)
% source = cell array of to-be-resliced files
%
% Rose Cooper, March 2020
%------------------------------------------------------------------------%%

clear matlabbatch
matlabbatch{1}.spm.spatial.coreg.write.ref = {ref};
matlabbatch{1}.spm.spatial.coreg.write.source = source;
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

spm_jobman('run', matlabbatch);

end