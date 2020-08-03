function threshold_conn(name, binary, type, weight)
% threshold connectivity map according to various options:

% name = string of .nii map
% binary = "binarized" or "weighted" to save out connections as binarized or not
% type = type of threshold to use ("quantile" = quantile of connections, e.g. top 10% (.9),
% or "r" = absolute weight, e.g. r > .2). 
% weight = quantile (0-1) or r value (0-1).

% saves out thresholded .nii map

% -------------------------------------------------------------------- %
v = spm_vol([name '.nii']);
conn = spm_read_vols(v);

% quantile or r?
if strcmp(type,'quantile')
    th = quantile(conn(~isnan(conn)),weight);
elseif strcmp(type,'r')
    th = weight;   
end

% threshold
conn(conn < th) = 0;

% binarise?
if strcmp(binary,'binarized')
    conn(conn > 0) = 1;
end

% save
v.fname = char(strcat(name,'_',binary,'_',type,'_',num2str(weight),'.nii'));
spm_write_vol(v,conn);

end
