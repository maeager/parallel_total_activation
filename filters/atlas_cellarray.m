function cellIND = atlas_cellarray(atlas)
% #codegen
% find indices for each region
for a=1:max(atlas(:)); % number of regions
    cellIND{a} = int32(find(atlas == a));
end