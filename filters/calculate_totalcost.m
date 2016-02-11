
function [cost] = calculate_totalcost(in,y,atlas,param)

IND = param.VoxelIdx;
SIZE = param.Dimension;

in_vol = zeros(SIZE);

for j = 1:length(IND(:,1));
    in_vol(IND(j,1),IND(j,2),IND(j,3),:) = in(:,j);
end

norml21 = find_l21_norm(in_vol,atlas,param);
norml1 = find_l1_norm(in,param);

cost = norml21+norml1+sum((in(:)-y(:)).^2)/2;
end

function norm_out = find_l1_norm(in,param)
norm_out = 0;
n = param.f_Analyze.num;
d = param.f_Analyze.den;

for ind_tem = 1:param.NbrVoxels,
    norm_out = norm_out + param.LambdaTempFin(ind_tem)*sum(abs(filter_boundary(n,d,in(:,ind_tem),'normal')));
end

end



function norm_out = find_l21_norm(in,atlas,param)
der_temp = zeros(param.Dimension(1:3));
norm_out = 0;

h = zeros(3,3,3);
h(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
h(:,:,2) = [0 1 0; 1 -6 1; 0 1 0];
h(:,:,3) = h(:,:,1); % laplace filter

H = fftn(h,param.Dimension(1:3));
for ind_vol=1:param.Dimension(4)
    der_temp = ifftn(H.*fftn(in(:,:,:,ind_vol),param.Dimension(1:3)));
    for a=1:max(atlas(:)); % number of regions
        IND_r = find(atlas == a);
        norm_out = norm_out + norm(der_temp(IND_r),2);
    end
end
end





