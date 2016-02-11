function [z] = strsparsity_par(z,y,h,lambda,maxeig,atlas)
%#codegen


z_in = abs(z + 1/(lambda*maxeig)*atlas_filter_par(y,h,0) - atlas_filter_par(z,h,1)./maxeig);

%function out = clip(in,atlas)
% find indices for each region
z = z.*0;

for a=1:max(atlas(:)); % number of regions
    IND_r = find(atlas == a);
    %     out(IND_r) = min(norm(in(IND_r),2),1).*in(IND_r)./norm(in(IND_r),2); % may have NAN problems...
    if (norm(z_in(IND_r),2) > 1)
        z(IND_r) =  (z_in(IND_r))./norm(z_in(IND_r),2);
        %         out(IND_r) =  sign(in(IND_r))./sqrt(length(IND_r));
    else
        z(IND_r) =  (z_in(IND_r));
    end
end