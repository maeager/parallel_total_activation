function [x] = strsparsity_nii(y,H,lambda,maxeig,atlas,Nit)
%#codegen

z = zeros(size(y));
z_in = zeros(size(y));

SIZE = [size(y,1),size(y,2),size(y,3)];
%H = fftn(h,size(y));

for k=1:Nit-1
    z_in = z + 1/(lambda*maxeig)*ifftn(H.*fftn(y,SIZE)) - ...
    ifftn(H.*conj(H).*fftn(z,SIZE))/maxeig;
    %k=k+1;
    
    %function out = clip(in,atlas)
    % find indices for each region
    z = z.*0;
    
    for a=1:max(atlas(:)); % number of regions
        IND_r = find(atlas == a);
        %     out(IND_r) = min(norm(in(IND_r),2),1).*in(IND_r)./norm(in(IND_r),2); % may have NAN problems...
        norm_z = norm(z_in(IND_r),2);
        if (norm_z > 1)
            z(IND_r) =  abs(z_in(IND_r))./norm_z;
            %         out(IND_r) =  sign(in(IND_r))./sqrt(length(IND_r));
        else
            z(IND_r) =  abs(z_in(IND_r));
        end
    end
end
x = y - lambda*ifftn(conj(H).*fftn(z,size(y)));
