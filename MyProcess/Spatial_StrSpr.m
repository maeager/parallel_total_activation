function [x_out,param] = Spatial_StrSpr(y,atlas,param)

%28.11.2011 isik, same as fgp_ss.m

% Structured sparsity (denoising)

%%%%%%---INPUTS---%%%%%%
%   x - initial value of x (may be output of total activation regularization)
%   y - noisy signal (row vector)
% atlas = functional atlas
% lambda = param.LambdaSpat; %regularization parameter
% Nit = param.NitSpat; %Number of iterations
% IND = param.VoxelIdx; %[xi,yi,zi] = nonzero indices to convert 2D x to 4D x_vol
% SIZE=param.Dimension; %the size of the 4D data [v1,v2,v3,t]
% order=param.OrderSpat; %order of the derivative filter 'first', 'second', 'zero'
% dim =param.dimStrSpr; %dimension (2,3 filter in 2D or 3D)
%%%%%%%

%%%%%%---OUTPUT---%%%%%%
% x_out - result of denoising,

%%% This code is another variant of fgp algorithm (fgp_general_v2.m) for
%%% STRUCTURED SPARSITY!!!

% solves 1/2 ||y-x||^2 + lambda*||D^order x||_{s,2,1}
% created on 3.10.2011

% sum_{i=1:T} sum_{j in R_i} ||D^order{x(j)}||_2 R_i: region i of atlas. T =
% number of regions in atlas


%%%%%----TO DO LIST----%%%%%
%%% INCLUDE DIFFERENT DIMENSIONS AND DIFFERENT ORDERS OF DERIVATIVE TO
%%% COMPARE

% disp('structured sparsity')

lambda = param.LambdaSpat; %regularization parameter
Nit = param.NitSpat; %Number of iterations
IND = param.VoxelIdx; %[xi,yi,zi] = nonzero indices to convert 2D x to 4D x_vol
SIZE=param.Dimension; %the size of the 4D data [v1,v2,v3,t]
order=param.OrderSpat; %order of the derivative filter 'first', 'second', 'zero'
dim =param.dimStrSpr; %dimension (2,3 filter in 2D or 3D)

if (length(SIZE) ~= 4)
    error('SIZE should have 4 dimensions.')
end

x_out = zeros(size(y)); % out

%%%% convert to 4D
x_vol = zeros(SIZE);
y_vol = zeros(SIZE);

if(strcmpi(param.File_EXT,'mat') && (param.COST_SAVE == 1))
    temp_derivative2 = zeros(SIZE(1:3));
    temp = zeros(SIZE(1:3));
    cost = zeros(Nit,SIZE(4));
end

for j = 1:length(IND(:,1));
    y_vol(IND(j,1),IND(j,2),IND(j,3),:) = y(:,j);
end
%%%%


z = zeros(size(y_vol));

% lambda = zeros(SIZE(4),1);
% noise_estimate = zeros(SIZE(4),1);

if (dim==3)
    h = zeros(3,3,3);
    h(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
    h(:,:,2) = [0 1 0; 1 -6 1; 0 1 0];
    h(:,:,3) = h(:,:,1); % laplace filter
    maxeig=144; % should be (Lx+Ly+Lz)^2 =>> for 1D it was (4)^2, so now (4+4+4)^2 =144;
    
    H = fftn(h,SIZE(1:3));
    
    if(strcmpi(param.File_EXT,'mat'))
        for i=1:SIZE(4)
            k=1;
            while(k<Nit)
                z(:,:,:,i) = clip(z(:,:,:,i) + 1/(lambda*maxeig)*atlas_filter(y_vol(:,:,:,i),h) - atlas_filter(z(:,:,:,i),h,'cconj')./maxeig,atlas);
                if(param.COST_SAVE == 1)
                    %                    temp= squeeze(y_vol(:,:,:,i) - lambda*ifftn(conj(H).*fftn(z(:,:,:,i),SIZE(1:3))));
                    temp = squeeze(y_vol(:,:,:,i) - lambda*atlas_filter(z(:,:,:,i),h,'conj'));
                    %                    temp_derivative2 = ifftn(H.*fftn(temp,SIZE(1:3)));
                    temp_derivative2 = atlas_filter(temp,h);
                    cost(k,i) = lambda*(sqrt(sum(sum(sum(((atlas ==1).*temp_derivative2).^2))))+sqrt(sum(sum(sum(((atlas ==2).*temp_derivative2).^2))))+ ...
                        + sqrt(sum(sum(sum(((atlas ==3).*temp_derivative2).^2))))+ sqrt(sum(sum(sum(((atlas ==4).*temp_derivative2).^2)))))+ ...
                        +sum(sum(sum((squeeze(y_vol(:,:,:,i)) - temp).^2)))/2;
                end
                k=k+1;
            end
            x_vol(:,:,:,i) = y_vol(:,:,:,i) - lambda*atlas_filter(z(:,:,:,i),h,'conj');
        end
        
        if(param.COST_SAVE == 1)
            param.cost_SPATIAL = [param.cost_SPATIAL; cost];
        end
        
    else
        for i=1:SIZE(4)
            k=1;
            while(k<Nit)
                z(:,:,:,i) = clip(z(:,:,:,i) + 1/(lambda*maxeig)*ifftn(H.*fftn(y_vol(:,:,:,i),SIZE(1:3))) - ifftn(H.*conj(H).*fftn(z(:,:,:,i),SIZE(1:3)))/maxeig,atlas);
                k=k+1;
            end
            x_vol(:,:,:,i) = y_vol(:,:,:,i) - lambda*ifftn(conj(H).*fftn(z(:,:,:,i),SIZE(1:3)));
        end
    end
end

for i=1:length(IND(:,1));
    x_out(:,i) = x_vol(IND(i,1),IND(i,2),IND(i,3),:);
end

end

function out = clip(in,atlas)
% find indices for each region
out= zeros(size(in));
for a=1:max(atlas(:)); % number of regions
    IND_r = find(atlas == a);
    %     out(IND_r) = min(norm(in(IND_r),2),1).*in(IND_r)./norm(in(IND_r),2); % may have NAN problems...
    if (norm(in(IND_r),2) > 1)
        out(IND_r) =  in(IND_r)./norm(in(IND_r),2);
        %         out(IND_r) =  sign(in(IND_r))./sqrt(length(IND_r));
    else
        out(IND_r) =  in(IND_r);
    end
end
end




%%% THIS APPLIES PERFECT FILTERING IF YOU HAVE RECTANGULAR REGIONS AND KNOW
%%% YOUR ATLAS, FOR PHANTOM DATA

function atlas_out = atlas_filter(in,h,type)

if (max(in(:))==0)
    return
end
atlas_out = zeros(size(in));

if (nargin == 2)
    for iter_region = 1:max(in(:)),
        temp = in == iter_region;
        r1 = find((~squeeze(all(all(temp==0,2),3))));
        r2 = find((~squeeze(all(all(temp==0,1),3))));
        r3 = find((~squeeze(all(all(temp==0,1),2))));
        if ~isempty(r1) && ~isempty(r2) && ~isempty(r3)
            atlas_out(r1,r2,r3) = ifftn(fftn(h,[length(r1),length(r2),length(r3)]).*fftn(in(r1,r2,r3)));
        end
    end
    %     atlas_out(1:3,1:10,1:10) = ifftn(fftn(h,[3,10,10]).*fftn(in(1:3,1:10,1:10)));
    %     atlas_out(4:10,1:3,1:10) = ifftn(fftn(h,[7,3,10]).*fftn(in(4:10,1:3,1:10)));
    %     atlas_out(4:10,4:10,1:5) = ifftn(fftn(h,[7,7,5]).*fftn(in(4:10,4:10,1:5)));
    %     atlas_out(4:10,4:10,6:10) = ifftn(fftn(h,[7,7,5]).*fftn(in(4:10,4:10,6:10)));
    
elseif (strcmp(type,'cconj'))
    for iter_region = 1:max(in(:)),
        temp = in == iter_region;
        r1 = find((~squeeze(all(all(temp==0,2),3))));
        r2 = find((~squeeze(all(all(temp==0,1),3))));
        r3 = find((~squeeze(all(all(temp==0,1),2))));
        if ~isempty(r1) && ~isempty(r2) && ~isempty(r3)
        atlas_out(r1,r2,r3) = ifftn(fftn(h,[length(r1),length(r2),length(r3)]).*conj(fftn(h,[length(r1),length(r2),length(r3)])).*fftn(in(r1,r2,r3)));
    end; end
elseif (strcmp(type,'conj'))
    for iter_region = 1:max(in(:)),
        temp = in == iter_region;
        r1 = find((~squeeze(all(all(temp==0,2),3))));
        r2 = find((~squeeze(all(all(temp==0,1),3))));
        r3 = find((~squeeze(all(all(temp==0,1),2))));
        if ~isempty(r1) && ~isempty(r2) && ~isempty(r3)
        atlas_out(r1,r2,r3) = ifftn(conj(fftn(h,[length(r1),length(r2),length(r3)])).*fftn(in(r1,r2,r3)));
    end; end
end
end







