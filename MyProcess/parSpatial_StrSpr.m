function [x_out,param] = parSpatial_StrSpr(y,atlas,param)

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

disp('In Spatial structured sparsity')

lambda = single(param.LambdaSpat); %regularization parameter
Nit = param.NitSpat; %Number of iterations
IND = param.VoxelIdx; %[xi,yi,zi] = nonzero indices to convert 2D x to 4D x_vol
SIZE=single(param.Dimension); %the size of the 4D data [v1,v2,v3,t]
order=single(param.OrderSpat); %order of the derivative filter 'first', 'second', 'zero'
dim =param.dimStrSpr; %dimension (2,3 filter in 2D or 3D)

if (length(SIZE) ~= 4)
    error('SIZE should have 4 dimensions.')
end

x_out = zeros(size(y),'single'); % out

%%%% convert to 4D
x_vol = zeros(SIZE,'single');
y_vol = zeros(SIZE,'single');

if(strcmpi(param.File_EXT,'mat') && (param.COST_SAVE == 1))
    temp_derivative2 = zeros(SIZE(1:3),'single');
    temp = zeros(SIZE(1:3),'single');
end

for j = 1:length(IND(:,1));
    y_vol(IND(j,1),IND(j,2),IND(j,3),:) = y(:,j);
end
%%%%

cost = zeros(Nit,SIZE(4),'single');
z=[];
%z = zeros(size(y_vol),'single');
cellIND = atlas_cellarray(atlas);

% lambda = zeros(SIZE(4),1);
% noise_estimate = zeros(SIZE(4),1);

if (dim==3)
    h = zeros(3,3,3);
    h(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
    h(:,:,2) = [0 1 0; 1 -6 1; 0 1 0];
    h(:,:,3) = h(:,:,1); % laplace filter
    h=single(h);
    maxeig=144; % should be (Lx+Ly+Lz)^2 =>> for 1D it was (4)^2, so now (4+4+4)^2 =144;
    
    H = single(fftn(h,SIZE(1:3)));
    HHbar = single(H.*conj(H));
    conjH=conj(H);
    
    if(strcmpi(param.File_EXT,'mat'))
        z = zeros(size(y_vol),'single');
        parfor i=1:SIZE(4)
                %tic;
                %z(:,:,:,i) = strsparsity_par_mex(z(:,:,:,i),y_vol(:,:,:,i),h,lambda,maxeig,atlas);
                %toc; disp('MEX sparsity done');
            
            %k=1;
            
            %while(k<Nit)
            for k=1:Nit-1
                tic;z(:,:,:,i) = clip(z(:,:,:,i) + 1/(lambda*maxeig)*atlas_filter(y_vol(:,:,:,i),h) - atlas_filter(z(:,:,:,i),h,'cconj')./maxeig,atlas);toc;

                if(param.COST_SAVE == 1)
                    %                    temp= squeeze(y_vol(:,:,:,i) - lambda*ifftn(conj(H).*fftn(z(:,:,:,i),SIZE(1:3))));
                    %temp = squeeze(y_vol(:,:,:,i) - lambda*atlas_filter_par_mex(z(:,:,:,i),h,2);  % 'conj'));
                    temp = squeeze(y_vol(:,:,:,i) - lambda*atlas_filter_par(z(:,:,:,i),h,2));  % 'conj'));
                    %                    temp_derivative2 = ifftn(H.*fftn(temp,SIZE(1:3)));
                    %temp_derivative2 = atlas_filter_par_mex(temp,h,0);
                    temp_derivative2 = atlas_filter_par(temp,h,0);
                    cost(k,i) = lambda*(sqrt(sum(sum(sum(((atlas == 1).*temp_derivative2).^2)))) + ...
                    sqrt(sum(sum(sum(((atlas == 2).*temp_derivative2).^2))))+ ...
                        + sqrt(sum(sum(sum(((atlas ==3).*temp_derivative2).^2))))+ sqrt(sum(sum(sum(((atlas ==4).*temp_derivative2).^2)))))+ ...
                        + sum(sum(sum((squeeze(y_vol(:,:,:,i)) - temp).^2)))/2;
                    %cost(k,i) = cost_Spatial_StrSpr_mex(y_vol(:,:,:,i),atlas,temp_derivative2,temp);
                end
                %k=k+1;
            end
            x_vol(:,:,:,i) = y_vol(:,:,:,i) - lambda*atlas_filter(z(:,:,:,i),h,'conj');
        end
        
        if(param.COST_SAVE == 1)
            param.cost_SPATIAL = [param.cost_SPATIAL; cost];
        end
        
    else
        parfor i=1:SIZE(4)
        %parfor i=1:8*4    
            %tic;
            %x_vol(:,:,:,i) = strsparsity_nii_mex(y_vol(:,:,:,i),H,lambda,maxeig,atlas);
            %toc; disp('MEX sparsity done');
            
            %k=1;
            %while(k<Nit)

            %tic
            %z = zeros(SIZE(1:3),'single');
            %load('wisdomTA.mat');fftw('wisdom',str);
            ytmp = y_vol(:,:,:,i);
            tmp1 = 1/(lambda*maxeig)*ifftn(H.*fftn(ytmp,SIZE(1:3)));
            z = clipcell(tmp1,cellIND);
            
            for k=2:Nit-1
                %tmp2 = fftn(z,SIZE(1:3));
                %tmp3 = ifftn(HHbar.*fftn(z,SIZE(1:3)))./maxeig;
                z = clipcell(z + tmp1 - ifftn(HHbar.*fftn(z,SIZE(1:3))./maxeig),cellIND);
                %k=k+1;
            end
            x_vol(:,:,:,i) = ytmp - lambda*ifftn(conjH.*fftn(z,SIZE(1:3)));
            %toc
            
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


atlas_out = zeros(size(in));

if (nargin ==2)
    for iter_region = 1:max(in(:)),
        temp = in == iter_region;
        r1 = find((~squeeze(all(all(temp==0,2),3))));
        r2 = find((~squeeze(all(all(temp==0,1),3))));
        r3 = find((~squeeze(all(all(temp==0,1),2))));
        atlas_out(r1,r2,r3) = ifftn(fftn(h,[length(r1),length(r2),length(r3)]).*fftn(in(r1,r2,r3)));
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
        atlas_out(r1,r2,r3) = ifftn(fftn(h,[length(r1),length(r2),length(r3)]).*conj(fftn(h,[length(r1),length(r2),length(r3)])).*fftn(in(r1,r2,r3)));
    end
elseif (strcmp(type,'conj'))
    for iter_region = 1:max(in(:)),
        temp = in == iter_region;
        r1 = find((~squeeze(all(all(temp==0,2),3))));
        r2 = find((~squeeze(all(all(temp==0,1),3))));
        r3 = find((~squeeze(all(all(temp==0,1),2))));
        atlas_out(r1,r2,r3) = ifftn(conj(fftn(h,[length(r1),length(r2),length(r3)])).*fftn(in(r1,r2,r3)));
    end
end
end







