%% Analysis routines for fMRI total activation
%
% - Oct 2013
% - (C) Michael Eager (michael.eager@monash,edu)
%
%%
addpath(genpath('../singleparTotalActivation')); % change the path

%% load the parameters, results and output used in analysis 
load ../Results/OUTS_BLOCK_LamCoef12_STRSPR_Lam5_rfMRI_REST1_LR_11-Oct-2013.mat

%% Correlation analysis of TC_D2_OUT

[R,P,RLO,RUP]=corrcoef(TC_D2_OUT(:,1:10000));

%% Binarising TC_D2_OUT

D2_OUT=zero_crossings(TC_D2_OUT,-0.01);
tic;D2_STD = zero_crossings_std(TC_D2_OUT);toc

%% Binarising D_OUT

D_OUT=zero_crossings(TC_D_OUT,-0.03);
tic;D_STD = zero_crossings_std(TC_D_OUT);toc

%% Setup ATLAS indexes

cellIND = atlas_cellarray(atlas);

%% Map atlas voxels to 2D voxel indexes: TC_D2_OUT 

map = MapAtlasToIND(cellIND,param.IND);

%%
TC_D2_OUT_4D=convert2Dto4D(TC_D2_OUT,param.Dimension,param.VoxelIdx);
D2_OUT_4D=convert2Dto4D(D2_OUT,param.Dimension,param.VoxelIdx);


%% Show cross-correlation coefficients and lags WITHIN atlas group : TC_D2_OUT
%  cross-correlation sequences for all combinations of the columns of x

% V1 left has atlas index 2005
[C lags] = xcorr(TC_D2_OUT(:,map{2005}(:)),100,'coeff');

% cross covariance on a 2D matrix does the same thing
[C] = xcov(TC_D2_OUT(:,map{2005}(:)),'coef');


for ii=1:length(map)
    map_D2_OUT{ii} = D2_OUT(:,map{ii}(:));
end
parfor ii=2000:2100  %1:length(map)
    if ~isempty(map_D2_OUT{ii}) && length(map_D2_OUT) > 2
        fprintf('map %d, voxels %d\n',ii, length(map{ii}))
        try
            tic;WITHIN_XCORR_D2_OUT{ii} = xcorr(map_D2_OUT{ii});toc
        catch ME
        end
    end
end
save('XCORR_D2_OUT_2000s.mat','WITHIN_XCORR_D2_OUT');

WITHIN_CORRCOEF_D2_OUT={};
WITHIN_PVAL_D2_OUT={};;
parfor ii=1:length(map)
    if ~isempty(map_D2_OUT{ii}) && length(map_D2_OUT{ii}) > 2
        fprintf('map %d, voxels %d\n',ii, length(map{ii}))
        %try
            tic;
            [WITHIN_CORRCOEF_D2_OUT{ii} WITHIN_PVAL_D2_OUT{ii}]= corrcoef(map_D2_OUT{ii});toc
        %catch ME
        %end
    end
end
save('CORRCOEF_D2_OUT_.mat','WITHIN_CORRCOEF_D2_OUT','WITHIN_PVAL_D2_OUT');


[R P] = corrcoef(D2_OUT(:,map{2005}(:)));



%% Estimated delay between signals using the information theoretic delay criterion 
% http://www.cs.rug.nl/~rudy/matlab/source/delay.m

addpath('./Analysis/rudy')

% [LAG,CRITERION] = delay(D2_OUT(:,map{2001}(1)),D2_OUT(:,map{2001}(2)));

VXINDEX=2005;  % 2005 = atlas V1
map_D2_LAG={};
for VXINDEX=2000:2100 %1:length(map)
    
    if isempty(map{VXINDEX})
        continue
    end
    num_voxels_ingroup=size(D2_OUT(:,map{VXINDEX}(:)),2);
    D2_OUT_LAG= zeros(num_voxels_ingroup,num_voxels_ingroup);
    parfor ii=1:num_voxels_ingroup
        for jj=2:num_voxels_ingroup
            if jj>ii
                fprintf('%d %d\n',ii, jj)
                [LAG,CRITERION] = delay(D2_OUT(:,map{VXINDEX}(ii))',D2_OUT(:,map{VXINDEX}(jj))');
                D2_OUT_LAG(ii,jj)=LAG;
            end
        end
    end
    map_D2_LAG{VXINDEX} = D2_OUT_LAG;
end

save('../Results/map_D2_LAG.mat','map_D2_LAG')

% note get the p-val of the correlation

            
            

%% Dynamic time-warp measure calculating min-path distance between temporal signals 
%  algorithm from Dan Ellis



VXINDEX=2005;  % 2005 = atlas V1
map_D2_DWT={};
for VXINDEX=2000:2100 %1:length(map)
    
    if isempty(map{VXINDEX})
        continue
    end
    num_voxels_ingroup=size(D2_OUT(:,map{VXINDEX}(:)),2);
    D2_OUT_DWT=zeros(num_voxels_ingroup,num_voxels_ingroup);
    for ii=1:num_voxels_ingroup
        for jj=2:num_voxels_ingroup
            if jj>ii
                fprintf('%d %d\n',ii, jj)
                x1 = find(D2_OUT(:,map{VXINDEX}(ii))>0);
                x2 = find(D2_OUT(:,map{VXINDEX}(jj))>0);
                M = zeros(length(x1),length(x2));
                for iM = 1:length(x1)
                        for jM = 1:length(x2)
                            M(iM,jM)=abs(x1(iM)-x2(jM));
                        end
                end
                [~,~,sum]=dp(M);
                D2_OUT_DWT(ii,jj)=sum;
            end
        end
    end
    map_D2_DWT{VXINDEX} = D2_OUT_DWT;
end

save('../Results/map_D2_DWT.mat','map_D2_DWT')

% note get the p-val of the correlation

            
%% Average,Sum, std   TCN, D and D2 per sections


mean_D2_OUT={};
sum_D2_OUT={};
mean_D_OUT={};
sum_D_OUT={};

mean_TCN_D2_OUT={};
sum_TCN_D2_OUT={};
mean_TCN_D_OUT={};
sum_TCN_D_OUT={};
for VXINDEX=1:length(map)
    
    if isempty(map{VXINDEX})
        continue
    end

    mean_D2_OUT{VXINDEX} = mean(D2_OUT(:,map{VXINDEX})');
    sum_D2_OUT{VXINDEX} = sum(D2_OUT(:,map{VXINDEX})');
    mean_D2_OUT{VXINDEX} = mean(D_OUT(:,map{VXINDEX})');
    sum_D2_OUT{VXINDEX} = sum(D_OUT(:,map{VXINDEX})');
    mean_TCN_D2_OUT{VXINDEX} = mean(TCN_D2_OUT(:,map{VXINDEX})');
    sum_TCN_D2_OUT{VXINDEX} = sum(TCN_D2_OUT(:,map{VXINDEX})');
    mean_TCN_D2_OUT{VXINDEX} = mean(TCN_D_OUT(:,map{VXINDEX})');
    sum_TCN_D2_OUT{VXINDEX} = sum(TCN_D_OUT(:,map{VXINDEX})');

end

save('../Results/mean_D2.mat','mean_D2_OUT','sum_D2_OUT','mean_D_OUT','sum_D_OUT','mean_TCN_D2_OUT','sum_TCN_D2_OUT','mean_TCN_D_OUT','sum_TCN_D_OUT')

% note get the p-val of the correlation

            
