function fastActivation(data_path)
%%%%%% Spatio-Temporal Regularization of FMRI Data %%%%%
%
%
% version 1.0, December 2012
% by Isik Karahanoglu 
%
% Modifications by Michael Eager (michael.eager@monash.edu)


% add TA path
addpath(genpath(pwd));
%% Enter the inputs

%% CHANGE THE FOLLOWING FILES ACCORDING TO YOUR DATA

param.data_path = data_path; 
param.File_EXT = 'nii.gz'; 
param.functfilename='rfMRI_REST1_LR'
param.atlasfilename = 'Atlas_wmparc.2'; % enter the name of the atlas

RealignParam = [];
% RealignParam = importdata(['rp_',filename,'.mat']); % add realignment parameters if any

param.path_results = param.data_path; % enter the path to write the results

param.filename = fullfile(param.data_path,[param.functfilename '.' param.File_EXT]); % enter the name of fuctional data 
param.atlasname = fullfile(param.data_path,[param.atlasfilename '.' param.File_EXT]); % enter the name of fuctional data 

%% Check path and files

if ~isdir(data_path)
    disp 'Data path is not valid'
    return;
end

if ~exist(param.filename)
    disp 'Funct filename does not exist'
    return;
end
if ~exist(param.atlasname)
    disp 'Atlas filename does not exist'
    return;
end

%% SELECT PROCESS STEPS;


param.METHOD_TEMP = 'B' ;% S:spikes, B:blocks, W:Wiener
param.METHOD_SPAT = 'StrSpr'; % Tik:tikhonov OR 'StrSpr:Structured Sparsity' OR 'NO': no spatial regularization
param.DETRENDING = 'dct'; % or 'normalize' (z normalization only) if the data is detrended for linear/polyomial trends or DCT coefficients before.

% select shape of HRF
param.HRF ='bold'; %%% or bold/spmhrf
param.TR = 0.8; % enter repetition time
param.DCT_TS = 125; % 250 % ENTER DCT CUT OFF PERIOD
param.LambdaTempCoef =1/0.8095; % mad coefficient for temporal regularization
param.COST_SAVE = 0; % do not save the costs


%param.NitTemp = 10
%param.NitSpat = 10

%% READ DATA

fprintf('Reading the data... \n');

if (strcmpi(param.File_EXT,'nii'))
    [data,hdr] = cbiReadNifti(param.filename); %read data
    [atlas,hdratlas] = cbiReadNifti(param.atlasname); %read atlas
elseif (strcmpi(param.File_EXT,'nii.gz'))
    nii = load_nii(param.filename);data = nii.img;
    anii = load_nii(param.atlasname);atlas = anii.img;
else
    data = importdata(param.filesname);
    atlas = importdata(param.atlasname);
end

%%% ignore the first 10 timecourses...
% data(:,:,:,1:10)=[];

if ~isdir(param.path_results)
    fprintf('Results path does not exist %s \n creating new directory \n', path_results);
    mkdir(param.path_results)
end
param.Dimension = size(data);

%%% Take only non-zeros values both in the atlas and in the data...
%%% REMEMBER: Atlas is gray matter only, but the data is GM+WM+CSF


param.IND = int32(find(atlas.*sum(data,4)));
[param.VoxelIdx(:,1),param.VoxelIdx(:,2),param.VoxelIdx(:,3)] = ind2sub(param.Dimension(1:3),param.IND);
param.VoxelIdx=int16(param.VoxelIdx);

fprintf('%d voxels out of %d voxels are allocated by functional atlas \n', length(find(atlas)), length(find(sum(data,4))));
fprintf('%d voxels are considered, %d voxels of the atlas are out of bound \n\n', length(param.IND),length(find(atlas))-length(param.IND));

param.NbrVoxels = length(param.VoxelIdx(:,1));
TC = zeros(param.NbrVoxels,param.Dimension(4),'single');
for i=1:length(param.VoxelIdx(:,1));
    TC(i,:) = single(data(param.VoxelIdx(i,1),param.VoxelIdx(i,2),param.VoxelIdx(i,3),:)); % timecourses
end

TC=single(TC);
data = single(data);
atlas= int16(atlas);



%% DETREND

MyDetrend
%load('ParneshDataDetrended.mat','TCN');
%fprintf(' Loaded detrended data.\n')



%% SOLVE

MyRegularization

% load ParneshData_Step2.mat
% load ParneshData_param.mat
% 
%        save_path = fullfile(path_results,['OUTS_',param.METHOD_TEMP,'_LamCoef' , ...
%             strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.METHOD_SPAT,...
%             '_Lam',num2str(param.LambdaSpat),'_',filename,'_',date,'.mat']);


%%
MyPostProc
% [TC_OUT,TC_D_OUT,TC_D_OUTN,paramOUT], COST,COSTN

%% SAVE
        save_path = fullfile(param.path_results,['OUTS_',param.METHOD_TEMP,'_LamCoef' , ...
            strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.METHOD_SPAT,...
            '_Lam',num2str(param.LambdaSpat),'_',param.functfilename,'_',date,'.mat']);

display('Saving...');
save(save_path);
display('DONE...');

%% WRITE VOLUME

TCN_vol =  zeros(param.Dimension); % NORMALIZED TIME COURSES
TC_OUT_vol = zeros(param.Dimension); % ACTIVITY-RELATED SIGNALS
TC_D_OUT_vol = zeros(param.Dimension); % ACTIVITY-INDUCED SIGNALS

for i = 1:length(param.VoxelIdx(:,1));
    TCN_vol(param.VoxelIdx(i,1),param.VoxelIdx(i,2),param.VoxelIdx(i,3),:) = TCN(:,i);
    TC_OUT_vol(param.VoxelIdx(i,1),param.VoxelIdx(i,2),param.VoxelIdx(i,3),:)  = TC_OUT(:,i);
    TC_D_OUT_vol(param.VoxelIdx(i,1),param.VoxelIdx(i,2),param.VoxelIdx(i,3),:) = TC_D_OUT(:,i);
end


if strcmpi(param.METHOD_TEMP,'block')
    TC_D2_OUT_vol = zeros(param.Dimension); %INNOVATION SIGNALS
    %innvoation signal
    for i = 1:length(param.VoxelIdx(:,1));
        TC_D2_OUT_vol(param.VoxelIdx(i,1),param.VoxelIdx(i,2),param.VoxelIdx(i,3),:) = TC_D2_OUT(:,i);
    end
    
end

if(strcmpi(METHOD_SPAT,'no'))
    outfilename_tcn = fullfile(param.path_results,['TCN_',param.METHOD_TEMP,'_LamCoef' , ...
        strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.functfilename,'.nii']);
    outfilename_tcout = fullfile(param.path_results,['TC_OUT_',param.METHOD_TEMP,'_LamCoef' , ...
        strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.functfilename,'.nii']);
    outfilename_tcdout = fullfile(param.path_results,['TC_D_OUT_',param.METHOD_TEMP,'_LamCoef' , ...
        strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.functfilename,'.nii']);
    if strcmpi(param.METHOD_TEMP,'block')
        outfilename_tcd2out= fullfile(param.path_results,['TC_D2_OUT_',param.METHOD_TEMP,'_LamCoef' , ...
            strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.functfilename,'.nii']);
    end
else
    outfilename_tcn = fullfile(param.path_results,['TCN_',param.METHOD_TEMP,'_LamCoef' , ...
        strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.METHOD_SPAT,'_Lam', ...
        num2str(param.LambdaSpat),'_',param.functfilename,'.nii']);
    outfilename_tcout = fullfile(param.path_results,['TC_OUT_',param.METHOD_TEMP,'_LamCoef' , ...
        strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.METHOD_SPAT,'_Lam', ...
        num2str(param.LambdaSpat),'_',param.functfilename,'.nii']);
    outfilename_tcdout = fullfile(param.path_results,['TC_D_OUT_',param.METHOD_TEMP,'_LamCoef' , ...
        strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.METHOD_SPAT,'_Lam', ...
        num2str(param.LambdaSpat),'_',param.functfilename,'.nii']);
    if strcmpi(param.METHOD_TEMP,'block')
        outfilename_tcd2out= fullfile(param.path_results,['TC_D2_OUT_',param.METHOD_TEMP,'_LamCoef' , ...
            strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.METHOD_SPAT,'_Lam', ...
            num2str(param.LambdaSpat),'_',param.functfilename,'.nii']);
    end
end

if(exist(outfilename_tcn,'file') || exist(outfilename_tcout,'file') || exist(outfilename_tcdout,'file'))
    if(exist(outfilename_tcn,'file'))
        fprintf('Filename %s already exists, overwriting... \n',outfilename_tcn);
    end
    if(exist(outfilename_tcout,'file'))
        fprintf('Filename %s already exists, overwriting...\n',outfilename_tcout);
    end
    if(exist(outfilename_tcdout,'file'))
        fprintf('Filename %s already exists, overwriting...\n',outfilename_tcdout);
    end
end

if(strcmpi(param.File_EXT,'nii'))
    hdr_pc = hdr;
    hdr_pc.hdr_name = outfilename_tcn;
    hdr_pc.img_name = hdr_pc.hdr_name;
    
    hdr_tcS= hdr;
    hdr_tcS.hdr_name = outfilename_tcout;
    hdr_tcS.img_name = hdr_tcS.hdr_name;
    
    hdr_tcdS = hdr;
    hdr_tcdS.hdr_name = outfilename_tcdout;
    hdr_tcdS.img_name = hdr_tcdS.hdr_name;
    
    if strcmpi(param.METHOD_TEMP,'block')
        hdr_tcd2 = hdr;
        hdr_tcd2.hdr_name = outfilename_tcd2out;
        hdr_tcd2.img_name = hdr_tcd2.hdr_name;
    end
    
    % write data
    %%% NORMALIZED
    [a,b] = cbiWriteNifti(hdr_pc.hdr_name,TCN_vol,hdr_pc,'float32');
    %%% DENOISED
    [a,b] = cbiWriteNifti(hdr_tcS.hdr_name,TC_OUT_vol,hdr_tcS,'float32');
    %%% DECONVOLVED
    [a,b] = cbiWriteNifti(hdr_tcdS.hdr_name,TC_D_OUT_vol,hdr_tcdS,'float32');
    
    if strcmpi(param.METHOD_TEMP,'block')
        [a,b] = cbiWriteNifti(hdr_tcd2.hdr_name,TC_D2_OUT_vol,hdr_tcd2,'float32'); %deconvolved, spikes
    end
    
    
else
    [byte_pc,hdr_pc] = cbiWriteNifti(outfilename_tcn,TCN_vol);
    [byte_tc,hdr_tc] = cbiWriteNifti(outfilename_tcout,TC_OUT_vol);
    [byte_tcd,hdr_tcd] = cbiWriteNifti(outfilename_tcdout,TC_D_OUT_vol);
    if strcmpi(param.METHOD_TEMP,'block')
        [byte1,hdr1] = cbiWriteNifti(outfilename_tcd2out,TC_D2_OUT_vol);
    end
    
end



%exit

