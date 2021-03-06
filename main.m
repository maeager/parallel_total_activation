
%%%%%% Spatio-Temporal Regularization of FMRI Data %%%%%
%
%
% version 1.0, December 2012
% by Isik Karahanoglu 
%
% Modifications by Michael Eager (michael.eager@monash.edu)

% add TA path
	addpath(genpath(pwd)); % change the path

%% Enter the inputs
MyInputsTest %For phantom data
%MyInputs 

%param.NitTemp = 100
%param.NitSpat = 100

%% READ DATA
path_data = fullfile(file_name,[filename,'.',param.File_EXT]);
path_atlas = fullfile(atlas_name,[atlasname,'.',param.File_EXT]); %functional atlas


fprintf('Reading the data... \n');

%%% CHECK PATHS
fprintf('%s \n\n',file_name);
%fprintf('Exist dir %d \n',exist(atlas_name,'dir'));

%if (~exist(file_name,'dir'))
%    error(['Path does not exist : ' file_name ]);
if (~exist(path_data,'file'));
    error([filename ,'.',param.File_EXT,' can not be found in  ' file_name]);
end

%if (~exist(atlas_name,'dir'))
%    error(['Path does not exist : ' atlas_name ]);
if (~exist(path_atlas,'file'))
    error([atlasname ,'.',param.File_EXT, ' can not be found in  ' atlas_name]);
end

if (strcmpi(param.File_EXT,'nii'))
    [data,hdr] = cbiReadNifti(path_data); %read data
    [atlas,hdratlas] = cbiReadNifti(path_atlas); %read atlas
else
    data = importdata(path_data);
    atlas = importdata(path_atlas);
end

%%% ignore the first 10 timecourses...
% data(:,:,:,1:10)=[];

if (~exist(path_results,'dir'))
    fprintf('Results path does not exist %s \n creating new directory \n', path_results);
    mkdir(path_results)
end
param.Dimension = size(data);

%%% Take only non-zeros values both in the atlas and in the data...
%%% REMEMBER: Atlas is gray matter only, but the data is GM+WM+CSF

param.IND = find(atlas.*sum(data,4));
[param.VoxelIdx(:,1),param.VoxelIdx(:,2),param.VoxelIdx(:,3)] = ind2sub(param.Dimension(1:3),param.IND);


fprintf('%d voxels out of %d voxels are allocated by functional atlas \n', length(find(atlas)), length(find(sum(data,4))));
fprintf('%d voxels are considered, %d voxels of the atlas are out of bound \n\n', length(param.IND),length(find(atlas))-length(param.IND));

param.NbrVoxels = length(param.VoxelIdx(:,1));
TC = zeros(param.NbrVoxels,param.Dimension(4));
for i=1:length(param.VoxelIdx(:,1));
    TC(i,:) = data(param.VoxelIdx(i,1),param.VoxelIdx(i,2),param.VoxelIdx(i,3),:); % timecourses
end


%% DETREND

MyDetrend

%% SOLVE
MyRegularization

%%
MyPostProc
% [TC_OUT,TC_D_OUT,TC_D_OUTN,paramOUT], COST,COSTN

%% SAVE

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
    outfilename_tcn = fullfile(path_results,['TCN_',param.METHOD_TEMP,'_LamCoef' , ...
        strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',filename,'.nii']);
    outfilename_tcout = fullfile(path_results,['TC_OUT_',param.METHOD_TEMP,'_LamCoef' , ...
        strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',filename,'.nii']);
    outfilename_tcdout = fullfile(path_results,['TC_D_OUT_',param.METHOD_TEMP,'_LamCoef' , ...
        strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',filename,'.nii']);
    if strcmpi(param.METHOD_TEMP,'block')
        outfilename_tcd2out= fullfile(path_results,['TC_D2_OUT_',param.METHOD_TEMP,'_LamCoef' , ...
            strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',filename,'.nii']);
    end
else
    outfilename_tcn = fullfile(path_results,['TCN_',param.METHOD_TEMP,'_LamCoef' , ...
        strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.METHOD_SPAT,'_Lam', ...
        num2str(param.LambdaSpat),'_',filename,'.nii']);
    outfilename_tcout = fullfile(path_results,['TC_OUT_',param.METHOD_TEMP,'_LamCoef' , ...
        strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.METHOD_SPAT,'_Lam', ...
        num2str(param.LambdaSpat),'_',filename,'.nii']);
    outfilename_tcdout = fullfile(path_results,['TC_D_OUT_',param.METHOD_TEMP,'_LamCoef' , ...
        strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.METHOD_SPAT,'_Lam', ...
        num2str(param.LambdaSpat),'_',filename,'.nii']);
    if strcmpi(param.METHOD_TEMP,'block')
        outfilename_tcd2out= fullfile(path_results,['TC_D2_OUT_',param.METHOD_TEMP,'_LamCoef' , ...
            strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.METHOD_SPAT,'_Lam', ...
            num2str(param.LambdaSpat),'_',filename,'.nii']);
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

