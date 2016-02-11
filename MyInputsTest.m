
%MyInputsTest.m

%%% Enter inputs for test data

clear;
close all;
clc;

%% ENTER FILE

file_name = 'test_data/phantom1db'; % enter the folder name of the data
atlas_name = 'test_data/phantom1db'; % enter the folder name of atlas

filename = 'data'; % enter the name of fuctional data
atlasname = 'atlas'; % enter the name of the atlas

path_results = 'test_data/phantom1db'; % enter the path to results

RealignParam = []; 
% RealignParam = importdata(['rp_',filename,'.mat']); % add realignment parameter

param.File_EXT = 'mat'; % nii

%% SELECT PROCESS STEPS;

METHOD_TEMP = 'B' ;% S:spikes, B:blocks, W:Wiener % Type of Temporal Regularization 
METHOD_SPAT = 'STRSPR'; % Tik:tikhonov OR 'StrSpr:Structured Sparsity' OR 'NO': no spatial regularization % Type of Temporal Regularization 
DETRENDING = 'normalize'; % or 'dct' data will be detrended for low frequency oscialltions (DCT).
% NOTE: normalization is voxel by voxel z-normalization 

% select shape of HRF
param.HRF ='bold'; %%% or bold/spmhrf
param.TR = 1; % enter repetition time
param.DCT_TS = 250; % enter dct cut off period (only for DETRENDING='dct')
param.LambdaTempCoef = 1/0.8095; % OR 1/0.6745; estimates higher noise % mad coefficient for temporal regularization
param.COST_SAVE = 1; % save the costs % (optional) for test data

