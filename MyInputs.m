
%MyInputs.m

clear;
close all;
clc;

%% CHANGE THE FOLLOWING FILEDS ACCORDING TO YOUR DATA

file_name = 'PATH TO DATA FOLDER'; % enter the path to data
atlas_name = 'PATH TO ATLAS FOLDER'; % enter the path to atlas

filename = 'DATA NAME'; % enter the name of fuctional data 
atlasname = 'ATLAS NAME'; % enter the name of the atlas

RealignParam = [];
% RealignParam = importdata(['rp_',filename,'.mat']); % add realignment parameters if any

path_results = 'PATH TO RESULTS'; % enter the path to write the results


param.File_EXT = 'nii'; 

%% SELECT PROCESS STEPS;


METHOD_TEMP = 'B' ;% S:spikes, B:blocks, W:Wiener
METHOD_SPAT = 'StrSpr'; % Tik:tikhonov OR 'StrSpr:Structured Sparsity' OR 'NO': no spatial regularization
DETRENDING = 'dct'; % or 'normalize' (z normalization only) if the data is detrended for linear/polyomial trends or DCT coefficients before.

% select shape of HRF
param.HRF ='bold'; %%% or bold/spmhrf
param.TR = 2; % enter repetition time
param.DCT_TS = 125; % 250 % ENTER DCT CUT OFF PERIOD
param.LambdaTempCoef =1/0.8095; % mad coefficient for temporal regularization
param.COST_SAVE = 0; % do not save the costs




