
% MyRegularization.m

switch lower(param.METHOD_TEMP)
    case 's'
        param.METHOD_TEMP='SPIKE';
        [param.f_Analyze,param.f_Recons,param.MaxEig] = hrf_filters(param.TR,'spike',param.HRF);
        param.NitTemp = 200;
        fprintf('Temporal Regularization for SPIKES\n');
    case 'b'
        param.METHOD_TEMP = 'BLOCK';
        [param.f_Analyze,param.f_Recons,param.MaxEig] = hrf_filters(param.TR,'block',param.HRF);
        param.NitTemp = 500;
        fprintf('Temporal Regularization for BLOCKS\n');
    case'poss'
        disp('NOT YET IMPLEMENTED!');
    case 'w'
        param.METHOD_TEMP='WIENER';
        [param.f_Analyze,param.f_Recons,param.MaxEig] = hrf_filters(param.TR,'block',param.HRF);
        param.NitTemp = 1;
        fprintf('Temporal Regularization: WIENER FILTER for BLOCKS\n');
    otherwise
        disp('Unknown method.');
end


switch lower(param.METHOD_SPAT)
    case 'no'
        param.METHOD_SPAT = 'no';
        save_path = fullfile(path_results,['OUTS_',param.METHOD_TEMP,'_LamCoef' , ...
            strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.functfilename,'_',date,'.mat']);
        fprintf('NO Spatial Constraint\n');
    case 'tik'
        param.METHOD_SPAT = 'TIK';
        fprintf('Spatial Regularization: Tikhonov 2nd order\n');
        param.Nit=5;
        param.NitSpat=100;
        param.LambdaSpat=1;
        param.stepsize=0.01;
        param.weights = [0.5 0.5]; % weights for Gen Back-Forward
        param.dimTik=3; % Tikhonov in 3d (or 2d).
        save_path = fullfile(param.path_results,['OUTS_',param.METHOD_TEMP,'_LamCoef' , ...
            strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.METHOD_SPAT,...
            '_Lam',num2str(param.LambdaSpat),'_',param.functfilename,'_',date,'.mat']);
    case 'strspr'
        param.METHOD_SPAT = 'STRSPR';
        fprintf('Spatial Regularization: Sptructured Sparsity l_{2-1}-norm\n');
		 % Number of outer iterations
         param.Nit=10;
         param.NitSpat=100;
         % Here Adjust the weight(LambdaSpat) of spatial regularization...
         param.LambdaSpat=5;
         % We assign equal weights for both solutions
         param.weights = [0.5 0.5]; % equal weights for Gen Back-Forward
         param.OrderSpat = 2; % use 2nd order derivative... for now
         param.dimStrSpr=3; %only 3 for now...
        save_path = fullfile(param.path_results,['OUTS_',param.METHOD_TEMP,'_LamCoef' , ...
            strrep(num2str(param.LambdaTempCoef,'%1.1f'),'.',''),'_',param.METHOD_SPAT,...
            '_Lam',num2str(param.LambdaSpat),'_',param.functfilename,'_',date,'.mat']);
        
    otherwise
        disp('Unknown method!');
end

fprintf('All the parameters are saved in param: Check the parameters before proceeding further \n');
param

%reply = input('parameters OK? Y to continue: ','s');
%if ~(strcmpi(reply,'y'))
%    fprintf('Exiting... \n')
%    break;
%end
    

% load('ParneshData_param.mat')

%param.Nit=2;
%param.NitTemp=2;
%param.NitSpat=2;

%% parameters set, START REGULARIZATION

tspat=tic;
[TC_OUT,param] = MySpatial(TCN,atlas,param);
time2 = toc(tspat);

disp(' ');
disp(['IT TOOK ', num2str(time2), ' SECONDS FOR SPATIO_TEMPORAL REGULARIZATION OF ', num2str(param.NbrVoxels), ' TIMECOURSES']);
disp(' ');
param.time=time2;




%load ParneshData_param.mat
%load ParneshData_Spat_99.mat
