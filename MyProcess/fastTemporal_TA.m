function [x,nv,Lambda] = fastTemporal_TA(y,ParametersIn, VoxelIdx,LambdaTempAtIdx)

% 28.11.2011 by isik
% the same as fgp_general_v2.m, only adapted for fmri.
% updated lamda as in (chambolle,Lecture notes on Sparsity)
% Total variation filtering (denoising)

%%%%%%---INPUTS---%%%%%%%

% y - noisy signal, one voxel time course
% n = ParametersIn.f_Analyze.num;
% d = ParametersIn.f_Analyze.den;
% maxeig = Param.MaxEig;
% N = ParametersIn.Dimension(4);
% Nit = ParametersIn.NitTemp;
% lambda = ParametersIn.LambdaTemp(ParametersIn.VxlNbr);
% ParametersIn.VxlNbr = voxel number

%%%%%%---OUTPUTS---%%%%%%

% x - result of denoising
% ParametersOut.NoiseEstimateInitial = noise_estimate;
% ParametersOut.NoiseEstimateFin = nv;
% ParametersOut.Cost = cost;
% ParametersOut.LambdasTempFin = Lambda;

%%%%%%-------------%%%%%%
%%% NOTE ON 'transpose'; 04.04.2011
%%% old version = filter_boundary(fliplr(n),d,s,'transpose')
%%% new version = filter_boundary(n,d,s,'transpose')
%%%%%%%%%%%%%%%%%

% NEW version 31.11.2012, check noise_estimate,
% noise_estimate has to be < 1, data is normalized...


n = single(ParametersIn.f_Analyze.num);
d = ParametersIn.f_Analyze.den;
d1 = single(d{1}); d2=single(d{2});
maxeig = single(ParametersIn.MaxEig);
N = single(ParametersIn.Dimension(4));
Nit = single(ParametersIn.NitTemp);

% if estimated before, take previous lambdas
if (isfield(ParametersIn,'NoiseEstimateFin') && (length(ParametersIn.NoiseEstimateFin)>=VoxelIdx))
        lambda = single(ParametersIn.NoiseEstimateFin(VoxelIdx));
else
    lambda = single(LambdaTempAtIdx);
end

% std for noise;
noise_estimate = single(LambdaTempAtIdx);
% std has to be < 1... only for initial guess...
%if(noise_estimate>0.95), noise_estimate = 0.95; end

nv = single(0);
Lambda = single(0);
precision = single(noise_estimate/100000);

% if (ParametersIn.COST_SAVE)
%     cost = zeros(Nit,1,'single');
% end


z = zeros(N,1,'single');
k = single(1);
t = single(1);
s = zeros(N,1,'single');

while (k <= Nit)
    
    %%%% estimate for y %%%%%%
    z_l = z;
    z = 1/(lambda*maxeig)*filter_boundary_single_mex(n,d1,d2,y,single(0)) + s - filter_boundary_single_mex(n,d1,d2,filter_boundary_single_mex(n,d1,d2,s,single(1)),single(0))/maxeig;
    % clipping
    z = max(min(z,1),-1);
    t_l = t;
    t = (1+sqrt(1+4*(t^2)))/2;
    s = z + (t_l - 1)/t*(z-z_l);
    if (ParametersIn.COST_SAVE)
        temp = y - (lambda)*filter_boundary_single_mex(n,d1,d2,z,single(1)); % update x
        cost(k) = sum((temp-y).^2)/2 + lambda*sum(abs(filter_boundary_single_mex(n,d1,d2,temp,single(0))));
        
        % All nv s are the same
        % nv(k) = sqrt(sum((lambda*filter_boundary(n,d,z,'transpose')).^2)/N);
        % nv(k) = std((temp-y));
        
        nv = sqrt(sum((temp-y).^2)/N);
        
    else
        nv = sqrt(sum((lambda*filter_boundary_single_mex(n,d1,d2,z,single(1))).^2)/N);
    end
    
    if(abs(nv - noise_estimate)>precision);
        lambda = lambda*noise_estimate/nv;
    end
    
    Lambda = lambda;
    k=k+1;
end

x = y - (lambda)*filter_boundary_single_mex(n,d1,d2,z,single(1)); % update x

% ParametersOut.NoiseEstimateIn = noise_estimate;
%ParametersOut.NoiseEstimateFin = nv;
%ParametersOut.LambdasTempFin = Lambda;


end








