function [x,ParametersOut] = Temporal_TA(y,ParametersIn)

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


n = ParametersIn.f_Analyze.num;
d = ParametersIn.f_Analyze.den;
maxeig = ParametersIn.MaxEig;
N = ParametersIn.Dimension(4);
Nit = ParametersIn.NitTemp;

% if estimated before, take previous lambdas
if (isfield(ParametersIn,'NoiseEstimateFin') && (length(ParametersIn.NoiseEstimateFin)>=ParametersIn.VxlInd))
        lambda = ParametersIn.NoiseEstimateFin(ParametersIn.VxlInd);
else
    lambda = ParametersIn.LambdaTemp(ParametersIn.VxlInd);
end

% std for noise;
noise_estimate = ParametersIn.LambdaTemp(ParametersIn.VxlInd);
% std has to be < 1... only for initial guess...
%if(noise_estimate>0.95), noise_estimate = 0.95; end

nv = zeros(Nit,1);
Lambda = zeros(Nit,1);
precision = noise_estimate/100000;

if (ParametersIn.COST_SAVE)
    cost = zeros(Nit,1);
end


z = zeros(N,1);
k = 1;
t = 1;
s = zeros(N,1);
while (k <= Nit)
    
    %%%% estimate for y %%%%%%
    z_l = z;
    z = 1/(lambda*maxeig)*filter_boundary(n,d,y,'normal') + s - filter_boundary(n,d,filter_boundary(n,d,s,'transpose'),'normal')/maxeig;
    % clipping
    z = max(min(z,1),-1);
    t_l = t;
    t = (1+sqrt(1+4*(t^2)))/2;
    s = z + (t_l - 1)/t*(z-z_l);
    if (ParametersIn.COST_SAVE)
        temp = y - (lambda)*filter_boundary(n,d,z,'transpose'); % update x
        cost(k) = sum((temp-y).^2)/2 + lambda*sum(abs(filter_boundary(n,d,temp,'normal')));
        
        % All nv s are the same
        % nv(k) = sqrt(sum((lambda*filter_boundary(n,d,z,'transpose')).^2)/N);
        % nv(k) = std((temp-y));
        
        nv(k) = sqrt(sum((temp-y).^2)/N);
        
    else
        nv(k) = sqrt(sum((lambda*filter_boundary(n,d,z,'transpose')).^2)/N);
    end
    
    if(abs(nv(k) - noise_estimate)>precision);
        lambda = lambda*noise_estimate/nv(k);
    end
    
    Lambda(k) = lambda;
    k=k+1;
end

x = y - (lambda)*filter_boundary(n,d,z,'transpose'); % update x

% ParametersOut.NoiseEstimateIn = noise_estimate;
ParametersOut.NoiseEstimateFin = nv(end);
ParametersOut.LambdasTempFin = Lambda(end);

if (ParametersIn.COST_SAVE)
    ParametersOut.CostTemp = cost;
end

end








