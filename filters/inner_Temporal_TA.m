function [x lambda cost] = inner_Temporal_TA(y, n, d1, d2, noise_estimate, lambda, maxeig,N, Nit) %#codegen

nv = zeros(Nit,1);
%Lambda = zeros(Nit,1);
precision = noise_estimate/100000;

%if (COST_SAVE)
    cost = zeros(Nit,1);
    %end


z = zeros(N,1);

t = 1;
s = zeros(N,1);
for  k = 1:Nit
    
    %%%% estimate for y %%%%%%
    z_l = z;
    z = 1/(lambda*maxeig)*filter_boundary_par_mex(n,d1,d2,y,logical(0)) + s - filter_boundary_par_mex(n,d1,d2,filter_boundary_par_mex(n,d1,d2,s,logical(1)),logical(0))/maxeig;
    % clipping
    z = max(min(z,1),-1);
    t_l = t;
    t = (1+sqrt(1+4*(t^2)))/2;
    s = z + (t_l - 1)/t*(z-z_l);
    
            nv(k) = sqrt(sum((lambda*filter_boundary_par_mex(n,d1,d2,z,1)).^2)/N);
    
    if(abs(nv(k) - noise_estimate)>precision);
        lambda = lambda*noise_estimate/nv(k);
    end
    
 %   Lambda(k) = lambda;
    % k=k+1;
end

x = y - (lambda)*filter_boundary_par_mex(n,d1,d2,z,1); % update x
