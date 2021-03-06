function [filter_analyze,filter_reconstruct,maxeig] = hrf_filters(TR,condition,condition2)

% THIS FUNCTION GENERATES THE FILTERS FOR ANALYSIS AND RECONSTRUTION OF HRF
% MODEL, CONDITION IS EITHER BLOCK/SPIKE
% SPIKE MODEL: TWO FILTERS ARE THE SAME
% BLOCK MODEL: RECONSTRCUTION MODEL HAS ONE MORE ZERO = 0.

% 20.12.2010

% INPUTS
% TR = repetition time of FMRI
% condition = block/spike

% OUTPUTS
% filter_analyze = analysis filter
% filter_reconstruct =reconstruction filter
% maxeig = maximum eigenvalue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%EDIT condition2, change cons_filter -->  cons_filter2
% EDIT: condition 2 is spmhrf/bold


if strcmpi(condition2,'bold')
    
    % Constants for BOLD (Friston et al.)
    
    eps = 0.54;
    ts = 1.54;
    tf = 2.46;
    t0 = 0.98;
    alpha = 0.33;
    E0 = 0.34;
    V0 = 1;
    k1 = 7*E0;
    k2 = 2;
    k3 = 2*E0 - 0.2;
    
    c = (1 + (1 - E0)*log(1 - E0)/E0)/t0;
    
    % zeros
    a1 = -1/t0;
    a2 = -1/(alpha*t0);
    a3 = -(1 + 1i*sqrt(4*ts^2/tf - 1))/(2*ts);
    a4 = -(1 - 1i*sqrt(4*ts^2/tf - 1))/(2*ts);
    % pole
    psi = -((k1+k2)*((1 - alpha)/alpha/t0 - c/alpha) - (k3 - k2)/t0)/(-(k1 + k2)*c*t0 - k3 + k2);
    % psti = 0;
    % constant = V0*eps/t0*(-(k1 + k2)*c*t0 - k3 + k2);
    % psi=[];
    
elseif strcmpi(condition2,'spmhrf')
    
    
    %%%%%%%%%%%%%%%%%%
    % CONSTRUCT MODEL
    % FilZeros = [a1,a2,a3,a4]*TR;
    % FilPoles = psi*TR;
    % FilPoles = [];
    
    
    % CONSTRUCT MODEL (SPM_HRF)
    % FilZeros = [a1*4,a1*4,real(a3)*0.7+i*imag(a3)*0.55,real(a4)*0.7+i*imag(a4)*0.55]*TR;
    % FilPoles = psi*TR;
    a1 = -0.27;
    a2 = -0.27;
    a3 =-0.4347-1i*0.3497;
    a4 = -0.4347+1i*0.3497;
    psi = -0.1336;
    %%%%%%%%%%%%%%%%%%
    
    
else
    error('Unknown filter')
end

FilZeros = [a1,a2,a3,a4]*TR;
FilPoles = psi*TR;

% normalize for HRF, analog/digital

% cons=1/(norm((1-exp(a1*TR))*(1-exp(a2*TR))*(1-exp(a3*TR))*(1-(exp(a4*TR)))/(1-(exp(psi*TR))))/norm(a1*a2*a3*a4/psi)); % first order
% cons=1/(norm((3/2-2*exp(a1*TR)+exp(2*a1*TR)/2)*(3/2-2*exp(a2*TR)+exp(2*a2*TR)/2)*(3/2-2*exp(a3*TR)+exp(2*a3*TR)/2)*(3/2-2*exp(a4*TR)+exp(2*a4*TR)/2)/(3/2-2*exp(psi*TR)+exp(2*psi*TR)/2))/norm(a1*a2*a3*a4/psi)); % first order
% cons=1/(norm((1-exp(a1*TR))*(1-exp(a2*TR))*(1-exp(a3*TR))*(1-(exp(a4*TR)))/(1-(exp(psi*TR))))/norm(a1*a2*a3*a4/psi)); %first order

cons=1;
hnum = cons_filter(FilZeros)*cons;
hden = cons_filter(FilPoles);

% second order approximation
%     cons=1/(norm((3/2-2*exp(a1*TR)+exp(2*a1*TR)/2)*(3/2-2*exp(a2*TR)+exp(2*a2*TR)/2)*(3/2-2*exp(a3*TR)+exp(2*a3*TR)/2)*(3/2-2*exp(a4*TR)+exp(2*a4*TR)/2)/(3/2-2*exp(psi*TR)+exp(2*psi*TR)/2))/norm(a1*a2*a3*a4/psi)); % first order
%     hnum = cons_filter2(FilZeros)*cons;
%     hden = cons_filter2(FilPoles);


causal = FilPoles(real(FilPoles)<0);
n_causal = FilPoles(real(FilPoles)>0);

% Shortest Filter, 1st order approximation
h_dc = cons_filter(causal);
h_dnc = cons_filter(n_causal);
%     % second order approximation
%     h_dc = cons_filter2(causal);
%     h_dnc = cons_filter2(n_causal);

h_d{1} = h_dc;
h_d{2} = h_dnc;

filter_reconstruct.num = hnum;
filter_reconstruct.den = h_d;


if( strcmpi(condition,'spike'))
    % if spike both filters are the same
    
    [d1,d2] = freqz(hnum,hden,1024);
    maxeig = max(abs(d1).^2);
    filter_analyze = filter_reconstruct;
    
elseif(strcmpi(condition,'block'))
    FilZeros2 = [FilZeros,0];
    
    %     hnum2 = cons_filter(FilZeros2)*cons*2*sqrt(2);
    
    % Shortest Filter, 1st order approximation
    hnum2 = cons_filter(FilZeros2)*cons;
    % second order approximation
    %         hnum2 = cons_filter2(FilZeros2)*cons;
    
    [d1,d2] = freqz(hnum2,hden,1024);
    maxeig = max(abs(d1).^2);
    
    filter_analyze.num = hnum2;
    filter_analyze.den = h_d;
    
else
    error('Unknown type of condition: Should be "SPIKE" or "BLOCK"')
end

% maxeig = maxeig./filternorm(hnum,h_d{1},2)^2;




