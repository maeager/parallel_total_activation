function [filter_analyze,filter_reconstruct,maxeig] = hrf_bold_spike(TR)
%#codegen

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

%h_d{1} = h_dc;
%h_d{2} = h_dnc;

%% Output setup
filter_reconstruct.num = hnum;
%filter_reconstruct.den = h_d;
filter_reconstruct.den = {h_dc, h_dnc};


    % if spike both filters are the same
    
    [d1,d2] = freqz(hnum,hden,1024);
    maxeig = max(abs(d1).^2);
    filter_analyze = filter_reconstruct;
    
