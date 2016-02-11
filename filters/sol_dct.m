function [dct_sol, c_dct] = sol_dct(in,TR,TS,realign_param)

% detrend realignment parameters and linear trend 
% 07.08.2011

% DCT basis = 198 sec threshold = order(6);
% order = 2*N*TR/TS + 1; (N = dimension, TR = repetition time in sec., TS = period in sec)

% in = acc_l_tc;
n = length(in);

k = round(2*n*TR/TS + 1);

dct = spm_dctmtx(n,k); %construct DCT matrix
lt = (1:n);
lt = lt./sqrt(sum(lt.^2)); % linear trend
dct = [dct , lt' , realign_param];

c_dct = (transpose(dct)*dct)\transpose(dct)*in;  % DCT coefficients

dct_sol = in - dct*c_dct; % detrended

% dct_sol = dct*c_dct;


% PLOT
% figure; plot(in);hold on; plot(dct*c_dct); legend('INPUT','DCT FITTING') % show fitting 

% figure; plot(dct_sol); legend('MINIMUM') %show detrended 

end
