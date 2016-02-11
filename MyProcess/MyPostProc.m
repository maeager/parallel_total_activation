% MyPostProc.m
% Generate the activity-inducing signal and innovation signal (if blocks)
% for spikes activity-inducing signal = innovation signal

initial=1;

TC_D_OUT = zeros(param.Dimension(4),param.NbrVoxels); % ACTIVITY-INDUCING SIGNAL

if strcmpi(param.METHOD_TEMP,'block') || strcmpi(param.METHOD_TEMP,'wiener')
    TC_D2_OUT = zeros(param.Dimension(4),param.NbrVoxels); % innovation signal
%    TC_D_OUT2 = zeros(param.Dimension(4),param.NbrVoxels);
end

for i=1:param.NbrVoxels,
	TC_D_OUT(:,i) = filter_boundary(param.f_Recons.num,param.f_Recons.den,TC_OUT(:,i),'normal');
    if strcmpi(param.METHOD_TEMP,'block') || strcmpi(param.METHOD_TEMP,'wiener')
        TC_D2_OUT(:,i) = [0;diff((TC_D_OUT(:,i)))];
%        TC_D_OUT2(:,i) = cumsum([zeros(5,1); TC_D2_OUT(6:end,i)]);  %Neglect the first 5 volumes?? sometimes shifts the response...
    end
end





