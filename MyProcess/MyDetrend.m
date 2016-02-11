
% MyDetrend.m

fprintf('Detrending the data... \n');

TCN = zeros(param.Dimension(4),param.NbrVoxels,'single'); % normalize by dividing to std (var=1).

tic;% start timing

if strcmpi(param.DETRENDING,'DCT')
    param.METHOD_DETREND = 'DCT and normalize';
    fprintf('Detrending method: DCT \n\n');
    for i=1:param.NbrVoxels;
%         [TCN(:,i), coef_tc(:,i)] = sol_dct(TC(i,:)',param.TR,param.DCT_TS,RealignParam); % subtract mean liner detrend + DCT
        [TCN(:,i), ~] = sol_dct(TC(i,:)',param.TR,param.DCT_TS,RealignParam); % subtract mean liner detrend + DCT
        TCN(:,i) = single(TCN(:,i)./std(TCN(:,i)));
    end
    
elseif strcmpi(param.DETRENDING,'normalize')
    param.METHOD_DETREND = 'normalize';
    
    fprintf('No Detrending: only normalize \n\n')
    if(strcmp(param.File_EXT,'mat'))
        for i=1:param.NbrVoxels;
            TCN(:,i) = single((TC(i,:))'./std(TC(i,:)));
        end
    else
        for i=1:param.NbrVoxels;
            TCN(:,i) = single((TC(i,:)-mean(TC(i,:)))'./std(TC(i,:)));
        end
    end
else
    error('Unknown detrending method...');
end
time_detrend = toc;
fprintf('It took %f seconds to detrend (method %s)  %d voxel timecouses \n\n', time_detrend, param.DETRENDING,param.NbrVoxels);
