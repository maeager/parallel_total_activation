function [TC_OUT,param] = MyTemporal(TCN,param)

% 28.11.2011 isik
% computes temporal regularization for all voxels

fprintf('In MyTemporal\n')
TC_OUT = zeros(param.Dimension(4),param.NbrVoxels,'single');
param.LambdaTemp = zeros(param.NbrVoxels,1,'single');

% paramOUT = cell(1,param.NbrVoxels);
switch lower(param.METHOD_TEMP)
    case{'spike','block'}
        if (param.COST_SAVE)
            costtemp = zeros(param.NitTemp,param.NbrVoxels,'single');
        end
        LambdaTemp=param.LambdaTemp;
        LambdaTempFin=zeros(1,param.NbrVoxels,'single');
        NoiseEstimateFin = zeros(1,param.NbrVoxels,'single');
        parfor i=1:param.NbrVoxels,
        %parfor i=1:1000
            ttemp = tic;
            [coef,len] = wavedec(TCN(:,i),1,'db3');
            coef(1:len(1)) = [];
            LambdaTemp(i) = mad(coef,1)*param.LambdaTempCoef;      
            
            [TC_OUT(:,i),nv,Lambda] = fastTemporal_TA(TCN(:,i),param, i, ...
                                                 LambdaTemp(i));
            LambdaTempFin(i) = Lambda;
            NoiseEstimateFin(i) = nv;
%             if (param.COST_SAVE)
%                 costtemp(:,i) = paramOUT.CostTemp;
%             end
            %toc(ttemp)      
            fprintf('|%d ' , i); toc(ttemp)
        end
        param.LambdaTemp=LambdaTemp;
        param.LambdaTempFin=LambdaTempFin;
        param.NoiseEstimateFin = NoiseEstimateFin;
        
        if (param.COST_SAVE)
            param.cost_TEMP = [param.cost_TEMP; costtemp];
        end
    case{'poss'}
        disp('NOT YET IMPLEMETED for positivity spike...')
    case{'wiener'}
        
        param.fftHnum = abs(fft(param.f_Analyze.num,param.Dimension(4))).^2;
        param.fftHden = abs(fft(param.f_Analyze.den{1},param.Dimension(4)).*fft(param.f_Analyze.den{2},param.Dimension(4)) ...
            .*(param.f_Analyze.den{2}(end)).*(exp((1:param.Dimension(4))*(length(param.f_Analyze.den{2})-1)/param.Dimension(4)))).^2;
        for i=1:param.NbrVoxels,
            [coef,len] = wavedec(TCN(:,i),1,'db3');
            coef(1:len(1)) = [];
            param.LambdaTemp(i) = (mad(coef,1)*param.LambdaTempCoef)^2*param.Dimension(4);
            param.VxlInd = i;
            TC_OUT(:,i) = ifft(fft(TCN(:,i)).*(param.fftHden'./(param.fftHden' + param.fftHnum'.*param.LambdaTemp(i))));
        end
        
        
    otherwise
        fprintf('Unknown temporal method.')
end

end



