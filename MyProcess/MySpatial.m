function [TC_OUT,param] = MySpatial(TCN,atlas,param)

% 28.11.2011 isik
% computes spatio-temporal regularization for all voxels

% 01.03.2012 isik
% cost computation -- for synthetic data


TC_OUT = zeros(param.Dimension(4),param.NbrVoxels,'single');
% paramOUT = cell(1,param.NbrVoxels);
switch lower(param.METHOD_SPAT)
    case{'no'}
        if (param.COST_SAVE)
            param.cost_TEMP = [];
        end
        param.NitTemp = param.NitTemp*5;
        [TC_OUT,param] = MyTemporal(TCN,param);
    case{'tik'} % fgp_L2.m
        xT = zeros(param.Dimension(4),param.NbrVoxels,'single');%temporal
        xS = zeros(param.Dimension(4),param.NbrVoxels,'single');%spatial
        k=1;
        while (k <= param.Nit)
            [temp,param] = MyTemporal(TC_OUT-xT+TCN,param);
            xT = xT + (temp - TC_OUT); %update temporal, stepsize=1; xT = xT + stepsize*(temp - TC_OUT);
            if(k<param.Nit)
                temp2 = Spatial_Tikhonov(TC_OUT,TC_OUT-xS+TCN,param); % calculates for the whole volume
                xS = xS+(temp2-TC_OUT);
            end
            TC_OUT = xT*param.weights(1)+param.weights(2)*xS;
            k = k+1;
        end
        
        
    case{'strspr'}
        xT = zeros(param.Dimension(4),param.NbrVoxels,'single');%temporal
        xS = zeros(param.Dimension(4),param.NbrVoxels,'single');%spatial
        k=1;
        if (param.COST_SAVE)
            param.cost_TEMP = [];
            param.cost_SPATIAL = [];
        end
        while (k <= param.Nit)
            ttemp = tic;
            param.NitTemp = param.NitTemp+100; % increase temporal step at each iteration for better convergence.
            [temp,param] = MyTemporal(TC_OUT-xT+TCN,param);
            
            xT = xT + (temp - TC_OUT); %update temporal, stepsize=1; xT = xT + stepsize*(temp - TC_OUT);
            if (size(TCN,1)>300)
                save(fullfile(param.path_results,['Temp_' num2str(k) '.mat']),'xT');
            end
            fprintf('Temporal ');toc(ttemp)
            tspat=tic;
            if(k <= param.Nit)
                [temp2,param] = parSpatial_StrSpr(TC_OUT-xS+TCN,atlas,param); % calculates for the whole volume
                xS = xS+(temp2-TC_OUT);
            end
            TC_OUT = xT*param.weights(1)+param.weights(2)*xS;
            if (size(TCN,1)>300)
            save(fullfile(param.path_results,['Spat_' num2str(k) '.mat']),'TC_OUT','xS');
            end
            fprintf('Spatial ');toc(tspat)

            if (param.COST_SAVE)
                param.cost_TOTAL(k) = calculate_totalcost(TC_OUT,TCN,atlas,param);
            end
            %        param.SOL{k} = TC_OUT;
            k = k+1;
        end
        %        param.temp =temp;
        
    otherwise
        error('Unknown Method.')
end

end
