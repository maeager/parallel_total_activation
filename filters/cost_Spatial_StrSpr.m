function val = cost_Spatial_StrSpr(y,atlas,lambda,temp_derivative2,temp)
%#codegen
val= lambda*(sqrt(sum(sum(sum(((atlas == 1).*temp_derivative2).^2)))) + ...
       sqrt(sum(sum(sum(((atlas == 2).*temp_derivative2).^2)))) + ...
       sqrt(sum(sum(sum(((atlas ==3).*temp_derivative2).^2)))) + ...
       sqrt(sum(sum(sum(((atlas ==4).*temp_derivative2).^2))))) + ...
       sum(sum(sum((squeeze(y) - temp).^2)))/2;
               