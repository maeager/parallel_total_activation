function y_vol = convert2Dto4D(y,SIZE,IND)

y_vol = zeros(SIZE,'single');

for j = 1:length(IND(:,1));
    y_vol(IND(j,1),IND(j,2),IND(j,3),:) = y(:,j);
end