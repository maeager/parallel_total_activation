function cost = diff_dwt(x1,x2)

cost = zeros(length(x1),length(x2));
for i=1:length(x1) 
    for j=1:length(x2)
        cost(i,j)=abs(x1(i)-x2(j)); 
    end
end