function D_OUT = zero_crossings(TC_D_OUT,thresh)

if nargin < 2
    thresh = 0.0;
end

D_OUT=zeros(size(TC_D_OUT));
dsize = size(TC_D_OUT,2);
parfor i=1:dsize
    y=TC_D_OUT(:,i);
    t1 = y(2:end);
    t2 = y(1:end-1);
    dt = t2-t1;
    ind = t2.*t1;
    
    %indx=find(ind<0);
    indx_up   = find( ind<thresh & (dt>0) );
    indx_down = find( ind<thresh & (dt<0) );
    D = zeros(1,length(y));
    D(indx_up)=1;
    D(indx_down)=-1;
    D_OUT(:,i) = D;
end