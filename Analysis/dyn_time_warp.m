function DWT = dyn_time_warp(D2_OUT)

DWT={};

%Lower triangular matrix
for ii=1:size(D2_OUT,1)
    for jj=2:size(D2_OUT,1)
        if ii<jj
            x1 = find(D2_OUT(:,ii)>0);
            x2 = find(D2_OUT(:,jj)>0);
            diffx12 = diff_dwt(x1,x2);
            
            [DWT_p{ii,jj} DWT_q{ii,jj}] = dpcore(diffx12);
            
        end
    end
end
