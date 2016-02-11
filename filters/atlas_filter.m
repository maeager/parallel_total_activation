function atlas_out = atlas_filter(in,h,type)


atlas_out = zeros(size(in));

if (nargin ==2)
    for iter_region = 1:max(in(:)),
        temp = in == iter_region;
        r1 = find((~squeeze(all(all(temp==0,2),3))));
        r2 = find((~squeeze(all(all(temp==0,1),3))));
        r3 = find((~squeeze(all(all(temp==0,1),2))));
        atlas_out(r1,r2,r3) = ifftn(fftn(h,[length(r1),length(r2),length(r3)]).*fftn(in(r1,r2,r3)));
    end
    %     atlas_out(1:3,1:10,1:10) = ifftn(fftn(h,[3,10,10]).*fftn(in(1:3,1:10,1:10)));
    %     atlas_out(4:10,1:3,1:10) = ifftn(fftn(h,[7,3,10]).*fftn(in(4:10,1:3,1:10)));
    %     atlas_out(4:10,4:10,1:5) = ifftn(fftn(h,[7,7,5]).*fftn(in(4:10,4:10,1:5)));
    %     atlas_out(4:10,4:10,6:10) = ifftn(fftn(h,[7,7,5]).*fftn(in(4:10,4:10,6:10)));
    
elseif (strcmp(type,'cconj'))
    for iter_region = 1:max(in(:)),
        temp = in == iter_region;
        r1 = find((~squeeze(all(all(temp==0,2),3))));
        r2 = find((~squeeze(all(all(temp==0,1),3))));
        r3 = find((~squeeze(all(all(temp==0,1),2))));
        atlas_out(r1,r2,r3) = ifftn(fftn(h,[length(r1),length(r2),length(r3)]).*conj(fftn(h,[length(r1),length(r2),length(r3)])).*fftn(in(r1,r2,r3)));
    end
elseif (strcmp(type,'conj'))
    for iter_region = 1:max(in(:)),
        temp = in == iter_region;
        r1 = find((~squeeze(all(all(temp==0,2),3))));
        r2 = find((~squeeze(all(all(temp==0,1),3))));
        r3 = find((~squeeze(all(all(temp==0,1),2))));
        atlas_out(r1,r2,r3) = ifftn(conj(fftn(h,[length(r1),length(r2),length(r3)])).*fftn(in(r1,r2,r3)));
    end
end
end







