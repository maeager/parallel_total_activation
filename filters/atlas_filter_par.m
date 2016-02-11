function atlas_out = atlas_filter_par(in,h,ftype)
%#codegen
coder.extrinsic('fprintf')

atlas_out = zeros(size(in)); %),zeros(size(in)));
if (max(in(:))==0)
    return
end

size(in)
%H = fftn(h,size(in));

if (ftype == 0)  %nargin ==2)
    for iter_region = 1:max(real(in(:))),
        temp = in == iter_region;
        r1 = find((~squeeze(all(all(temp==0,2),3))));
        r2 = find((~squeeze(all(all(temp==0,1),3))));
        r3 = find((~squeeze(all(all(temp==0,1),2))));
        if ~isempty(r1) && ~isempty(r2) && ~isempty(r3)
            fprintf('%d ',length(r1))
            fprintf('%d ',length(r2))
            fprintf('%d ',length(r3))
            atlas_out(r1,r2,r3) = abs(ifftn(fftn(h,[length(r1),length(r2),length(r3)]).*fftn(in(r1,r2,r3))));
%            atlas_out(r1,r2,r3) = abs(ifftn(H.*fftn(in(r1,r2,r3))));
        end
    end
    %     atlas_out(1:3,1:10,1:10) = ifftn(fftn(h,[3,10,10]).*fftn(in(1:3,1:10,1:10)));
    %     atlas_out(4:10,1:3,1:10) = ifftn(fftn(h,[7,3,10]).*fftn(in(4:10,1:3,1:10)));
    %     atlas_out(4:10,4:10,1:5) = ifftn(fftn(h,[7,7,5]).*fftn(in(4:10,4:10,1:5)));
    %     atlas_out(4:10,4:10,6:10) = ifftn(fftn(h,[7,7,5]).*fftn(in(4:10,4:10,6:10)));
    
elseif (ftype == 1)  %strcmp(ftype,'cconj'))
    for iter_region = 1:max(real(in(:))),
        temp = in == iter_region;
        r1 = find((~squeeze(all(all(temp==0,2),3))));
        r2 = find((~squeeze(all(all(temp==0,1),3))));
        r3 = find((~squeeze(all(all(temp==0,1),2))));
        if ~isempty(r1) && ~isempty(r2) && ~isempty(r3)
            length(r1)
            length(r2)
            length(r3)
         atlas_out(r1,r2,r3) = abs(ifftn(fftn(h,[length(r1),length(r2),length(r3)]).*conj(fftn(h,[length(r1),length(r2),length(r3)])).*fftn(in(r1,r2,r3))));
%         atlas_out(r1,r2,r3) = abs(ifftn(H.*conj(H).*fftn(in(r1,r2,r3))));
        end
    end
elseif (ftype == 2) %strcmp(ftype,'conj'))
    for iter_region = 1:max(real(in(:))),
        temp = in == iter_region;
        r1 = find((~squeeze(all(all(temp==0,2),3))));
        r2 = find((~squeeze(all(all(temp==0,1),3))));
        r3 = find((~squeeze(all(all(temp==0,1),2))));
        if ~isempty(r1) && ~isempty(r2) && ~isempty(r3)
            length(r1)
            length(r2)
            length(r3)
           atlas_out(r1,r2,r3) = abs(ifftn(conj(fftn(h,[length(r1),length(r2),length(r3)])).*fftn(in(r1,r2,r3))));
%           atlas_out(r1,r2,r3) = abs(ifftn(conj(H).*fftn(in(r1,r2,r3))));
        end
    end
    
end
end


