function out = clip(in,atlas)
%#codegen
% find indices for each region
out = in;
for a=1:max(atlas(:)); % number of regions
    IND_r = find(atlas == a);
    %     out(IND_r) = min(norm(in(IND_r),2),1).*in(IND_r)./norm(in(IND_r),2); % may have NAN problems...
    if ~isempty(IND_r)
        norm_in = norm(in(IND_r),2);
        if (norm_in > 1)
            out(IND_r) =  in(IND_r)./norm_in;
            %         out(IND_r) =  sign(in(IND_r))./sqrt(length(IND_r));
        end
    end
end
end