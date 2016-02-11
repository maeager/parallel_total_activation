function out = clipcell(in,cellIND)
% # codegen
% find indices for each region
out = in;
for a=1:length(cellIND); % number of regions
    %     out(IND_r) = min(norm(in(IND_r),2),1).*in(IND_r)./norm(in(IND_r),2); % may have NAN problems...
    if ~isempty(cellIND{a})
        %IND_r = cellIND{a};
        norm_in = norm(in(cellIND{a}),2);
        if (norm_in > 1)
            out(cellIND{a}) =  in(cellIND{a})./norm_in;
            %         out(IND_r) =  sign(in(IND_r))./sqrt(length(IND_r));
        end
    end
end
end