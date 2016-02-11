function map = MapAtlasToIND(cellIND,IND)

map = {};%cellIND;

for icell=1:length(cellIND)
    jmapindex=1;
    for jindex=1:length(cellIND{icell})
        ind = find(IND==cellIND{icell}(jindex));
        if isempty(ind)
            fprintf('Ind empty: %d\t%d\n',icell,jindex);
        else
            if length(ind)==1
                map{icell}(jmapindex) = ind;
                jmapindex=jmapindex+1;
            else
                fprintf('Ind not equal to 1: %d\t%d\t%d\n',icell,jindex,length(ind));
            end
        end
    end

end
