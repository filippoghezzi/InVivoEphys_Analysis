function sulayer=findSingleUnitLayer(s,ops)
        
    for unit=1:numel(s.suid)
        unitChannel=s.cch(ismember(s.cids,s.suid(unit)));
        
        if ismember(unitChannel,ops.L4)
            sulayer(unit,1)=2;
        elseif unitChannel < min(ops.L4)
            sulayer(unit,1)=1;
        elseif unitChannel > max(ops.L4)
            sulayer(unit,1)=3;
        end
    end
    
end    
    
    
 
