function sulayer=findSingleUnitLayer(spike,folder)
    
    load(fullfile(folder,'CSDinfo.mat'))
    
    for unit=1:size(spike.suid,1)
        su=spike.spikeTimes(ismember(spike.spikeTimes(:,1),spike.suid(unit)),:);
        if ismember(su(unit,3),CSDinfo.L4)
            sulayer(unit,1)=2;
        elseif su(unit,3)<min(CSDinfo.L4)
            sulayer(unit,1)=1;
        elseif su(unit,3)>max(CSDinfo.L4)
            sulayer(unit,1)=3;
        end

    end
    
end    
    
    
 
