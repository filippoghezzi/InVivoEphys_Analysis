function sulayer=findSingleUnitLayer(spike,folder)
    
    load(fullfile(folder,'CSDinfo.mat'))
    
    for unit=1:size(spike.suid,1)
        su=spike.spikeTimes(ismember(spike.spikeTimes(:,1),spike.suid(unit)),:);
        if ismember(su(1,3),CSDinfo.L4)
            sulayer(unit,1)=2;
        elseif su(1,3)<min(CSDinfo.L4)
            sulayer(unit,1)=3;
        elseif su(1,3)>max(CSDinfo.L4)
            sulayer(unit,1)=1;
        end

    end
    
end    
    
    
 
