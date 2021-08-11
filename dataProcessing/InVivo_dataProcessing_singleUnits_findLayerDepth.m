function s=InVivo_dataProcessing_singleUnits_findLayerDepth(s,ops)
        
    for unit=1:numel(s.suid)
        unitChannel=s.cch(ismember(s.cids,s.suid(unit)));
        
        %% Find layer
        if ismember(unitChannel,ops.L4)
            s.sulayer(unit,1)=2;
        elseif unitChannel < min(ops.L4)
            s.sulayer(unit,1)=1;
        elseif unitChannel > max(ops.L4)
            s.sulayer(unit,1)=3;
        end
        
        %% Find depth
        if strcmp(ops.electrodeType,'Poly2')
            s.suDepth(unit,1) = (ops.L4best-unitChannel)*25;
        elseif strcmp (ops.electrodeType,'Poly3')
%             leftLine=[3,6,9,12,15,18,21,24,27,30];
%             centralLine=[1,leftLine-1,32];
%             rightLine=leftLine+1;
%             
%             if ismember(ops.L4best,centralLine)
%                 l4ChannelIdx=find(ops.L4best==centralLine)-1;
%             elseif ismember(ops.L4best,leftLine)
%                 l4ChannelIdx=find(ops.L4best==leftLine);
%             elseif ismember(ops.L4best,rightLine)
%                 l4ChannelIdx=find(ops.L4best==rightLine);
%             end
%             
%             if ismember(unitChannel,centralLine)
%                 unitChannelIdx=find(unitChannel==centralLine);
%             elseif ismember(unitChannel,leftLine)
%                 unitChannelIdx=find(unitChannel==leftLine);
%             elseif ismember(unitChannel,rightLine)
%                 unitChannelIdx=find(unitChannel==rightLine);
%             end
%             
%             dIdx=l4ChannelIdx-unitChannelIdx;
%             if dIdx==0
            s.suDepth(unit,1)=NaN;
    end
    
end    
    
    
 
