function plotPSTH(PSTH,PSTHbins,L4)
   

figure('units','normalized','outerposition',[0 0 1 1]);

%         if ismember(unit(1,3),L4)
%             Layer='L4';
%         elseif unit(1,3)<min(L4)
%             Layer='Ingrafranular';         
%         elseif unit(1,3)>max(L4)
%             Layer='Supragranular';         
%         end

for suIdx=1:size(PSTH,1)
    subplot(7,8,suIdx)
    
    plot(PSTHbins,PSTH(suIdx,:),'k')

        
        
    


%         if ~isempty(eventIdx{3}) 
            patch([0 100 100 0],[0 0 max(PSTH(suIdx,:)) max(PSTH(suIdx,:))],'y','FaceAlpha',.3,'EdgeColor','none')
%         end
%         if opto
%             patch([-50 150 150 -50],[0 0 50 50],'c','FaceAlpha',.1,'EdgeColor','none')
%         end
        hold off
%         title(strcat('N',int2str(suid(suIdx)),' - ',Layer),'FontSize',12)
        ax=gca;
        ax.XLim=[-1000, 5000];
%         ax.YLim=[0, size(PSTH,2)];
        ax.YTick=[];
        ax.YAxis.Color='none'; 
        box off
    
end
    
end