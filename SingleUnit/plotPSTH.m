function plotPSTH(PSTH,PSTHbins,spike,varargin)
  
    figure('units','normalized','outerposition',[0 0 1 1]);
    for suIdx=1:size(PSTH,1)
        subplot(ceil(sqrt(size(PSTH,1))),ceil(sqrt(size(PSTH,1))),suIdx)
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
        ax.XLim=[-500, 1000];
        if ~max(PSTH(suIdx,:))==0
            ax.YLim=[0, max(PSTH(suIdx,:))];
        end
%         ax.YTick=[];
%         ax.YAxis.Color='none'; 
        box off
        
        if spike.suLayer(suIdx)==1
            Layer='2/3';
        elseif spike.suLayer(suIdx)==2
            Layer='4';
        elseif spike.suLayer(suIdx)==3
            Layer='5/6';
        end
        
        title(strcat('Unit-',int2str(spike.suid(suIdx)),' - L',Layer))
        

    end
    
    if ~size(varargin,1)==0
        folder=varargin{1,1};
        ID=varargin{1,2};
        age=varargin{1,3};
        figname=strcat(ID,'- P',int2str(age),' - Visual Flash Units - Laser OFF - PSTH - 1');
        sgtitle(figname)
        export_fig(fullfile(folder,figname),'-tiff','-transparent')
        close
    end
     
end   

