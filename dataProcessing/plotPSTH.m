function plotPSTH(s,condition,varargin)
    
    if strcmp(condition,'Visual')
        PSTH = s.PSTHvisual;
    elseif strcmp(condition,'Optotagging')
        PSTH = s.PSTHoptotagging;
    elseif strcmp(condition,'VisualOpto')
        PSTH = s.PSTHvisualOpto;
    elseif strcmp(condition,'LaserOnly')
        PSTH = s.PSTHlaser;
    elseif strcmp(condition,'Visual_K')
        PSTH = s.PSTHvisual_K;
    elseif strcmp(condition,'WhiskerStim')
        PSTH = s.PSTHwhisker;
    elseif strcmp(condition,'WhiskerStim_K')
        PSTH = s.PSTHwhisker_K;
    end
    
    PSTHbins = s.PSTHbins;
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    for suIdx=1:size(PSTH,1)
        subplot(ceil(sqrt(size(PSTH,1))),ceil(sqrt(size(PSTH,1))),suIdx)
        
        plot(PSTHbins+14,smooth(PSTH(suIdx,:),5),'k')
        ax=gca;
        
        if strcmp(condition,'Visual') || strcmp(condition,'Visual_K')       
%             patch([0 100 100 0],[0 0 max(PSTH(suIdx,:)) max(PSTH(suIdx,:))],'y','FaceAlpha',.3,'EdgeColor','none')
            ax.XLim=[-500,3000];
            
        elseif strcmp(condition,'Optotagging')
            patch([0 50 50 0],[0 0 max(PSTH(suIdx,:)) max(PSTH(suIdx,:))],'c','FaceAlpha',.1,'EdgeColor','none')
            ax.XLim=[-100, 200];
        
        elseif strcmp(condition,'VisualOpto')
%             patch([0 100 100 0],[0 0 max(PSTH(suIdx,:)) max(PSTH(suIdx,:))],'y','FaceAlpha',.3,'EdgeColor','none')
%             patch([-50 150 150 -50],[0 0 max(PSTH(suIdx,:)) max(PSTH(suIdx,:))],'c','FaceAlpha',.1,'EdgeColor','none')
            ax.XLim=[-500,3000];
        
        elseif strcmp(condition,'LaserOnly')
%             patch([-50 150 150 -50],[0 0 max(PSTH(suIdx,:)) max(PSTH(suIdx,:))],'c','FaceAlpha',.1,'EdgeColor','none')
            ax.XLim=[-1000,1000];
        end
        
        ax.XLim=[-500,3000];

        ax.YLim=[0, 40];
        ax.XTick=[0:1000:3000];
        ax.XTickLabel=[0:3];
%         ax.XAxis.Visible='off';
        ax.YAxis.TickLength=[0,0];
        ax.YTick=[0,40];
        hold off
        box off
        
        if s.sulayer(suIdx)==1
            Layer='2/3';
        elseif s.sulayer(suIdx)==2
            Layer='4';
        elseif s.sulayer(suIdx)==3
            Layer='5/6';
        end
      
        title(strcat('Unit-',int2str(s.suid(suIdx)),' - L',Layer))
    end
    
    if ~size(varargin,1)==0
        folder=varargin{1,1};
        figname=strcat('PSTH - ',condition);
        sgtitle(figname)
        export_fig(fullfile(folder,figname),'-tiff','-transparent')
        close
    end
     
end   

