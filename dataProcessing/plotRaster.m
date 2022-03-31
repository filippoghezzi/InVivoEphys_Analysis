function plotRaster(raster,s,condition,varargin)
   
    figure('units','normalized','outerposition',[0 0 1 1]);
    for suIdx=1:size(raster,1)
        subplot(ceil(sqrt(size(raster,1))),ceil(sqrt(size(raster,1))),suIdx)
        for trial=1:size(raster,2)
            scatter(raster{suIdx,trial},ones(size(raster{suIdx,trial},1),1)*trial,'.k')
            hold on 
        end
        
        ax=gca;        
        if strcmp(condition,'Visual') || strcmp(condition,'Visual_K')      
            patch([0 100 100 0],[0 0 size(raster,2) size(raster,2)],'y','FaceAlpha',.3,'EdgeColor','none')
            ax.XLim=[-1000, 4000];
            
        elseif strcmp(condition,'Optotagging')
            patch([0 50 50 0],[0 0 size(raster,2) size(raster,2)],'c','FaceAlpha',.1,'EdgeColor','none')
            ax.XLim=[-100, 200];
        
        elseif strcmp(condition,'VisualOpto')
            patch([0 100 100 0],[0 0 size(raster,2) size(raster,2)],'y','FaceAlpha',.3,'EdgeColor','none')
            patch([-50 150 150 -50],[0 0 size(raster,2) size(raster,2)],'c','FaceAlpha',.1,'EdgeColor','none')
            ax.XLim=[-1000, 4000];
        
        elseif strcmp(condition,'LaserOnly')
            patch([-50 150 150 -50],[0 0 size(raster,2) size(raster,2)],'c','FaceAlpha',.1,'EdgeColor','none')
            ax.XLim=[-1000, 4000];
        end
%             ax.XLim=[-1000, 4000];

        hold off
        ax.YLim=[0, size(raster,2)];
        ax.YTick=[];
        ax.YAxis.Color='none'; 
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
        figname=strcat('Raster - ',condition);
        sgtitle(figname)
        export_fig(fullfile(folder,figname),'-tiff','-transparent')
        close
    end
     
end