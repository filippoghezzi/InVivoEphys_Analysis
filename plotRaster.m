function plotRaster(raster,spike,condition,varargin)
   
    figure('units','normalized','outerposition',[0 0 1 1]);
    for suIdx=1:size(raster,1)
        subplot(ceil(sqrt(size(raster,1))),ceil(sqrt(size(raster,1))),suIdx)
        for trial=1:size(raster,2)
            scatter(raster{suIdx,trial},ones(size(raster{suIdx,trial},1),1)*trial,'.k')
            hold on 
        end
        
        ax=gca;
        if strcmp(condition,'Visual')
            patch([0 100 100 0],[0 0 size(raster,2) size(raster,2)],'y','FaceAlpha',.3,'EdgeColor','none')
            ax.XLim=[-1000, 5000];

        elseif strcmp(condition,'Optotagging')
            if size(raster,2)>=100
                patch([0 50 50 0],[0 0 size(raster,2) size(raster,2)],'c','FaceAlpha',.1,'EdgeColor','none')
                ax.XLim=[-100, 200];
            else
                patch([0 3000 3000 0],[0 0 size(raster,2) size(raster,2)],'c','FaceAlpha',.1,'EdgeColor','none')
                ax.XLim=[-2000, 4000];
            end
        end

        hold off
        ax.YLim=[0, size(raster,2)];
        ax.YTick=[];
        ax.YAxis.Color='none'; 
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
        figname=strcat(ID,'- P',int2str(age),' - ',condition,'- 1');
        sgtitle(figname)
        export_fig(fullfile(folder,figname),'-tiff','-transparent')
        close
    end
     
end