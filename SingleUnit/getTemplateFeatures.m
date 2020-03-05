function [halfWidth,troughPeakTime]=getTemplateFeatures(template,suid,sr,varargin)

    if size(varargin,2)==1
        verbose = varargin{1,1};
    elseif size(varargin,2)==2
        verbose=varargin{1,1};
        folder=varargin{1,2};
    else
        verbose=0;
    end
    if verbose
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    
    for suIdx=1:size(suid,1)
        [~,trough,width,p]=findpeaks(-template(suIdx,:));
        [~,maxProminence]=max(p);
        trough=trough(maxProminence);
        width=width(maxProminence);
        
        [~,peakTroughSample]=max(template(suIdx,trough:end));
        
        if peakTroughSample+trough>size(template,2)
            peakTroughSample=peakTroughSample-1;
        end
        
        halfWidth(suIdx,1)=width/sr*1000;
        troughPeakTime(suIdx,1)=peakTroughSample/sr*1000;
        


        if verbose
            subplot(ceil(sqrt(size(suid,1))),ceil(sqrt(size(suid,1))),suIdx)
            hold on
            plot(template(suIdx,:))
            scatter(trough,template(suIdx,trough),'r*')
            scatter(peakTroughSample+trough,template(suIdx,peakTroughSample+trough),'r*')
            title(strcat('Unit-',int2str(suid(suIdx))))
            hold off 
        end
    end 
    
    if verbose
        figname=fullfile(folder,'Templates');
        export_fig(figname,'-tiff','-transparent')
        close 
    end
end