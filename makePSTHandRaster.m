function [PSTH,PSTHbins,raster]=makePSTHandRaster(spikes,suid,stimulus,sr,window,binSize)

    for suIdx=1:size(suid,1)
        unit=spikes(ismember(spikes(:,1),suid(suIdx)),:);
        allTrialRaster = [];
        
        for trial=1:size(stimulus,1)
            tmpSu=unit(:,4)-stimulus(trial,1);  % Single unit spike time (sample) aligned to stimulus
            raster(suIdx,trial)={tmpSu(tmpSu>-window(1)*sr&tmpSu<window(2)*sr)};
            raster{suIdx,trial}=double(raster{suIdx,trial}).*(1000/sr);       % Convert from samples to time (ms)     
            allTrialRaster=[allTrialRaster;raster{suIdx,trial}]; 
        end
        
        PSTHbins=-window(1)*1000:binSize*1000:window(2)*1000;
        PSTH(suIdx,:)=histcounts(allTrialRaster,PSTHbins);        
        PSTH(suIdx,:)=PSTH(suIdx,:)./binSize./size(stimulus,1);
        PSTH(suIdx,:)=smooth(PSTH(suIdx,:),5);
    end
    PSTHbins=PSTHbins(1:end-1);
end