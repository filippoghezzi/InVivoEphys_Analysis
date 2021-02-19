function [PSTH,PSTHbins,raster]=makePSTHandRaster(spiketimes,stimulus,fs,window,binSize,artefactRemoval)
% function [PSTH,PSTHbins,raster]=makePSTHandRaster(spiketimes,stimulus,fs,window,binSize,artefactRemoval)
% Obtain PSTH and raster from spiketimes aligned to stimulus given as
% input; PSTH window and bin size are given in inpu as well as windows to
% remove artefact due to light stimuli. 

        allTrialRaster = [];
        
        for trial=1:size(stimulus,1)
            tmpSu=double(spiketimes)-stimulus(trial,1);  % Single unit spike time (sample) aligned to stimulus
            
            if ~isempty(artefactRemoval)
                tmpSu=tmpSu(~(tmpSu>=artefactRemoval(1,1)*fs & tmpSu<=artefactRemoval(1,2)*fs));
                tmpSu=tmpSu(~(tmpSu>=artefactRemoval(2,1)*fs & tmpSu<=artefactRemoval(2,2)*fs));
            end
            
            raster(trial)={tmpSu(tmpSu>-window(1)*fs&tmpSu<window(2)*fs)};
            raster{trial}=raster{trial}.*(1000/fs);       % Convert from samples to time (ms)     
            allTrialRaster=[allTrialRaster;raster{trial}]; 
        end
        
        PSTHbins=-window(1)*1000:binSize*1000:window(2)*1000;
        PSTH=histcounts(allTrialRaster,PSTHbins);        
        PSTH=PSTH./binSize./size(stimulus,1);
        PSTH=smooth(PSTH,5);
    
        PSTHbins=PSTHbins(1:end-1);
end