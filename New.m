clear 
close all
clc

data=readtable('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\V1_InVivo_SST;Ai32.csv'); 
spikeFolder='C:\Users\Butt Lab\Documents\SpikeSorting';
load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\A1x32_Map.mat');

data=data(data.Use~=0,:);
data=loadSampleDuration(data,spikeFolder);

recToAnalyse=unique(data.MouseID);
allPSTH=[];
for recording=1:length(recToAnalyse)
    if ~isempty(data(strcmp(data.MouseID,recToAnalyse{recording}),:))
        PSTH=doAnalysis(data(strcmp(data.MouseID,recToAnalyse{recording}),:),ElectrodeMap,spikeFolder);
    end  
    allPSTH=[allPSTH;PSTH];
end

doClustering(allPSTH);


function PSTH=doAnalysis(tab,ElectrodeMap,spikeFolder)
    savingFolder=tab.Folder{1};
    [spike.spikeTimes,spike.templates,spike.suid]=LoadSpikes(fullfile(spikeFolder,tab.MouseID{1}),ElectrodeMap);
    if ~isempty(spike.suid)
        spike.suLayer=findSingleUnitLayer(spike,savingFolder);
        [~,~,~,sr]=loadEventsForSpikes(fullfile(tab.Folder{1},tab.Experiment{1}));
        load(fullfile(savingFolder,'Stimuli.mat'))
        
        rasterWindow=[1,5];
        binSize=0.03;
        [PSTH,PSTHbins,raster]=makePSTHandRaster(spike.spikeTimes,spike.suid,stimuli.Visual,sr,rasterWindow,binSize);
%         plotRaster(raster)
%         plotPSTH(PSTH,PSTHbins)
%               findVisualResponsiveUnits(raster)   
    else
        PSTH=[];
    end
end
