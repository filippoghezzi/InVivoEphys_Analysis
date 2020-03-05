clear 
close all
clc

data=readtable('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\V1_InVivo_SST;Ai32.csv'); 
spikeFolder='C:\Users\Butt Lab\Documents\SpikeSorting';
load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\A1x32_Map.mat');

data=data(data.Sorting~=0,:);
data=loadSampleDuration(data,spikeFolder);
data=data(data.Use~=0,:);

recToAnalyse=unique(data.MouseID);
su.PSTH=[];
su.layer=[];
su.age=[];
su.templates=[];
su.suid=[];
su.halfWidth=[];
su.troughPeakTime=[];
su.response=[];
su.optotagging=[];

% for recording=1:size(recToAnalyse,1)
for recording=19
    if ~isempty(data(strcmp(data.MouseID,recToAnalyse{recording}),:))
        spike=doAnalysis(data(strcmp(data.MouseID,recToAnalyse{recording}),:),ElectrodeMap,spikeFolder,1);
    end  
    su.PSTH=[su.PSTH;spike.PSTH];
    su.layer=[su.layer;spike.suLayer];
    su.age=[su.age;spike.age];
    su.templates=[su.templates;spike.templates];
    su.suid=[su.suid;spike.suid];
    su.halfWidth=[su.halfWidth;spike.halfWidth];
    su.troughPeakTime=[su.troughPeakTime;spike.troughPeakTime];
    su.response=[su.response;spike.response];
    su.optotagging=[su.optotagging;spike.optotagging];
end

% doClustering(allPSTH);


function spike=doAnalysis(tab,ElectrodeMap,spikeFolder,verbose)
    savingFolder=tab.Folder{1};
    [spike.spikeTimes,spike.templates,spike.suid]=LoadSpikes(fullfile(spikeFolder,tab.MouseID{1}),ElectrodeMap);
    
    if ~isempty(spike.suid)
        spike.age=ones(size(spike.suid,1),1)*tab.Age(1);
        spike.suLayer=findSingleUnitLayer(spike,savingFolder);
        [~,~,~,sr]=loadEventsForSpikes(fullfile(tab.Folder{1},tab.Experiment{1}));
        load(fullfile(savingFolder,'Stimuli.mat'))
        
        [spike.halfWidth,spike.troughPeakTime]=getTemplateFeatures(spike.templates,spike.suid,sr,verbose,savingFolder);
        
        %% Contralateral visual stimulation
        rasterWindow=[1,5];
        binSize=0.003;
        [spike.PSTH,PSTHbins,raster]=makePSTHandRaster(spike.spikeTimes,spike.suid,stimuli.Visual,sr,rasterWindow,binSize);
        spike.response=findVisualResponsiveUnits(raster);
        if verbose
            plotRaster(raster,spike,'Visual',savingFolder,tab.MouseID{1},tab.Age(1))
            plotPSTH(spike.PSTH,PSTHbins,spike,savingFolder,tab.MouseID{1},tab.Age(1))
        end
        
        %% SST Optotagging
        [~,~,rasterOptotagging]=makePSTHandRaster(spike.spikeTimes,spike.suid,stimuli.Optotagging,sr,[2,5],binSize);
        spike.optotagging=findVisualResponsiveUnits(rasterOptotagging);
        if verbose
            plotRaster(rasterOptotagging,spike,'Optotagging',savingFolder,tab.MouseID{1},tab.Age(1))
        end
        
        
        
    else
        spike.suid=[];
        spike.PSTH=[];
        spike.suLayer=[];
        spike.age=[];
        spike.response=[];
        spike.halfWidth=[];
        spike.troughPeakTime=[];
        spike.response=[];
        spike.optotagging=[];
        
    end
end
