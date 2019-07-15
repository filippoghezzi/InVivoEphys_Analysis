%%%%% Analysis to add: - effects of blue light on LFP and MUA activity.
%                      - single units: rasters, responsive, percentage of responsive neurons across development
%                      - repeat across cortical depth: see single unit according to depth



clear 
close all
clc

load('InVivo_V1_SST;Ai32_Home.mat')
load('A1x32_Map.mat');
spikeFolder='C:\Users\Filippo\Desktop\INVIVO\SpikeSorting';

%% Load sample duration
tab.SampleLength(1,1)=0;
tab.startSample(1,1)=0;
tab.endSample(1,1)=0;
for i=1:unique(tab.RecID)
    load(fullfile(spikeFolder,tab(tab.RecID==i,:).MouseID{1},'SampleDuration.mat'))
    tab(tab.RecID==i,:).SampleLength=samplesToSave;
    for j=1:height(tab(tab.RecID==i,:))
        if j==1
            tab.startSample(j)=1;
            tab.endSample(j)=tab(tab.RecID==i,:).SampleLength(j);
        else
            tab.startSample(j)=tab(tab.RecID==i,:).endSample(j-1)+1;
            tab.endSample(j)=tab(tab.RecID==i,:).SampleLength(j)+tab.startSample(j);
        end
    end
    
end
tab = removevars(tab,'SampleLength');

%% Main data analysis
for i=1:height(tab)
    dataFolder=fullfile(tab.Folder{i},tab.Experiment{i});
    dataFolderContent=dir(dataFolder);
    
    %% Load LFP data
    dataPresence=[];
    for j=1:length(dataFolderContent)
        dataPresence=[dataPresence;strcmp('LFP_Data.mat',dataFolderContent(j).name)];
    end
    if ~any(dataPresence)
        [LFP]=loadLFPData(dataFolder,ElectrodeMap);
        save(fullfile(dataFolder,'LFP_Data.mat'),'LFP')
    else
        load(fullfile(dataFolder,'LFP_Data.mat'))
    end
    
    %% Load spike data
    [spikes,templates,suid]=LoadSpikes(fullfile(spikeFolder,tab.MouseID{i}),ElectrodeMap);
    eventIdx=loadEventsForSpikes(dataFolder);
    
    %% Choose analysis
    switch tab.Protocol{i}
        case 'VisualFlash'
%             CSDinfo=eLFP_Analysis(LFP.data,LFP.dataTime,LFP.eventArray,0,[]);
% %             Spike_Analysis(spikes,templates,suid,eventIdx,CSDinfo.L4,tab.startSample(i))
%             sgtitle(strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - Laser:OFF'))
        case 'VisualFlash_Opto'
            CSDinfo=eLFP_Analysis(LFP.data,LFP.dataTime,LFP.eventArray,0,[]);
%             Spike_Analysis(spikes,templates,suid,eventIdx,CSDinfo.L4,tab.startSample(i))
            sgtitle(strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - Laser:OFF'))
            eLFP_Analysis(LFP.data,LFP.dataTime,LFP.eventArray,1,CSDinfo);
            sgtitle(strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - Laser:ON'))
        case 'Baseline'
            Baseline_Analysis(LFP.data,LFP.dataTime,LFP.eventArray)
    end
    
end

%% Analysis functions
function Spike_Analysis(spikes,templates,suid,eventIdx,L4,start)
%% Plot multi unit activity histogram
    mua=spikes(ismember(spikes(:,2),[1,2]),:); % Select both single units and multi units
    mua=mua(ismember(mua(:,3),L4),:); % Only in Layer 4
    
    LED=eventIdx{3}+start; %3 for LED
    muaTimes=[];
    for i=1:size(LED,1)
        tmpSpikes(:,4)=mua(:,4)-LED(i);
        muaTimes=[muaTimes;tmpSpikes(tmpSpikes(:,4)>-20000&tmpSpikes(:,4)<5*20000,4)];
    end
    muaTimes=double(muaTimes)./20; % Convert to double and divide by 20 to get time in ms
 
    subplot(4,1,4)
    histogram(muaTimes,600)
    xlim([-1000 5000])
    
%% Plot single units
    for i=1:size(suid,1)
        su=spikes(ismember(spikes(:,1),suid(i)),:);
        for j=1:size(LED,1)
            tmpSu=su(:,4)-LED(j);
            raster={tmpSu(tmpSu>-20000&tmpSu<5*20000)};
        end
% % % % %         plot rasters
    end
end

function CSDinfo=eLFP_Analysis(data,dataTime,eventArray,opto,CSDinfo)
    %% Current source density
    LED=3;          %3, index of eventArray line for LED visual stimulation;
    laser=1;        %1, index of eventArray line for 470nm laser;
    
    % Align on LED onset
    LEDOnsetIdx=find([0,diff(eventArray(LED,:))]>0);
    LEDOffsetIdx=find([0,diff(eventArray(LED,:))]<0)-1;

	% % Align on Laser
    % stimulusOnsetIdx=find([0,diff(eventArray(laser,:))]>0); 
    % stimulusOffsetIdx=find([0,diff(eventArray(laser,:))]<0)-1;
    
    %Select trial in which laser is OFF or ON

    LEDOnsetIdx=LEDOnsetIdx(eventArray(laser,LEDOnsetIdx)==opto);
    LEDOffsetIdx=LEDOffsetIdx(eventArray(laser,LEDOffsetIdx)==opto);
    stimuliIdx=[LEDOnsetIdx',LEDOffsetIdx'];
    
    if opto
        [CSD,~]=getAverageCSD(data,[stimuliIdx(:,1) stimuliIdx(:,2)+5100],1000,800,25,1,1,0);
    else
        [CSD,CSDinfo]=getAverageCSD(data,[stimuliIdx(:,1) stimuliIdx(:,2)+5100],1000,800,25,1,1,1);
    end

    %% Evoked LFP
    factor = 1; % the time factor (20 or 1 per msecond depending on sampling (MAKE IT AUTOMATIC ACCORDING TO DOWNSAMPLING)

    evokedLFP = [];
    i=0;
    for channel=1:size(data,1)
        if any(channel==CSDinfo.L4)
            i=i+1;
            for trial=1:size(stimuliIdx,1)
                LFPwindow=(LEDOnsetIdx(trial)-1000*factor:LEDOnsetIdx(trial)+100+5000);    
                evokedLFP(i,trial,1:length(LFPwindow))=data(channel,LFPwindow);
            end
        end
    end 

    %% Time-frequency plot and all plots
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(4,1,1)
    PlotAllData(mean(squeeze(mean(evokedLFP(:,:,:),1)),1));
    xlim([-0, 6000])
    subplot(4,1,2)
    imagesc(CSD)
    xlim([-0000, 6000])
    subplot(4,1,3)
    timeFrequencyAnalysis(squeeze(mean(evokedLFP,1)))
    xlim([-1000, 5000])
    colorbar('off')


end

function Baseline_Analysis(data,dataTime,eventArray)

end