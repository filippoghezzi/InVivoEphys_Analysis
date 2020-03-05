clear 
close all
clc

tab=readtable('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\V1_InVivo_SST;Ai32.csv'); 
load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\A1x32_Map.mat');
tab=tab(tab.Use~=0,:);
spikeFolder='C:\Users\Butt Lab\Documents\SpikeSorting';

tab=loadSampleDuration(tab,spikeFolder);



for i=20:max(unique(tab.RecID))
    if ~isempty(tab(tab.RecID==i,:))
        CSDinfo = Analysis_SingleAnimal(tab(tab.RecID==i,:),ElectrodeMap,spikeFolder);
        save(fullfile(tab(tab.RecID==i,:).Folder{1},'CSDinfo.mat'),'CSDinfo')
    end
end


%% Main data analysis
function CSDinfo = Analysis_SingleAnimal(tab,ElectrodeMap,spikeFolder)
    savingFolder=tab.Folder{1};
    
    optotaggingSession=1;
    visualSession=1;
    baselineSession=1;
    tab=sortrows(tab,'Use');
    [spikes,templates,suid]=LoadSpikes(fullfile(spikeFolder,tab.MouseID{1}),ElectrodeMap);

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
        eventIdx=loadEventsForSpikes(dataFolder);

        %% Choose analysis
        switch tab.Protocol{i}
            case 'VisualFlash'
                CSDinfo=eLFP_Analysis(LFP.data,LFP.dataTime,LFP.eventArray,0,[],'LED', spikes(ismember(spikes(:,2),[1,2]),:), eventIdx,tab.startSample(i));
                figname=strcat(tab.MouseID{i},'- P',int2str(tab.Age(i)),' - Visual Flash LFP - Laser OFF - 1');
%                 sgtitle(figname)
%                 export_fig(fullfile(savingFolder,figname),'-tiff','-transparent')
                close
                visualSession=visualSession+1;
            case 'VisualFlash_Opto'
                CSDinfo=eLFP_Analysis(LFP.data,LFP.dataTime,LFP.eventArray,0,[],'LED', spikes(ismember(spikes(:,2),[1,2]),:),eventIdx,tab.startSample(i));
                figname=strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - Visual Flash LFP - Laser OFF');
%                 sgtitle(figname)
%                 export_fig(fullfile(savingFolder,figname),'-tiff','-transparent')
                close
                eLFP_Analysis(LFP.data,LFP.dataTime,LFP.eventArray,1,CSDinfo,'LED', spikes(ismember(spikes(:,2),[1,2]),:),eventIdx,tab.startSample(i));
                figname=strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - Visual Flash LFP - Laser ON');
%                 sgtitle(figname)
%                 export_fig(fullfile(savingFolder,figname),'-tiff','-transparent')
                close
            case 'Baseline'
%                 Baseline_Analysis(LFP.data,LFP.dataTime,LFP.eventArray)
            case 'Baseline_Opto'
                eLFP_Analysis(LFP.data,LFP.dataTime,LFP.eventArray,1,CSDinfo,'Laser', spikes(ismember(spikes(:,2),[1,2]),:),eventIdx,tab.startSample(i));
                figname=strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - LFP Laser only');
%                 sgtitle(figname)
%                 export_fig(fullfile(savingFolder,figname),'-tiff','-transparent')
                close
                baselineSession =baselineSession +1;
            case 'Optotagging'
                figname=strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - Optotagging',int2str(optotaggingSession));
%                 sgtitle(figname)
%                 export_fig(fullfile(savingFolder,figname),'-tiff','-transparent')
                close
                optotaggingSession=optotaggingSession+1;
            case 'ComplexSession'
                
        end

    end
    

end

%% Analysis functions

function CSDinfo=eLFP_Analysis(data,dataTime,eventArray,opto,CSDinfo,alignTrials, mua,eventIdx,start)
    %Event for spikes, used for MUA
    %% Find events with laser ON or OFF
    if ~isempty(eventIdx{1})
        if ~isempty(eventIdx{3}) %Case both laser and LED
            for i=1:size(eventIdx{3},1)
                LaserON(i,1)=any(eventIdx{3}(i,1)-eventIdx{1}(:,1)==1000);
            end
            if opto
                stimulus=eventIdx{3}(LaserON,:)+start;
            else
                stimulus=eventIdx{3}(~LaserON,:)+start;
            end
        else %Case laser only
            stimulus=eventIdx{1}+start; %3 for LED
        end
    else %Case LED only
        stimulus=eventIdx{3}+start; %3 for LED
    end
    

    %% Current source density
    LED=3;          %3, index of eventArray line for LED visual stimulation;
    laser=1;        %1, index of eventArray line for 470nm laser;
    
    if strcmp(alignTrials,'LED')
        % Align on LED onset
        LEDOnsetIdx=find([0,diff(eventArray(LED,:))]>0);
        LEDOffsetIdx=find([0,diff(eventArray(LED,:))]<0)-1;
        
        %Select trial in which laser is OFF or ON
        LEDOnsetIdx=LEDOnsetIdx(eventArray(laser,LEDOnsetIdx)==opto);
        LEDOffsetIdx=LEDOffsetIdx(eventArray(laser,LEDOffsetIdx)==opto);
        stimuliIdx=[LEDOnsetIdx',LEDOffsetIdx'];
        
    elseif strcmp(alignTrials,'Laser')
        % Align on Laser
        LaserOnsetIdx=find([0,diff(eventArray(laser,:))]>0); 
        LaserOffsetIdx=find([0,diff(eventArray(laser,:))]<0)-1;
        stimuliIdx=[LaserOnsetIdx',LaserOffsetIdx'];
        
    end

    if opto
        [CSD,~]=getAverageCSD(data,[stimuliIdx(:,1) stimuliIdx(:,2)+5100],1000,800,25,1,1,0, mua,stimulus);
    else
        [CSD,CSDinfo]=getAverageCSD(data,[stimuliIdx(:,1) stimuliIdx(:,2)+5100],1000,800,25,1,1,1, mua,stimulus);
    end

    %% Evoked LFP
    factor = 1; % the time factor (20 or 1 per msecond depending on sampling (MAKE IT AUTOMATIC ACCORDING TO DOWNSAMPLING)

    evokedLFP = [];
    i=0;
    for channel=1:size(data,1)
        if any(channel==CSDinfo.L4)
            i=i+1;
            for trial=1:size(stimuliIdx,1)
                if ~(stimuliIdx(trial,1)+100+5000 >size(data,2))
                    LFPwindow=(stimuliIdx(trial,1)-1000*factor:stimuliIdx(trial,1)+100+5000);    
                    evokedLFP(i,trial,1:length(LFPwindow))=data(channel,LFPwindow);
                end
            end
        end
    end 

    %% Time-frequency plot and all plots
    figure('units','normalized','outerposition',[0 0 1 1]);
    meanLFP=mean(squeeze(mean(evokedLFP(:,:,:),1)),1);
    meanSweep=squeeze(mean(evokedLFP,1));
    subplot(4,1,1)
    hold on
    for sweep=1:size(evokedLFP,2)
        PlotAllData(meanSweep(sweep,:),1000,'Grey')
    end
    PlotAllData(meanLFP,1000,'k');
    ax=gca;
    ax.XLim=[-1000, 2000];
    ax.YLim=[min(meanLFP,[],'all')-20 max(meanLFP,[],'all')+20];
    ax.YLabel.String='Voltage (\muV)';
    ax.YLabel.FontSize=15;
    ax.XAxis.Color='none';
    
    if strcmp(alignTrials,'LED')
        patch([0 100 100 0],[min(meanLFP)-20 min(meanLFP)-20 max(meanLFP)+20 max(meanLFP)+20],'y','FaceAlpha',.2,'EdgeColor','none')
    end
    if opto
        if strcmp(alignTrials,'Laser')
            patch([0 200 200 0],[min(meanLFP)-20 min(meanLFP)-20 max(meanLFP)+20 max(meanLFP)+20],'c','FaceAlpha',.1,'EdgeColor','none')
        else
            patch([-50 150 150 -50],[min(meanLFP)-20 min(meanLFP)-20 max(meanLFP)+20 max(meanLFP)+20],'c','FaceAlpha',.1,'EdgeColor','none')
        end
    end
    box off
    
    
    subplot(4,1,4)
    imagesc([-1000 5000],[0 750],CSD,[min(CSD(:))/5,max(CSD(:))/5])
    colormap('jet')
    ax=gca;
    ax.XLim=[-1000, 2000];
    ax.YLabel.String='Cortical depth (\mum)';
    ax.YLabel.FontSize=15;
    ax.XLabel.String='Time (ms)';
    ax.XLabel.FontSize=15;
    box off
   
    subplot(4,1,2)
    timeFrequencyAnalysis(squeeze(mean(evokedLFP,1)))
    colorbar('off')
    ax=gca;
    ax.XLim=[-1000, 2000];
    ax.XTick=[];
    ax.YLabel.String='Frequency (Hz)';
    ax.YLabel.FontSize=15;
    ax.XAxis.Color='none'; 
    box off


end