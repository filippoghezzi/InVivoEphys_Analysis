clear 
close all
clc
plotUnits=1;
foldername='D:\6.InVivo_V1_SST;Ai32\SC28\2020-02-24_13-20-31';

load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\A1x32_Map.mat');
[spike.spikeTimes,spike.templates,spike.suid]=LoadSpikes('C:\Users\Butt Lab\Documents\SpikeSorting\SC28',ElectrodeMap);
start=55137281;
spike.spikeTimes(:,4)=spike.spikeTimes(:,4)-start;

    %% Load data and downsample
    j=1;
    
    for i=min(ElectrodeMap):max(ElectrodeMap)
       [tmp,dataTime,info]=load_open_ephys_data_faster(strcat(foldername,'\100_CH',int2str(i),'.continuous'));
       sr=info.header.sampleRate;

       [b,a]=butter(5,400/(sr/2),'high');
       tmp=filtfilt(b,a,tmp);
       
       data(j,:)=tmp';
       j=j+1;
    end
    
    data=data(ElectrodeMap-16,:); %Rearrange data according to electrode map

    stop=(max(dataTime)-min(dataTime))*sr;
    %% Load events 
    [events,eventTime,~]=load_open_ephys_data_faster(strcat(foldername,'\all_channels.events'));
    j=1;
    for i = 1:2:size(events,1)
        onset(j,1)=(eventTime(i)-min(dataTime))*sr;
        offset(j,1)=(eventTime(i+1)-min(dataTime))*sr;
        j=j+1;
    end

    spike.spikeTimes=spike.spikeTimes(spike.spikeTimes(:,2)==2,:); %Only single unit
    
    for i=1:size(data,1)
        hold on
%         plot(dataTime',data(i,:)-(i*200),'k')
        plot(data(i,:)-(i*200),'k')
%         xlim([min(dataTime) max(dataTime)])
        if plotUnits
            su=spike.spikeTimes(ismember(spike.spikeTimes(:,3),i),:);  %Select channel
            suid=unique(su(:,1));
            if ~isempty(suid)
                for j=1:size(suid,1)
                    su=su(ismember(su(:,1),suid(j,1)),:);
                    su=su(su(:,4)>0 & su(:,4)<stop-10000,:);
                    scatter(su(:,4),data(i,su(:,4))-(i*200),'r*')
                end
            end
        end
    end
    
    for i=1:size(onset,1)
            patch([onset(i) offset(i) offset(i) onset(i)],[0 0 -7000 -7000],'c','FaceAlpha',.3,'EdgeColor','none')
    end
    
   
    

