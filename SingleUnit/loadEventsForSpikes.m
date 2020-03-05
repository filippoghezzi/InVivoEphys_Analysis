function [stim1,stim2,stim12,sr]=loadEventsForSpikes(foldername,start)
% events: 0 -> laser
%         2 -> contralateral LED


% stim1 -> laser only
% stim2 -> LED only
% stim12 -> laser+ LED


 %% Load datatime
    [~,dataTime,info]=load_open_ephys_data_faster(strcat(foldername,'\100_CH17','.continuous'));
    sr=info.header.sampleRate;
    
    %% Load events 
    [events,eventTime,~]=load_open_ephys_data_faster(strcat(foldername,'\all_channels.events'));
    [events,eventTime]=runExceptions(foldername,events,eventTime);
    
    stim1=[];
    stim2=[];
    stim12=[];
    i=1;
    while i <= size(events,1)
        if events(i)==0 && events(i+1)==0 %Case laser only
            [~,onset] = min(abs(dataTime-eventTime(i)));
            stim1=[stim1;onset+start];
            i=i+2;
        
        elseif events(i)==2 && events(i+1)==2 %Case contralateral LED only
            [~,onset] = min(abs(dataTime-eventTime(i)));
            stim2=[stim2;onset+start];
            i=i+2;
            
        elseif events(i)==0 && events(i+1)==2 %Case LED + laser
            [~,onset] = min(abs(dataTime-eventTime(i+1)));
            stim12=[stim12;onset+start];
            i=i+4;   %%To fix in case multiple stimuli are presented at the same time.
        end
    end    
end

function [events,eventTime]=runExceptions(foldername,events,eventTime)

if strcmp(foldername,'D:\6.InVivo_V1_SST;Ai32\SC1\2019-03-29_14-21-43')
    events=[events(1:188);events(191:end)];
    eventTime=[eventTime(1:188);eventTime(191:end)];
elseif strcmp(foldername,'D:\6.InVivo_V1_SST;Ai32\SC14\2019-07-30_15-35-30')
    events=[events(1:40);events(43:end)];
    eventTime=[eventTime(1:40);eventTime(43:end)];    
elseif strcmp(foldername,'D:\6.InVivo_V1_SST;Ai32\SC15\2019-07-30_16-57-07')
    events=[];
    eventTime=[];
end
end