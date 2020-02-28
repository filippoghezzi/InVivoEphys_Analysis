function [LFP]=loadLFPData(foldername,ElectrodeMap)


    %% Load data and downsample
    j=1;
    for i=max(ElectrodeMap):max(ElectrodeMap)
       [tmp,dataTime,~]=load_open_ephys_data_faster(strcat(foldername,'\100_CH',int2str(i),'.continuous'));
       tmp=tmp(1:20:end); %Downsampling to 1000 Hz
       data(j,:)=tmp';
       j=j+1;
    end

    dataTime=dataTime(1:20:end)'; %Downsampling to 1000 Hz
    data=data(ElectrodeMap-16,:); %Rearrange data according to electrode map

    %% Load events 
    [events,eventTime,~]=load_open_ephys_data_faster(strcat(foldername,'\all_channels.events'));
    eventArray=zeros(8,length(dataTime));
    
    if ~isempty(events)
        for j=1:8
            singleEventTime=eventTime(events==j-1); 
            if ~isempty(singleEventTime)
                % Build event logic array
                for i=1:2:length(singleEventTime)
                    [~,onset] = min(abs(dataTime-singleEventTime(i)));
                    [~,offset] = min(abs(dataTime-singleEventTime(i+1)));   
                    eventArray(j,onset:offset)=1;       
                end
            end
        end
    end
    
    %Output variable
    LFP.data=data;
    LFP.dataTime=dataTime;
    LFP.eventArray=eventArray;
    
    
    
end
