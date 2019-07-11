function outputEvents=loadEventsForSpikes(foldername)

 %% Load datatime
    [~,dataTime,~]=load_open_ephys_data_faster(strcat(foldername,'\100_CH17','.continuous'));


    %% Load events 
    [events,eventTime,~]=load_open_ephys_data_faster(strcat(foldername,'\all_channels.events'));
    eventArray=zeros(8,length(dataTime));
    
    if isempty(events)
        outputEvents=[];
    else
        for j=1:8
            eventIdx=[];
            singleEventTime=eventTime(events==j-1); 
            if ~isempty(singleEventTime)
                % Build event logic array
                for i=1:2:length(singleEventTime)
                    [~,onset] = min(abs(dataTime-singleEventTime(i)));
                    [~,offset] = min(abs(dataTime-singleEventTime(i+1)));   
                    eventIdx=[eventIdx;onset,offset];    
                end
            end
            outputEvents(j)={eventIdx};
        end
    end
    
    %Output variable
    
    
    
    
end