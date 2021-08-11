function [laser,ledR,laserAndLedR,ledL]=loadOpenEphysEvents(foldername,start)
% events: 0 -> laser
%         2 -> contralateral LED
    
    addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\analysis-tools')) % path to kilosort folder

    %% Load datatime and events
    if exist(fullfile(foldername,'100_CH17_2.continuous'),'file')
        [~,dataTime,~]=load_open_ephys_data_faster(strcat(foldername,'\100_CH17_2','.continuous'));
        [events,eventTime,~]=load_open_ephys_data_faster(strcat(foldername,'\all_channels_2.events'));
    elseif exist(fullfile(foldername,'100_CH17.continuous'),'file')  %%%% This condition account for problem in pausing then restarting the recording. 
        [~,dataTime,~]=load_open_ephys_data_faster(strcat(foldername,'\100_CH17','.continuous'));
        [events,eventTime,~]=load_open_ephys_data_faster(strcat(foldername,'\all_channels.events'));
    elseif exist(fullfile(foldername,'107_CH17.continuous'),'file')
        [~,dataTime,~]=load_open_ephys_data_faster(strcat(foldername,'\107_CH17','.continuous'));
        [events,eventTime,~]=load_open_ephys_data_faster(strcat(foldername,'\all_channels.events'));
    end
    laser=[];
    ledR=[];
    laserAndLedR=[];
    ledL=[];
    i=1;
    while i <= size(events,1)
        if events(i)==0 && events(i+1)==0 %Case laser only
            [~,onset] = min(abs(dataTime-eventTime(i)));
            [~,offset] = min(abs(dataTime-eventTime(i+1)));
            if ~(offset-onset > 30000) % Discard all 3-s long optotagging ramp 
                laser=[laser;onset+start];
            end
            i=i+2;
        
        elseif events(i)==2 && events(i+1)==2 %Case contralateral LED only
            [~,onset] = min(abs(dataTime-eventTime(i)));
            ledR=[ledR;onset+start];
            i=i+2;
            
        elseif events(i)==0 && events(i+1)==2 %Case LED + laser
            if (eventTime(i+1)-eventTime(i))>0.0510
                error('Delay between laser and LED bigger than 50 ms; check events'); 
            end
            [~,onset] = min(abs(dataTime-eventTime(i+1)));
            laserAndLedR=[laserAndLedR;onset+start];
            i=i+4;   %%To fix in case multiple stimuli are presented at the same time.
        
        elseif events(i)==6 && events(i+1)==6 %Case ipsilateral LED only
            [~,onset] = min(abs(dataTime-eventTime(i+1)));
            ledL=[ledL;onset+start];
            i=i+2;

        else
            i=i+2;
        end
        

    end    
end

