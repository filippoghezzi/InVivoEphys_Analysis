clear 
close all
clc

%% Load data
foldername='C:\Users\Butt Lab\Desktop\In vivo data\2019-01-16_12-53-00';
load(strcat(foldername,'\Data.mat'))

%% Current source density
stimulusIdx=3; %1 for blue laser; 3 for LED visual stimulation.

stimulusOnsetIdx=find([0,diff(eventArray(stimulusIdx,:))]>0);
stimulusOffsetIdx=find([0,diff(eventArray(stimulusIdx,:))]<0);
stimuli=[stimulusOnsetIdx',stimulusOffsetIdx'];

[avgCsd1, cdInfo1] = CalculateAverageCSD(data,[stimuli(:,1) stimuli(:,2)+500],500,800,25,1,1);
imagesc(avgCsd1)

% Determine L4
[~,L4Idx]=max(cdInfo1.Sink(:,1));
L4Idx=[L4Idx*2-1,L4Idx*2];
%%
factor = 1; % the time factor (20 or 1 per msecond depending on sampling (MAKE IT AUTOMATIC ACCORDING TO DOWNSAMPLING)

evokedLFP = [];
i=0;
for channel=1:size(data,1)
    if any(channel==L4Idx)
        i=i+1;
        for trial=1:length(stimulusOnsetIdx)
            LFPwindow = (stimulusOnsetIdx(trial)-1000*factor:stimulusOnsetIdx(trial)+100+5000);    
            evockedLFP(i,trial,1:length(LFPwindow))=data(channel,LFPwindow);
        end
    end
end 

% save('Data.mat','stimTrials','stimTraces','stimT','data','ElectrodeMap','-v7.3');
PlotAllData(squeeze(mean(evokedLFP(:,:,:),2)));

timeFrequencyAnalysis(evokedLFP)
