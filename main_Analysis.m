clear 
close all
clc

%% Load data
foldername='Z:\In vivo data\2019-01-16_12-53-00';
load(strcat(foldername,'\Data.mat'))

%% Current source density
stimulusIdx=3; %1 for blue laser; 3 for LED visual stimulation.

stimulusOnsetIdx=find([0,diff(eventArray(stimulusIdx,:))]>0);
stimulusOffsetIdx=find([0,diff(eventArray(stimulusIdx,:))]<0);
stimuliIdx=[stimulusOnsetIdx',stimulusOffsetIdx'];

[CSD,CSDinfo]=getAverageCSD(data,[stimuliIdx(:,1) stimuliIdx(:,2)+500],500,800,25,1,1,1);

%%
factor = 1; % the time factor (20 or 1 per msecond depending on sampling (MAKE IT AUTOMATIC ACCORDING TO DOWNSAMPLING)

evokedLFP = [];
i=0;
for channel=1:size(data,1)
    if any(channel==CSDinfo.L4)
        i=i+1;
        for trial=1:length(stimulusOnsetIdx)
            LFPwindow = (stimulusOnsetIdx(trial)-1000*factor:stimulusOnsetIdx(trial)+100+5000);    
            evokedLFP(i,trial,1:length(LFPwindow))=data(channel,LFPwindow);
        end
    end
end 

% save('Data.mat','stimTrials','stimTraces','stimT','data','ElectrodeMap','-v7.3');
PlotAllData(squeeze(mean(evokedLFP(:,:,:),2)));

timeFrequencyAnalysis(evokedLFP)
