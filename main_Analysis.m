clear 
close all
clc

%% Load data
foldername='C:\Users\Butt Lab\Desktop\In vivo data\8.Lpar1.4\2019-02-06_12-39-02';
load(strcat(foldername,'\Data.mat'))

%% Current source density
LED=3; %1 for blue laser; 3 for LED visual stimulation.
laser=1;

stimulusOnsetIdx=find([0,diff(eventArray(LED,:))]>0);
stimulusOffsetIdx=find([0,diff(eventArray(LED,:))]<0)-1;

% Align on Laser
% stimulusOnsetIdx=find([0,diff(eventArray(laser,:))]>0); 
% stimulusOffsetIdx=find([0,diff(eventArray(laser,:))]<0)-1;

stimulusOnsetIdx=stimulusOnsetIdx(eventArray(laser,stimulusOnsetIdx)==0); %Only the trial in which laser was off. 
stimulusOffsetIdx=stimulusOffsetIdx(eventArray(laser,stimulusOffsetIdx)==0); %Only the trial in which laser was off. 

stimuliIdx=[stimulusOnsetIdx',stimulusOffsetIdx'];

[CSD,CSDinfo]=getAverageCSD(data,[stimuliIdx(:,1) stimuliIdx(:,2)+5100],1000,800,25,1,1,1);

%%
factor = 1; % the time factor (20 or 1 per msecond depending on sampling (MAKE IT AUTOMATIC ACCORDING TO DOWNSAMPLING)

evokedLFP = [];
i=0;
for channel=1:size(data,1)
    if any(channel==CSDinfo.L4)
        i=i+1;
        for trial=1:size(stimuliIdx,1)
            LFPwindow=(stimulusOnsetIdx(trial)-1000*factor:stimulusOnsetIdx(trial)+100+5000);    
            evokedLFP(i,trial,1:length(LFPwindow))=data(channel,LFPwindow);
        end
    end
end 

% save('Data.mat','stimTrials','stimTraces','stimT','data','ElectrodeMap','-v7.3');
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(3,1,1)
PlotAllData(mean(squeeze(mean(evokedLFP(:,:,:),1)),1));
subplot(3,1,2)
imagesc(CSD)
subplot(3,1,3)
timeFrequencyAnalysis(squeeze(mean(evokedLFP,1)))
