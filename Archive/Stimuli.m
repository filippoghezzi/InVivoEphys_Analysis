clear 
close all
clc
% Stimuli values are samples from start of first file (Baseline)

data=readtable('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\V1_InVivo_SST;Ai32.csv'); 
spikeFolder='C:\Users\Butt Lab\Documents\SpikeSorting';
data=data(data.Sorting~=0,:);
data=loadSampleDuration(data,spikeFolder);

recToAnalyse=unique(data.MouseID);
for recording=20:20
    subset=data(strcmp(data.MouseID,recToAnalyse{recording}),:);
    if ~isempty(subset)
        savingFolder=subset.Folder{1};
        led=[];
        laser=[];
        ledAndLaser=[];
        optotagging=[];
        
        for i=1:height(subset)
            subset.Protocol{i}
            [stim1,stim2,stim12,sr]=loadEventsForSpikes(fullfile(subset.Folder{i},subset.Experiment{i}),subset.startSample(i));             
            if strcmp(subset.Protocol{i},'VisualFlash_Opto')
                led=[led;stim2];
                laser=[laser;stim1];
                ledAndLaser=[ledAndLaser;stim12];
            elseif strcmp(subset.Protocol{i},'Baseline_Opto')
                laser=[laser;stim1];
            elseif strcmp(subset.Protocol{i},'Optotagging')
                optotagging=[optotagging;stim1];
            elseif strcmp(subset.Protocol{i},'ComplexSession')
                led=[led;stim2];
                optotagging=[laser;stim1];
                ledAndLaser=[ledAndLaser;stim12];                
            end
        end  
    end
    
    
    stimuli.Visual=led;
    stimuli.Laser=laser;
    stimuli.VisualOpto=ledAndLaser;
    
    stimuli.Optotagging=optotagging;
    save(fullfile(savingFolder,'Stimuli.mat'),'stimuli');
end
