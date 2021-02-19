clear 
clc
close all

dataFolder='C:\Users\Butt Lab\Documents\SpikeSorting';

data=readtable('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\V1_InVivo_SST;Ai32.csv'); 
data=data(data.Sorting~=0,:);


rec='SC37';

Stimuli_Main(fullfile(dataFolder,rec),data(strcmp(data.MouseID,rec),:));



function Stimuli_Main(dir,data)
    load(fullfile(dir,'rez.mat'),'rez')
  
    samples=rez.ops.nSamplesBlocks;
    
    if numel(samples) ~= height(data)
        error('Discrepancy between recordings and sample blocks size')
    end
    
    for file=1:height(data)
        if file==1
            data.startSamples(file)=1;
            data.endSamples(file)=samples(file);
        else
            data.startSamples(file)=data.endSamples(file-1)+1;
            data.endSamples(file)=data.startSamples(file)+samples(file);
        end     
    end
    
    
    stim.ledR=[];
    stim.laser=[];
    stim.laserAndLedR=[];
    stim.optotagging=[];
    for file=1:height(data)
        [laser,ledR,laserAndLedR]=loadOpenEphysEvents(fullfile(data.Folder{file},data.Experiment{file}),data.startSamples(file));
       
        if strcmp(data.Protocol{file},'VisualFlash_Opto')
            stim.ledR=[stim.ledR;ledR];
            stim.laser=[stim.laser;laser];
            stim.laserAndLedR=[stim.laserAndLedR;laserAndLedR];
            
        elseif strcmp(data.Protocol{file},'Baseline_Opto')
            stim.laser=[stim.laser;laser];
            
        elseif strcmp(data.Protocol{file},'Optotagging')
            stim.optotagging=[stim.optotagging;laser];
            
        elseif strcmp(data.Protocol{file},'ComplexSession')
            stim.ledR=[stim.ledR;ledR];
            stim.optotagging=[stim.optotagging;laser];
            stim.laserAndLedR=[stim.laserAndLedR;laserAndLedR];                
        end
        
        
    end
    
    save(fullfile(dir,'stimuli.mat'),'stim')
end
   