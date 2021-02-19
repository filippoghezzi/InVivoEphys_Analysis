function stim=getStimuli(ops,data)

    if numel(ops.nSamplesBlocks) ~= height(data)
        error('Discrepancy between recordings and sample blocks size')
    end 
    
    stim.ledR=[];
    stim.laser=[];
    stim.laserAndLedR=[];
    stim.optotagging=[];
    stim.ledL=[];
    
    for file=1:height(data)
        if file==1
            startSamples=1;
            endSamples=ops.nSamplesBlocks(file);
        else
            startSamples=endSamples+1; 
            endSamples=startSamples+ops.nSamplesBlocks(file);
        end 
        
        [laser,ledR,laserAndLedR,ledL]=loadOpenEphysEvents(data.Experiment{file},startSamples);
       
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
        
        elseif strcmp(data.Protocol{file},'VisualFlash')
            stim.ledL=[stim.ledL;ledL];
            
        end
        
    end
    
end