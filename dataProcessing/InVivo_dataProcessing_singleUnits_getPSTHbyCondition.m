function s = InVivo_dataProcessing_singleUnits_getPSTHbyCondition(s,stim,ops,verbose) 

    rasterWindow=ops.LFPwindow; %s
    binSize=ops.PSTHbinSize; %s
    
    if ops.SalB
        if strcmp(ops.brainArea, 'V1')
            conditions={'Visual','Visual_K'};
        elseif strcmp(ops.brainArea, 'S1BF')
            conditions={'WhiskerStim','WhiskerStim_K'};
        end
    else
        conditions={'Visual','Optotagging','VisualOpto','LaserOnly'};
    end

    for cond=1:numel(conditions)
        switch conditions{cond}
            case 'Visual'
                stimulus=stim.ledR;
                artefactRemoval=[];
                [s.PSTHvisual,s.PSTHbins,raster]=getPSTHbyUnit(s,stimulus,ops.fs,rasterWindow,binSize,artefactRemoval);
                s.response.visual=getUnitResponsiveness(raster,conditions{cond});
                if verbose; plotPSTH(s,conditions{cond},ops.dirOUT); end
            
            case 'Visual_K'
                stimulus=stim.ledR_SalB;
                artefactRemoval=[];
                [s.PSTHvisual_K,s.PSTHbins,raster]=getPSTHbyUnit(s,stimulus,ops.fs,rasterWindow,binSize,artefactRemoval);
                s.response.visual_K=getUnitResponsiveness(raster,conditions{cond});
                if verbose; plotPSTH(s,conditions{cond},ops.dirOUT); end
            
            case 'WhiskerStim'
                stimulus=stim.whiskerStim(:,1);
                artefactRemoval=[];
                [s.PSTHwhisker,s.PSTHbins,raster]=getPSTHbyUnit(s,stimulus,ops.fs,rasterWindow,binSize,artefactRemoval);
                s.response.whisker=getUnitResponsiveness(raster,conditions{cond});
                if verbose; plotPSTH(s,conditions{cond},ops.dirOUT); end
            
            case 'WhiskerStim_K'
                stimulus=stim.whiskerStim_SalB(:,1);
                artefactRemoval=[];
                [s.PSTHwhisker_K,s.PSTHbins,raster]=getPSTHbyUnit(s,stimulus,ops.fs,rasterWindow,binSize,artefactRemoval);
                s.response.whisker_K=getUnitResponsiveness(raster,conditions{cond});
                if verbose; plotPSTH(s,conditions{cond},ops.dirOUT); end
                
            case 'Optotagging'
                stimulus=stim.optotagging;
                artefactRemoval=[-1,1;49,53]*10^-3; %s
                [s.PSTHoptotagging,s.PSTHbins,raster]=getPSTHbyUnit(s,stimulus,ops.fs,rasterWindow,binSize,artefactRemoval);
                s.response.optotagging=getUnitResponsiveness(raster,conditions{cond});
                if verbose; plotPSTH(s,conditions{cond},ops.dirOUT); end

            case 'VisualOpto'
                if ~isempty(stim.laserAndLedR)
                    stimulus=stim.laserAndLedR;
                    artefactRemoval=[-51,49;149,153]*10^-3; %s
                    [s.PSTHvisualOpto,s.PSTHbins,raster]=getPSTHbyUnit(s,stimulus,ops.fs,rasterWindow,binSize,artefactRemoval);
                    s.response.visualOpto=getUnitResponsiveness(raster,conditions{cond});
                    if verbose; plotPSTH(s,conditions{cond},ops.dirOUT); end
                else
                    s.PSTHvisualOpto=[];
                    raster=[];
                    s.response.visualOpto=[];
                end
                
            case 'LaserOnly'
                if ~isempty(stim.laser)
                    stimulus=stim.laser;
                    artefactRemoval=[-51,49;149,153]*10^-3; %s
                    [s.PSTHlaser,s.PSTHbins,raster]=getPSTHbyUnit(s,stimulus+0.05*ops.fs,ops.fs,rasterWindow,binSize,artefactRemoval);
                    s.response.laser=getUnitResponsiveness(raster,conditions{cond});
                    if verbose; plotPSTH(s,conditions{cond},ops.dirOUT); end
                else
                    s.PSTHlaser=[];
                    raster=[];
                    s.response.laser=[];
                end
        end

        if verbose
            plotRaster(raster,s,conditions{cond},ops.dirOUT)
        end
    end      
    
end


function [PSTH,PSTHbins,raster]=getPSTHbyUnit(s,stimulus,sr,window,binSize,artefactRemoval)

    for suIdx=1:numel(s.suid)
        spiketimes=s.st(ismember(s.sclu,s.suid(suIdx)));
        [PSTH(suIdx,:),PSTHbins,raster(suIdx,:)]=getPSTH(spiketimes,stimulus,sr,window,binSize,artefactRemoval);
    end
end