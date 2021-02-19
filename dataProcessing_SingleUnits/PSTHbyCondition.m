function s = PSTHbyCondition(s,stim,ops,verbose) 

    rasterWindow=[1,5]; %s
    binSize=0.003; %s

    conditions={'Visual','Optotagging','VisualOpto','LaserOnly'};

    for cond=1:numel(conditions)
        switch conditions{cond}
            case 'Visual'
                stimulus=stim.ledR;
                artefactRemoval=[];
                [s.PSTHvisual,s.PSTHbins,raster]=su_PSTH(s,stimulus,ops.fs,rasterWindow,binSize,artefactRemoval);
                s.response.visual=getUnitResponsiveness(raster,conditions{cond});
                if verbose; plotPSTH(s,conditions{cond},ops.dirOUT); end

            case 'Optotagging'
                stimulus=stim.optotagging;
                artefactRemoval=[-1,1;49,53]*10^-3; %s
                [s.PSTHoptotagging,s.PSTHbins,raster]=su_PSTH(s,stimulus,ops.fs,rasterWindow,binSize,artefactRemoval);
                s.response.optotagging=getUnitResponsiveness(raster,conditions{cond});
                if verbose; plotPSTH(s,conditions{cond},ops.dirOUT); end

            case 'VisualOpto'
                stimulus=stim.laserAndLedR;
                artefactRemoval=[-51,49;149,153]*10^-3; %s
                [s.PSTHvisualOpto,s.PSTHbins,raster]=su_PSTH(s,stimulus,ops.fs,rasterWindow,binSize,artefactRemoval);
                s.response.visualOpto=getUnitResponsiveness(raster,conditions{cond});
                if verbose; plotPSTH(s,conditions{cond},ops.dirOUT); end

            case 'LaserOnly'
                if ~isempty(stim.laser)
                    stimulus=stim.laser;
                    artefactRemoval=[-51,49;149,153]*10^-3; %s
                    [s.PSTHlaser,s.PSTHbins,raster]=su_PSTH(s,stimulus+0.05*ops.fs,ops.fs,rasterWindow,binSize,artefactRemoval);
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


function [PSTH,PSTHbins,raster]=su_PSTH(s,stimulus,sr,window,binSize,artefactRemoval)

    for suIdx=1:numel(s.suid)
        spiketimes=s.st(ismember(s.sclu,s.suid(suIdx)));
        [PSTH(suIdx,:),PSTHbins,raster(suIdx,:)]=makePSTHandRaster(spiketimes,stimulus,sr,window,binSize,artefactRemoval);
    end
end