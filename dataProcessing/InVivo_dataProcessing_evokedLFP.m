function [results,ops]=InVivo_dataProcessing_evokedLFP(ops,s,stim,condition)
%
% Inputs: ops -> structure of recording info and results;
%         s -> structure of spike data;    
    
    results=struct;
    ops.condition=condition;

    eLFP = loadLFP(ops.fbinary,stim,ops.fs,ops.NchanTOT,ops.LFPwindow,'LFP');
    
    CSD = getCSD(eLFP(ops.electrodeChannelsForCSD,:,:),ops.fs,ops.electrodeSpacing,ops.electrodeLength,ops.electrodeType);
    
    results.MUA = getMUA(s,stim,ops.fs,ops.LFPwindow,ops.PSTHbinSize,ops.NchanTOT);
    
    ops=InVivo_dataProcessing_evokedLFP_findL4(eLFP,CSD,results.MUA,ops);
    
    [results.LFP, ops]=InVivo_dataProcessing_evokedLFP_getFeaturesLFP(eLFP,ops);
    
    results.rw=InVivo_dataProcessing_evokedLFP_retinalWave(ops,s,eLFP,stim);
    
    results.spectrogram=InVivo_dataProcessing_evokedLFP_spectrogram(squeeze(eLFP(ops.L4best,:,:)),ops, ops.LFPwindow);
end