function ops=analyseBaseline(ops,s)
%Comprehensive baseline analysis function. Wrapper for powerSpectrum,
%coherenceSpikeLFP and phaseLockingSpikeLFP functions.
% Inputs: ops -> structure of recording info and results;
%         s -> structure of spike data;

    %% Read baseline LFP data and select spikes
    endBaseline=ops.nSamplesBlocks(1); %Samples
    if endBaseline > 30*60*ops.fs; endBaseline = 30*60*ops.fs; end %Select only first 30 min of baseline recording if longer
    LFP = get_eLFP(ops.fbinary,1,ops.fs,ops.NchanTOT,[0 endBaseline/ops.fs]);  
    
    %% Do analysis
    spectral=analyseBaseline_powerSpectrum(ops,s,'LFP',LFP,'endBaseline',endBaseline,'Plotting',1);
    coherence=analyseBaseline_coherenceSpikeLFP(ops,s,'LFP',LFP,'endBaseline',endBaseline,'Plotting',1);
    phaseLocking=analyseBaseline_phaseLockingSpikeLFP(ops,s,'LFP',LFP,'endBaseline',endBaseline,'Plotting',1);

    %% Output
    ops.coherence=coherence;
    ops.phaseLocking=phaseLocking;
    ops.spectral=spectral;
end