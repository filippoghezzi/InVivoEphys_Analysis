function results=InVivo_dataProcessing_baseline(ops,s,varargin)
%Comprehensive baseline analysis function. Wrapper for powerSpectrum,
%coherenceSpikeLFP and phaseLockingSpikeLFP functions.
% Inputs: ops -> structure of recording info and results;
%         s -> structure of spike data;
    
    %% Input parser
    p=inputParser;
    addRequired(p,'ops',@(x) isstruct(x));
    addRequired(p,'s', @(x) isstruct(x));

    addOptional(p,'Condition',[]);
    parse(p,ops,s,varargin{:});
    
    %% Set start and end baseline according to condition
    if isempty(p.Results.Condition) || strcmp(p.Results.Condition,'Control')
        startBaseline=0;
        endBaseline=ops.nSamplesBlocks(1); %Samples
        if endBaseline > 30*60*ops.fs; endBaseline = 30*60*ops.fs; end %Select only first 30 min of baseline recording if longer

    elseif strcmp(p.Results.Condition,'SalB')
        idxBaselineSalB=find(strcmp(ops.Protocol,'Baseline_K'));
        samplesProtocol=cumsum(ops.nSamplesBlocks);
        startBaseline=samplesProtocol(idxBaselineSalB-1)+1;
        endBaseline=ops.nSamplesBlocks(idxBaselineSalB);
        if endBaseline > 20*60*ops.fs; endBaseline = 20*60*ops.fs; end %Select only first 30 min of baseline recording if longer
        endBaseline=endBaseline+startBaseline;
    end


    %% Read baseline LFP data and select spikes
    fprintf(strcat('Analysing baseline...','\n'))
    LFP = loadLFP_baseline(ops.fbinary,ops.fs,ops.NchanTOT,startBaseline,endBaseline,'LFP');  
    
    %% Do analysis
    results=struct;
    results.durationBaseline=(endBaseline-startBaseline)/ops.fs;
    results.singleUnitFiringFrequency=getSpontaneousFiringSingleUnit(s, startBaseline, startBaseline+endBaseline, ops.fs);
    results.spectral = InVivo_dataProcessing_baseline_powerSpectrum(ops,s,'LFP',LFP(ops.L4best,:),'BaselineWindow',[startBaseline, endBaseline],'Plotting',1);
    results.spindleBurst = InVivo_dataProcessing_baseline_spindleBurst(ops,s,'LFP',LFP(ops.L4best,:),'endBaseline',endBaseline,'Plotting',1);
%     results.coherence = InVivo_dataProcessing_baseline_coherenceSpikeLFP(ops,s,'LFP',LFP,'endBaseline',endBaseline,'Plotting',1);
    results.phaseLocking = InVivo_dataProcessing_baseline_phaseLockingSpikeLFP(ops,s,'LFP',LFP,'endBaseline',endBaseline,'Plotting',1);
end

function singleUnitFiringFrequency=getSpontaneousFiringSingleUnit(s, startBaseline, endBaseline, fs)
    spikeTimes=s.st(s.st>startBaseline & s.st<endBaseline);
    spikeClu=s.sclu(s.st>startBaseline & s.st<endBaseline); 
    singleUnitFiringFrequency=nan(numel(s.suid),1);
    for singleUnit=1:numel(s.suid)
        singleUnitSpikeTimes=spikeTimes(spikeClu==s.suid(singleUnit));
        singleUnitFiringFrequency(singleUnit)=numel(singleUnitSpikeTimes);
    end
    baselineDuration=(endBaseline-startBaseline)/fs;
    singleUnitFiringFrequency=singleUnitFiringFrequency./baselineDuration;
end