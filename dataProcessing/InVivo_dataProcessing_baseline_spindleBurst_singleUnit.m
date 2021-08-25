function spindleBurst = InVivo_dataProcessing_baseline_spindleBurst_singleUnit(spindleBurst,s,endBaseline)

            
    
    
    for singleUnit=1:numel(s.suid)
        tmp_spikeLatency=[];
        tmp_spikeLatency_norm=[];
        singleUnitSpikeTimes=double(s.st(s.sclu==s.suid(singleUnit)));
        totSpikesBaseline(singleUnit)=nnz(singleUnitSpikeTimes<endBaseline);
        for burst=1:size(spindleBurst.start,2)
            singleUnitSpikeTimesThisBurst=singleUnitSpikeTimes-spindleBurst.start(burst);
            nsampleBurst=spindleBurst.end(burst)-spindleBurst.start(burst);
            spikeTimesThisBurst=singleUnitSpikeTimesThisBurst(singleUnitSpikeTimesThisBurst>0 & singleUnitSpikeTimesThisBurst<nsampleBurst);
            tmp_spikeLatency=[tmp_spikeLatency;spikeTimesThisBurst];
            tmp_spikeLatency_norm=[tmp_spikeLatency_norm;spikeTimesThisBurst/(spindleBurst.end(burst)-spindleBurst.start(burst))];
            spikeLatency{burst,singleUnit}=spikeTimesThisBurst;
        end
        spikeLatency_norm{singleUnit}=tmp_spikeLatency_norm;

    end  
    spikeProbability=cellfun(@numel,spikeLatency_norm)./totSpikesBaseline;
    
    
    %% Output
    spindleBurst.spikeProbability=spikeProbability;
    spindleBurst.spikeLatency=spikeLatency;
    spindleBurst.spikeLatency_norm=spikeLatency_norm;

    
end