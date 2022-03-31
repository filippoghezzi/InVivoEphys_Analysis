function spindleBurst = InVivo_dataProcessing_baseline_spindleBurst_singleUnit(spindleBurst,s,baseline,folder)
% function spindleBurst = InVivo_dataProcessing_baseline_spindleBurst_singleUnit(spindleBurst,s,endBaseline,folder)
% Analyse entrainment of single units spiking with spindle bursts. Outputs 
% the probability of spikes occurring within spindle bursts over total 
% spikes during the whole baseline, spike latency and latency normalized 
% related to start of every burst. Finally, measures the entrainment of 
% spike with spindle bursts through a ISI shuffling method.
% Inputs: spindleBurst -> struct, obtained from previous functions. 
%         s -> struct, containing spike data.  
%         endBaseline -> numerical, end of baseline (in sample).
%         folder -> string, folder where to save the figure.
% Outputs: spindleBurst -> struct, updated with new variables.
%
    %% Initialize variables
    spikeProbability=nan(numel(s.suid),1);
    spikeLatency=cell(numel(s.suid),1);
    spikeLatency_norm=cell(numel(s.suid),1);
    unitEntrained=zeros(numel(s.suid),1);
    unitEntrained_p=nan(numel(s.suid),1);
    
    %% Main analysis
    figure('units','normalized','outerposition',[0 0 1 1]);
    for su=1:numel(s.suid)  
        subplot(ceil(sqrt(numel(s.suid))),ceil(sqrt(numel(s.suid))),su)
        singleUnitSpikeTimes=double(s.st(s.sclu==s.suid(su)));
        singleUnitSpikeTimes=singleUnitSpikeTimes(singleUnitSpikeTimes>baseline(1) & singleUnitSpikeTimes<baseline(2));
        if numel(singleUnitSpikeTimes)>0
            [spikeProbability(su),spikeLatency{su},spikeLatency_norm{su}]=getBurstSpikeProb(spindleBurst,singleUnitSpikeTimes);
            [unitEntrained(su),unitEntrained_p(su)]=getBurstUnitEntrained(spindleBurst,singleUnitSpikeTimes,spikeProbability(su));
        end
    end  
    export_fig(fullfile(folder,'SpindleBurstEntrainment'),'-tiff','-transparent')
    close

    %% Output
    spindleBurst.unitEntrained=unitEntrained;
    spindleBurst.unitEntrained_p=unitEntrained_p;
    spindleBurst.spikeProbability=spikeProbability;
    spindleBurst.spikeLatency=spikeLatency;
    spindleBurst.spikeLatency_norm=spikeLatency_norm;
end


function [spikeProbability,spikeLatency,spikeLatency_norm]=getBurstSpikeProb(SB,spiketimes)
% Measure ratio between the number of spikes within spindle bursts and
% total spikes in baseline (spikeProbability). Records also the latency of
% spikes within each spindle bursts as raw latency(spikeLatency) or
% normalized on the total duration of the burst (spikeLatency_norm).

    spikeLatency=[];
    spikeLatency_norm=[];
    nSampleBurst=SB.end-SB.start;
    totSpikesBaseline=nnz(spiketimes);
    
    for burst=1:size(SB.start,2)
        spiketimes_aligned=spiketimes-SB.start(burst);
        spikeTimesThisBurst=spiketimes_aligned(spiketimes_aligned>0 & spiketimes_aligned<nSampleBurst(burst));
        spikeLatency=[spikeLatency;spikeTimesThisBurst];
        spikeLatency_norm=[spikeLatency_norm;spikeTimesThisBurst/nSampleBurst(burst)];
    end
    spikeProbability=numel(spikeLatency)./totSpikesBaseline;
end

function [unitEntrained,p]=getBurstUnitEntrained(SB,spiketimes,true_spikeProbability)
% Perform statistical test for entrainment of single units within bursts
% using an ISI shuffling method (Narayanan, N. S., & Laubach, M. (2009). 
% Methods for studying functional interactions among neuronal populations. 
% In Dynamic Brain Imaging (pp. 135-165). Humana Press). Original ISI are 
% shuffled 1000 times and the probability of spikes within spindle bursts
% is measure to obtain a distribution of shuffled spike probability. If the
% true spike probability is always higher than shuffled, then the single
% unit is entrained with p < 0.001.

    ISI=diff(spiketimes);
    
    for i=1:1000
        shuffledISI=ISI(randperm(length(ISI)));
        shuffledSpikeTimes=cumsum([spiketimes(1);shuffledISI]);
        [shuffled_spikeProbability(i),~,~]=getBurstSpikeProb(SB,shuffledSpikeTimes);
    end
    p=nnz(shuffled_spikeProbability>=true_spikeProbability)/1000;
    if p<0.001; unitEntrained=1; else; unitEntrained=0; end
    
    %% Plotting
    hist(shuffled_spikeProbability,30)
    hold on
    plot([true_spikeProbability,true_spikeProbability],[0 120],'r--') 
end