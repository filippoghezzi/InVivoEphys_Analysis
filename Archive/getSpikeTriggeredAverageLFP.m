function getSpikeTriggeredAverageLFP(spike,LFP,clusterID,ops)

    % Select only clusterID in input
    spikeTimes=spike.st(spike.sclu==clusterID);
    spikeTimes=spikeTimes(spikeTimes<size(LFP,2));
    spikeChannel=spike.cch(spike.cids==clusterID)+1;
    
    [b,a]=butter(2,[5,15]/(ops.fs/2),'bandpass');
    filtering=1;
    window=0.05;
    
    for i=1:numel(spikeTimes)
        if filtering
            stLFP(i,:)=filtfilt(b,a,LFP(spikeChannel,(spikeTimes(i)-(window*ops.fs):spikeTimes(i)+(window*ops.fs)))); 
        else
            stLFP(i,:)=LFP(spikeChannel,(spikeTimes(i)-(window*ops.fs):spikeTimes(i)+(window*ops.fs))); 
        end
    end
    
    time=linspace(-window,window,size(stLFP,2));
%     figure
    plot(time,mean(stLFP))
end