function [envelope,perc]=InVivo_dataProcessing_baseline_spindleBurst_getEnvelope(Signal,SR,peakPSD)
    % Function to perform the Hilbert transform on a filtered signal, compute
    % its absolute value to find the Hilbert envelope of the original signal, 
    % smooth the envelope and find out the nth percentile of the signal. 
    %Inputs:
    %filteredSignal -> a given signal properly filtered around a frequency
    %band;
    %SR -> sampling rate of the filtered signal.
    %Outputs:
    %envelope -> vector containing the time-series of the original signal;
    %perc -> scalar value of the nth percentile.

    smoothSpan=SR*0.2; %Window span (0.2 s) of the smoothing function
    percentile=(55:5:90); %nth percentile
    [b,a] = butter(2,[peakPSD-3,peakPSD+3]/(SR/2),'bandpass');
    
    x = filtfilt(b,a,Signal);
%     gpuDevice(1)
%     x=gpuArray(x);
    x=abs(hilbert(x)); %Calculate the Hilber transform of the filtered signal
    envelope=smooth(x,smoothSpan)'; %Smooth the Hilbert envelope with smoothspan
    perc=prctile(envelope,percentile); %Calculate nth percentile of smoothened Hilbert transform
%     envelope=gather(envelope);
%     perc=gather(perc);
end