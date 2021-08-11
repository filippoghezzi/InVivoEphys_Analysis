function spindleBurst=InVivo_dataProcessing_baseline_spindleBurst_BurstFeatures(envelope,SR,threshold,fx)    
%Function to calculate mean burst amplitude and duration for the envelope
%signal given as input according to several threhsold values.
%Inputs:
%envelope -> signal to analyze;
%SR -> sampling rate of the signal;
%threshold -> threshold for burst determination
%Output:
%MeanBurstAmplitude
%MeanBurstDuration

count=0;
nBurst=0;
BurstAmp=0;
Amplitude=[];
Duration=[];


for k=1:length(envelope) %K for scanning each envelope time series
    if envelope(k)>=threshold %Comparison of single time point with threshold
        count=count+1; %Index duration in samples
        if count==1 %For first time count increases count nBurst
            nBurst=nBurst+1; %Index burst
            onset=k; %Index
        end
        if envelope(k) >= BurstAmp %Udapte BurstAmp as trace of higher point in the burst
            BurstAmp=envelope(k);
        end
    elseif count>0 %If burst just finished
        %% count cycles

        if count/SR>=0.1 
            cycles=findpeaks(-zscore(fx(onset:k)),'MinPeakHeight',1);
            numCycles=numel(cycles);
            if numCycles>=3 %Take only burst with duration longer than 0.1 s
                Amplitude(nBurst)=BurstAmp;
                Duration(nBurst)=count/SR;
                burstStart(nBurst)=onset;
                burstEnd(nBurst)=k-1;
            else 
                nBurst=nBurst-1;
                onset=[];
            end
        else 
            nBurst=nBurst-1;
            onset=[];
        end
        count=0;
        BurstAmp=0;
    end
end

%% Output
spindleBurst.amplitude = Amplitude;
spindleBurst.duration = Duration;
spindleBurst.start = burstStart;
spindleBurst.end = burstEnd;


end
