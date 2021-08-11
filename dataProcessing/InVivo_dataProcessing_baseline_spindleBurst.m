function spindleBurst=InVivo_dataProcessing_baseline_spindleBurst(ops,s,varargin)

    %% Input parser
    p=inputParser;
    addRequired(p,'ops',@(x) isstruct(x));
    addRequired(p,'s', @(x) isstruct(x));

    addOptional(p,'LFP', [], @(x) isnumeric(x));
    addOptional(p,'endBaseline', [], @(x) isnumeric(x));
    addOptional(p,'Plotting', 0, @(x) isnumeric(x));

    parse(p,ops,s,varargin{:});
    
    ops=p.Results.ops;
    plotFlag=p.Results.Plotting;
    if ~isempty(p.Results.endBaseline)
        endBaseline=p.Results.endBaseline; 
    else
        endBaseline=ops.nSamplesBlocks(1); %Samples
        if endBaseline > 30*60*ops.fs; endBaseline = 30*60*ops.fs; end %Select only first 30 min of baseline recording if longer
    end
    
    if ~isempty(p.Results.LFP)
        LFP=p.Results.LFP;
    else
        LFP = loadLFP(ops.fbinary,1,ops.fs,ops.NchanTOT,[0 endBaseline/ops.fs]);  
    end
    
    
    fprintf(strcat('Analysing spindle bursts...','\n'))
    
    %% Get envelope of signal
    smoothSpan=ops.fs*0.5; %Window span (0.2 s) of the smoothing function
    percentile=(55:5:90); %nth percentile
    [b,a] = butter(2,[10,35]/(ops.fs/2),'bandpass');
    
    filt_LFP = filtfilt(b,a,LFP);
    envelope=abs(hilbert(filt_LFP)); %Calculate the Hilber transform of the filtered signal
    envelope=smooth(envelope,smoothSpan)'; %Smooth the Hilbert envelope with smoothspan
    perc=prctile(envelope,percentile); %Calculate nth percentile of smoothened Hilbert transform
    
    %% Find spindle bursts
    threshold = mean(envelope)+std(envelope)*2;
    spindleBurst=InVivo_dataProcessing_baseline_spindleBurst_BurstFeatures(envelope,ops.fs,threshold,filt_LFP);    
    
    %% Plotting
    if plotFlag
        t = (0:1:numel(LFP)-1)/ops.fs;
    
        figure
        % plot(filt_x)
        hold on
        plot(t,envelope,'k')
        plot([0 2000],[std(envelope),std(envelope)]*2,'g')
        for i=1:numel(spindleBurst.start)
            hold on
            plot(t(spindleBurst.start(i):spindleBurst.end(i)),envelope(spindleBurst.start(i):spindleBurst.end(i)),'r')
        end
    end
    
end
    
    
    