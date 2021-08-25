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
    smoothSpan=ops.fs*0.4; %Window span (0.2 s) of the smoothing function
    [b,a] = butter(2,[10,35]/(ops.fs/2),'bandpass');
    
    filt_LFP = filtfilt(b,a,LFP);
    envelope=abs(hilbert(filt_LFP)); %Calculate the Hilber transform of the filtered signal
    envelope=smooth(envelope,smoothSpan)'; %Smooth the Hilbert envelope with smoothspan
    
    %% Find spindle bursts
    threshold = median(envelope)+std(envelope)*2;
    spindleBurst=InVivo_dataProcessing_baseline_spindleBurst_BurstFeatures(envelope,ops.fs,threshold,filt_LFP);    
    spindleBurst=InVivo_dataProcessing_baseline_spindleBurst_singleUnit(spindleBurst,s,endBaseline);
    
    %% Plotting
    if plotFlag
        t = (0:1:numel(LFP)-1)/ops.fs;
    
        figure
        % plot(t,filt_LFP)
        hold on
        plot(t,LFP,'k')
%         plot([0 2000],[threshold,threshold],'g')
        for i=1:numel(spindleBurst.start)
            hold on
            plot(t(spindleBurst.start(i):spindleBurst.end(i)),LFP(spindleBurst.start(i):spindleBurst.end(i)),'r')
        end
    end
    
end
    
    
    