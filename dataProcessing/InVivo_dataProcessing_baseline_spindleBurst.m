function spindleBurst=InVivo_dataProcessing_baseline_spindleBurst(ops,s,LFP,baseline,varargin)

    %% Input parser
    p=inputParser;
    addRequired(p,'ops',@(x) isstruct(x));
    addRequired(p,'s', @(x) isstruct(x));
    addRequired(p,'LFP', @(x) isnumeric(x));
    addRequired(p,'endBaseline', @(x) isnumeric(x));
    addOptional(p,'Plotting', 0, @(x) isnumeric(x));
    addOptional(p,'Condition','Control');
    parse(p,ops,s,LFP,baseline,varargin{:});
    plotFlag=p.Results.Plotting;
    fprintf(strcat('Analysing spindle bursts...','\n'))
    
    %% Get envelope of signal
    smoothSpan=ops.fs*0.4; %Window span (0.2 s) of the smoothing function
    [b,a] = butter(2,[10,35]/(ops.fs/2),'bandpass');
    useGPU=0;
    if useGPU
        LFP=gpuArray(LFP);
        filt_LFP = filter(b, a, LFP);
        filt_LFP = flipud(filt_LFP);
        filt_LFP = filter(b, a, filt_LFP);
        filt_LFP = flipud(filt_LFP);
        envelope=abs(hilbert(filt_LFP)); %Calculate the Hilber transform of the filtered signal
        envelope=smooth(envelope,smoothSpan)'; %Smooth the Hilbert envelope with smoothspan
        envelope=gather(envelope);
    else
        filt_LFP = filtfilt(b,a,LFP);
        envelope=abs(hilbert(filt_LFP)); %Calculate the Hilber transform of the filtered signal
        envelope=smooth(envelope,smoothSpan)'; %Smooth the Hilbert envelope with smoothspan
    end
    
    
    
    %% Find spindle bursts
    threshold = median(envelope)+std(envelope)*2;
    spindleBurst=InVivo_dataProcessing_baseline_spindleBurst_BurstFeatures(envelope,ops.fs,threshold,filt_LFP);  
    spindleBurst.start=spindleBurst.start+baseline(1);
    spindleBurst.end=spindleBurst.end+baseline(1);
    spindleBurst=InVivo_dataProcessing_baseline_spindleBurst_singleUnit(spindleBurst,s,baseline,ops.dirOUT);
    
%     %% Plotting
%     if plotFlag
%         t = (0:1:numel(LFP)-1)/ops.fs;
%     
%         figure
%         % plot(t,filt_LFP)
%         hold on
%         plot(t,LFP,'k')
% %         plot([0 2000],[threshold,threshold],'g')
%         for i=1:numel(spindleBurst.start)
%             hold on
%             plot(t(spindleBurst.start(i):spindleBurst.end(i)),LFP(spindleBurst.start(i):spindleBurst.end(i)),'r')
%         end
%     end
    
    
end