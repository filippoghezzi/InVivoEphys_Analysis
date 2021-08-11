function s=InVivo_dataProcessing_singleUnits_getTemplateFeatures(s,ops,varargin)
    
    wf=s.suWf;
    suid=s.suid;
    fs=ops.fs;
    fshigh=ops.fshigh;
    
    if size(varargin,2)==1
        verbose = varargin{1,1};
    elseif size(varargin,2)==2
        verbose=varargin{1,1};
        folder=varargin{1,2};
    else
        verbose=0;
    end
    if verbose
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    
    %Add an high-pass filter to optimise trough-to-peak detection in
    %presence of slow frequency.
    [b,a]=butter(3,fshigh/fs/2,'high');
    filt_wfs = filtfilt(b,a,wf')';
    for suIdx=1:numel(suid)
        filt_wf=filt_wfs(suIdx,:);
        
        filt_wf_ii=diff(diff(filt_wf));
        [~,onsetWf]=findpeaks(-filt_wf_ii,'MinPeakProminence',20);
        if isempty(onsetWf)
            [~,onsetWf]=min(diff(diff(filt_wf)));
        else
            onsetWf=onsetWf(1);
        end
%         baseline=mean(filt_wf(onsetWf-0.0005*fs:onsetWf));
        baseline=mean(filt_wf(1:0.0015*fs));
        
        [troughValue,troughIndex,width,prom]=findpeaks(-filt_wf);
        [~,maxProminence]=max(prom);
        troughValue=-troughValue(maxProminence);
        troughIndex=troughIndex(maxProminence);
        width=width(maxProminence);
        
        [peakValue,peakTroughSample]=max(filt_wf(troughIndex+1:end));
        
        x=[peakTroughSample+troughIndex,peakTroughSample+troughIndex+0.0005*fs];
        
        if x(2)>numel(filt_wf)
            x(2)=numel(filt_wf);
        end
        lineCoefficients=polyfit(x,filt_wf(x),1);
        slope=lineCoefficients (1);
        
        halfWidth(suIdx,1)=width/fs*1000;
        troughPeakTime(suIdx,1)=peakTroughSample/fs*1000;
        peakTroughRatio(suIdx,1)=(peakValue-baseline)/abs(troughValue-baseline);
        endSlope(suIdx,1)=slope;

        if verbose
            subplot(ceil(sqrt(numel(suid))),ceil(sqrt(numel(suid))),suIdx)
            hold on
            plot(wf(suIdx,:),'r--','LineWidth',2)
            plot(filt_wf,'k','LineWidth',1)
            scatter(troughIndex,filt_wf(troughIndex),'r*')
            scatter(x(1),filt_wf(x(1)),'r*')
            scatter(x(2),filt_wf(x(2)),'r*')
            scatter(onsetWf,filt_wf(onsetWf),'r*')
            title(strcat('Unit-',int2str(suid(suIdx))))
            hold off 
        end
    end 
    
    if verbose && size(varargin,2)==2
        figname=fullfile(folder,'Waveforms');
        export_fig(figname,'-tiff','-transparent')
        close 
    end
    
    %% Set output
    s.filt_suWf=filt_wfs;
    s.halfWidth = halfWidth;
    s.troughPeakTime = troughPeakTime;
    s.peakTroughRatio = peakTroughRatio;
    s.endSlope = endSlope;
end