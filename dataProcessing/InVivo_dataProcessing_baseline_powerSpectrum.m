function spectral=InVivo_dataProcessing_baseline_powerSpectrum(ops,s,varargin)
% Provide basic analysis of baseline LFP and MUA data. 
% It calculates power spectral density (PSD), MUA power, average MUA rate
% and continuity for each electrode contact; provides spike rate and
% continuity values for single units. 
% Inputs: ops -> structure of recording info and results;
%         s -> structure of spike data;
%         'LFP' -> matrix(ops.NchanTOT,endBaseline) of LFP data
%         'endBaseline' -> numeric value of baseline end in samples
% Output: ops -> struct, containing function outputs.

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
    
    fprintf(strcat('Running spectral analysis...','\n'))

    %% Initialize variable  
    spikeTimes=s.st(s.st<endBaseline);
    spikeClu=s.sclu(s.st<endBaseline);
    numWindows=floor((endBaseline/ops.fs)*2);
    MUA_continuity=zeros(ops.NchanTOT ,1);
    MUA_spikeRate=zeros(ops.NchanTOT ,1);
    LFP_nChannels=size(LFP,1);
    LFP_frequency=1:0.01:49;
    LFP_PSD=zeros(LFP_nChannels,numel(LFP_frequency));
    MUA = loadLFP(ops.fbinary,1,ops.fs,ops.NchanTOT,[0 60*5],'MUA');
    MUA_frequency=400:10:4000;
    MUApower=zeros(ops.NchanTOT,numel(MUA_frequency));
    SU_spikeRate=zeros(numel(s.suid),1);
    SU_continuity=zeros(numel(s.suid),1);
    
    
    %% Get LFP power spectrum 
    for channel=1:LFP_nChannels
        data=gpuArray(LFP(channel,:)); % Use GPU to speed up
        [pxx,LFP_PSD_f] = pwelch(data,ops.fs*3,[],LFP_frequency,ops.fs);
        pxx=pxx.*LFP_PSD_f; %Reduce 1/f component
        pxx=pxx./mean(pxx);  %Normalize PSD to mean power
        LFP_PSD(channel,:)=gather(pxx);
    end
    if LFP_nChannels==32
        LFP_PSD=mean(LFP_PSD(ops.L4,:));
    end
    
    %% MUA and LFP
    for channel=1:ops.NchanTOT
        %% Get spike rate and continuity each contact - MUA
        IDs=s.cids(s.cch==channel);
        st=spikeTimes(ismember(spikeClu,IDs));
        MUA_spikeRate(channel)=numel(st)/(endBaseline/ops.fs);

        silent=0;
        for batchWindow=1:numWindows
            if nnz(st>(batchWindow-1)*0.5*ops.fs & st<(batchWindow*0.5)*ops.fs)<2
                silent=silent+1;
            end
        end
        MUA_continuity(channel)=1-(silent/numWindows);
   
        %% Get MUA power
        data=gpuArray(MUA(channel,:)); % Use GPU to speed up
        [pxx,~] = pwelch(data,ops.fs*3,[],MUA_frequency,ops.fs);
        pxx=mean(pxx);
        MUApower(channel,:)=gather(pxx);
    end 

    %% Single Units
    for singleUnit = 1:numel(s.suid)   
        %% Get spike rate and continuity each SU
        st=spikeTimes(spikeClu==s.suid(singleUnit));
        SU_spikeRate(singleUnit)=numel(st)/(endBaseline/ops.fs);

        silent=0;
        for batchWindow=1:numWindows
            if nnz(st>(batchWindow-1)*0.5*ops.fs & st<(batchWindow*0.5)*ops.fs)<2
                silent=silent+1;
            end
        end
        SU_continuity(singleUnit)=1-(silent/numWindows);
    end    
    
    %% Plot general spectral analysis
    if plotFlag
        figure('units','normalized','outerposition',[0 0 1 1]);
        sg=sgtitle(strcat('Baseline analysis - ',ops.recID,' - P',num2str(ops.age)));
        sg.FontSize=30;

        ax=subplot(2,3,1:2);
        PSDplot=plot(LFP_PSD_f,LFP_PSD);
        PSDplot.LineWidth=2;
        PSDplot.Color='k';
        ax.XLabel.String='Frequency (Hz)';
        ax.YLabel.String='Normalized PSD';
        ax.FontSize=20;
        ax.FontName='Arial';
        ax.Box='off';
        ax.LineWidth=1;
        ax.Title.String='LFP Power spectral density';

        ax=subplot(2,3,4);
        plot(MUA_spikeRate,1:32,'k','LineWidth',2);
        hold on
        plot([0,max(MUA_spikeRate)],[min(ops.L4),min(ops.L4)],'r--','LineWidth',2)
        plot([0,max(MUA_spikeRate)],[max(ops.L4),max(ops.L4)],'r--','LineWidth',2)
        ax.XLabel.String='MUA rate (Hz)';
        ax.YLabel.String='Channel';
        ax.FontSize=20;
        ax.FontName='Arial';
        ax.Box='off';
        ax.YDir='reverse';
        ax.YLim=[1,32];
        ax.LineWidth=1;
        ax.Title.String='MUA spike rate';

        ax=subplot(2,3,5);
        plot(MUA_continuity,1:32,'k','LineWidth',2);
        hold on
        plot([0,max(MUA_continuity)],[min(ops.L4),min(ops.L4)],'r--','LineWidth',2)
        plot([0,max(MUA_continuity)],[max(ops.L4),max(ops.L4)],'r--','LineWidth',2)
        ax.XLabel.String='MUA continuity';
        ax.YLabel.String='Channel';
        ax.FontSize=20;
        ax.FontName='Arial';
        ax.Box='off';
        ax.YDir='reverse';
        ax.YLim=[1,32];
        ax.LineWidth=1;
        ax.Title.String='MUA continuity';

        ax=subplot(2,3,[3,6]);
        plot(MUApower./max(MUApower),1:32,'k','LineWidth',2);
        hold on
        plot([0,1],[min(ops.L4),min(ops.L4)],'r--','LineWidth',2)
        plot([0,1],[max(ops.L4),max(ops.L4)],'r--','LineWidth',2)
        ax.XLabel.String='Normalised MUA (400-4000 Hz) power';
        ax.YLabel.String='Channel';
        ax.FontSize=20;
        ax.FontName='Arial';
        ax.Box='off';
        ax.YDir='reverse';
        ax.YLim=[1,32];
        ax.LineWidth=1;
        ax.Title.String='Average MUA power';

        figname='Baseline_powerSpectrumAndContinuity';
        export_fig(fullfile(ops.dirOUT,figname),'-tiff','-transparent')
        close
    end
    
    %% Output
    spectral.MUA_spikeRate = MUA_spikeRate;
    spectral.MUA_continuity = MUA_continuity;
    spectral.SU_spikeRate = SU_spikeRate;
    spectral.SU_continuity = SU_continuity;
    spectral.LFP_PSD = LFP_PSD;
    spectral.LFP_PSD_f=LFP_PSD_f;
    spectral.MUApower = MUApower;
end