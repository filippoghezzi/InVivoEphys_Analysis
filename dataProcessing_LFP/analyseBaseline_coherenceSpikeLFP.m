function coherence=analyseBaseline_coherenceSpikeLFP(ops,s,varargin)
% Provide spike LFP coherence analysis for single units. 
% Wrapper for getCoherence function:
% Timothy Olsen (2020). Spike-signal coherence (https://www.mathworks.com/matlabcentral/fileexchange/72755-spike-signal-coherence), MATLAB Central File Exchange. Retrieved October 16, 2020.

% Inputs: ops -> structure of recording info and results;
%         s -> structure of spike data;
%         'LFP' -> matrix(ops.NchanTOT,endBaseline) of LFP data
%         'endBaseline' -> numeric value of baseline end in samples
%         'Plotting' -> 1 to plot coherence curve each single unit
% Output: coherence -> struct, containing function outputs.

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
        LFP = get_eLFP(ops.fbinary,1,ops.fs,ops.NchanTOT,[0 endBaseline/ops.fs]);  
    end
    
    %% Initialize variable  
    spikeTimes=s.st(s.st<size(LFP,2));
    spikeClu=s.sclu(s.st<size(LFP,2));
    coherenceWindowSize=0.1; % seconds
    iterations=100;
    
    %% Get coherence
    for singleUnit=1:numel(s.suid)
        singleUnitChannel=s.cch(s.cids==s.suid(singleUnit));
        singleUnitSpikeTimes=spikeTimes(spikeClu==s.suid(singleUnit));
        LFPChannel=singleUnitChannel+1;
        if LFPChannel>32; LFPChannel=32; end
        if numel(singleUnitSpikeTimes)>=5                
           [~,coherenceCurve(singleUnit,:),~,~,coherenceCurve_shuff_mean(singleUnit,:),coherenceCurve_shuff_SD(singleUnit,:),coherenceCurvefreqs]...
               =getCoherence(double(singleUnitSpikeTimes),LFP(LFPChannel,:),ops.fs,coherenceWindowSize,iterations,0); 
        elseif singleUnit==numel(s.suid)
            coherenceCurve(singleUnit,:)=0;
            coherenceCurve_shuff_mean(singleUnit,:)=0;
            coherenceCurve_shuff_SD(singleUnit,:)=0;
        end
    end

    %% Plotting
    if plotFlag    
        figure('units','normalized','outerposition',[0 0 1 1]);
        sg=sgtitle(strcat('Coherence -',ops.binaryRoot(15:end),' - P',num2str(ops.age)));
        sg.FontSize=30; 
        for singleUnit=1:numel(s.suid)
            ax=subplot(ceil(sqrt(numel(s.suid))),ceil(sqrt(numel(s.suid))),singleUnit); 
            hold on
            plot(coherenceCurvefreqs,coherenceCurve(singleUnit,:),'k-'); 
            plot(coherenceCurvefreqs,coherenceCurve_shuff_mean(singleUnit,:),'r');
            upperSD = coherenceCurve_shuff_mean(singleUnit,:) + coherenceCurve_shuff_SD(singleUnit,:);
            lowerSD = coherenceCurve_shuff_mean(singleUnit,:) - coherenceCurve_shuff_SD(singleUnit,:);
            xAxisSD = [coherenceCurvefreqs fliplr(coherenceCurvefreqs) coherenceCurvefreqs(1)];
            fill(xAxisSD,[lowerSD fliplr(upperSD) lowerSD(1)],'r','faceAlpha',0.1,'edgeAlpha',0); 
            ax.XLim=[0, 50];
            ax.YLim=[0, 1]; 
            ax.XLabel.String='Frequency (Hz)';
            ax.YLabel.String='Coherence';
            ax.FontName='Arial';
            ax.Box='off';
            ax.LineWidth=1;
            ax.Title.String=strcat('Unit-',int2str(s.suid(singleUnit)),'- nspikes=',int2str(numel(spikeTimes(spikeClu==s.suid(singleUnit)))));
        end
        legend('P>0.05','P<=0.05')
        figname='Coherence';
        export_fig(fullfile(ops.dirOUT,figname),'-tiff','-transparent')
        close
    end
        
    %% Output
    coherence=struct;
    coherence.curve=coherenceCurve;
    coherence.curve_shuff_mean=coherenceCurve_shuff_mean;
    coherence.curve_suff_SD=coherenceCurve_shuff_SD;
    coherence.freqs=coherenceCurvefreqs;
end