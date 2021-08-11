function phaseLocking=InVivo_dataProcessing_baseline_phaseLockingSpikeLFP(ops,s,varargin)
% Provide spike LFP phase coupling analysis for single units. 
% Wrapper for getPPC function. See there for details.

% Inputs: ops -> structure of recording info and results;
%         s -> structure of spike data;
%         'LFP' -> matrix(ops.NchanTOT,endBaseline) of LFP data
%         'endBaseline' -> numeric value of baseline end in samples
%         'Plotting' -> 1 to plot PPC, vectorLength and vectorAngle for each single unit
% Output: phaseLocking -> struct, containing function outputs.    

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
    
    fprintf(strcat('Obtaining spike-LFP phase locking (PPC)...','\n'))
    %% Initialize variable  
    spikeTimes=s.st(s.st<size(LFP,2));
    spikeClu=s.sclu(s.st<size(LFP,2));
    frequencyBins=[1, 5; 10, 30; 50, 80];
    PPC=nan(numel(s.suid),size(frequencyBins,1));
    vectorLength=nan(numel(s.suid),size(frequencyBins,1));
    vectorAngle=nan(numel(s.suid),size(frequencyBins,1));
    pValue=nan(numel(s.suid),size(frequencyBins,1));
    
    %% Get pairwise phase consistency
    for frequency=1:size(frequencyBins,1)
        filterFrequency=frequencyBins(frequency,:);
        [b,a]=butter(2,filterFrequency/(ops.fs/2),'bandpass');
        
        for singleUnit=1:numel(s.suid)
            singleUnitChannel=s.cch(s.cids==s.suid(singleUnit));
            singleUnitSpikeTimes=spikeTimes(spikeClu==s.suid(singleUnit));
            selectedChannel=singleUnitChannel+1;
            if selectedChannel>32; selectedChannel=32; end
                
            if numel(singleUnitSpikeTimes)>=5
                filteredLFP=filtfilt(b,a,LFP(selectedChannel,:));
                [PPC(singleUnit,frequency), vectorAngle(singleUnit,frequency), vectorLength(singleUnit,frequency), pValue(singleUnit,frequency)]...
                    =getPPC(singleUnitSpikeTimes,filteredLFP);
            end
        end
    end
    
    %% Plotting
    if plotFlag       
        figure('units','normalized','outerposition',[0 0 1 1]);
        sg=sgtitle(strcat('PhaseLocking - VectorAngle -',ops.binaryRoot(15:end),' - P',num2str(ops.age)));
        sg.FontSize=30; 
        frequencyToPlot=frequencyBins(:,2);
        for i=1:numel(s.suid)
            entrained=pValue(i,:)<=0.05;
            ax=subplot(ceil(sqrt(numel(s.suid))),ceil(sqrt(numel(s.suid))),i); 
            hold on
            plot(frequencyToPlot,rad2deg(vectorAngle(i,:)),'k-o'); 
            scatter(frequencyToPlot(entrained),rad2deg(vectorAngle(i,entrained)),'ko','filled');
            ax.YLim=[-180, 180]; 
            ax.XLabel.String='Frequency (Hz)';
            ax.YLabel.String='Phase (degrees)';
            ax.FontName='Arial';
            ax.Box='off';
            ax.LineWidth=1;
            ax.Title.String=strcat('Unit-',int2str(s.suid(i)));
            ax.YTick=(-180:45:180);
            ax.YGrid='on';
        end
        legend('P>0.05','P<=0.05')
        figname='PhaseLocking_VectorAngle_Reduced';
        export_fig(fullfile(ops.dirOUT,figname),'-tiff','-transparent')
        close
    
        figure('units','normalized','outerposition',[0 0 1 1]);
        sg=sgtitle(strcat('PhaseLocking - VectorLength -',ops.binaryRoot(15:end),' - P',num2str(ops.age)));
        sg.FontSize=30; 
        for i=1:numel(s.suid)
            entrained=pValue(i,:)<=0.05;
            ax=subplot(ceil(sqrt(numel(s.suid))),ceil(sqrt(numel(s.suid))),i); 
            hold on;
            plot(frequencyToPlot,vectorLength(i,:),'k-o');
            scatter(frequencyToPlot(entrained),vectorLength(i,entrained),'ko','filled');
            ax.YLim=[0, 0.6]; 
            ax.XLabel.String='Frequency (Hz)';
            ax.YLabel.String='Vector length';
            ax.FontName='Arial';
            ax.Box='off';
            ax.LineWidth=1;
            ax.Title.String=strcat('Unit-',int2str(s.suid(i)));
        end
        legend('P>0.05','P<=0.05')
        figname='PhaseLocking_VectorLength_Reduced';
        export_fig(fullfile(ops.dirOUT,figname),'-tiff','-transparent')
        close
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        sg=sgtitle(strcat('PhaseLocking - Pairwise Phase-Consistency (PPC) -',ops.binaryRoot(15:end),' - P',num2str(ops.age)));
        sg.FontSize=30; 
        for i=1:numel(s.suid)
            entrained=pValue(i,:)<=0.05;
            ax=subplot(ceil(sqrt(numel(s.suid))),ceil(sqrt(numel(s.suid))),i); 
            hold on
            plot(frequencyToPlot,PPC(i,:),'k-o');
            scatter(frequencyToPlot(entrained),PPC(i,entrained),'ko','filled');
            ax.YLim=[-0.1, .6]; 
            ax.XLabel.String='Frequency (Hz)';
            ax.YLabel.String='PPC';
            ax.FontName='Arial';
            ax.Box='off';
            ax.LineWidth=1;
            ax.Title.String=strcat('Unit-',int2str(s.suid(i)));
        end
        legend('P>0.05','P<=0.05')
        figname='PhaseLocking_PPC_Reduced';
        export_fig(fullfile(ops.dirOUT,figname),'-tiff','-transparent')
        close        
    end    
    
    %% Output
    phaseLocking.vectorLength = vectorLength;
    phaseLocking.vectorAngle = vectorAngle;
    phaseLocking.pValue = pValue;
    phaseLocking.PPC = PPC;
end