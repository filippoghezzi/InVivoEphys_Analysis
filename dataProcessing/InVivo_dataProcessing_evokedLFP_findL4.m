function ops=InVivo_dataProcessing_evokedLFP_findL4(eLFP,CSD,MUA,ops)

    figure('units','normalized','outerposition',[0 0 1 1]);
    
    meanLFP=mean(eLFP,3);
    time=(-ops.fs*ops.LFPwindow(1):ops.fs*ops.LFPwindow(2));
    time=time/ops.fs*1000;
    if ops.age<10
        spacingLFP=200;
        csdXlim=650;
        muaXlim=4950;
    elseif ops.age>10 && ops.age<12
        spacingLFP=500;
        csdXlim=300;
        muaXlim=300;
    elseif ops.age >=12 && ops.age<16
        spacingLFP=700;
        csdXlim=300;
        muaXlim=300;
    else
        spacingLFP=2000;
        csdXlim=300;
        muaXlim=300;    
    end
            
    %LFP plotting and calculate cortical response
    baselineWindow=time<0 & time>=-100;
    conditionWindow=time>5 & time<=1000;
    
    ax_lfp=subplot(1,3,1);
    hold on
    patch([0 100 100 0],[1 1 -7 -7]*10^4,[1 1 0],'EdgeColor','none') 
    for channel=1:size(meanLFP,1) 
        LFP_to_plot=detrend(meanLFP(channel,:))-spacingLFP*(channel-1);
        plot(time(1:end-1),LFP_to_plot,'k','LineWidth',2)
        
        meanBaseline=mean(meanLFP(channel,baselineWindow),2);
        stdBaseline=std(meanLFP(channel,baselineWindow),[],2);
        cutoff=[meanBaseline-3*stdBaseline,meanBaseline+3*stdBaseline];
        responseOnsetIdx = find(meanLFP(channel,conditionWindow)<cutoff(1),1,'first')+find(time==5);
        if ~isempty(responseOnsetIdx)
            V_response(channel,1)=LFP_to_plot(responseOnsetIdx);
            t_response(channel,1)=time(responseOnsetIdx);
        else
            V_response(channel,1)=NaN;
            t_response(channel,1)=NaN;
        end
        [~,V_min_Idx]=min(meanLFP(channel,conditionWindow));
        V_min_Idx=V_min_Idx+find(time==5);
        V_min(channel,1)=LFP_to_plot(V_min_Idx);
        t_min(channel,1)=time(V_min_Idx);
    end
%     plot(t_response,V_response,'r','LineWidth',2)
%     plot(t_min,V_min,'r','LineWidth',2)
    
    minLFP=min(detrend(meanLFP(end,:)))-spacingLFP*(channel-1)-100;
    maxLFP=max(detrend(meanLFP(1,:)))+100;
    
    ax_lfp.YLim=[minLFP maxLFP];
    ax_lfp.XLim=[-50 3000];
    ax_lfp.XLabel.String='Time (ms)';
    ax_lfp.YLabel.String='Voltage (\muV)';
    
    
    %CSD
    ax_csd=subplot(1,3,2);
    imagesc(time,(1:size(CSD.raw,1)),CSD.raw,[min(CSD.raw(:))/5,max(CSD.raw(:))/5]);
    hold on
    plot(CSD.Source/ops.fs*1000-1000,(1:size(CSD.raw,1)),'b')
    plot(CSD.Sink/ops.fs*1000-1000,(1:size(CSD.raw,1)),'r')
    hold off
    ax_csd.XLim=[-50 csdXlim];
    colormap(flipud(jet));    
    ax_csd.XLabel.String='Time (ms)';
    ax_csd.YLabel.String='Channel';
    colorbar('Ticks',[min(CSD.raw(:))/5 max(CSD.raw(:))/5],'TickLabels',{'Sink','Source'},'Location','eastoutside');

    %MUA
    ax_mua=subplot(1,3,3);
    maxMUA=max(MUA.raw(:))/2;
    if maxMUA==0
        maxMUA=1;    
    end
    imagesc(MUA.bins,(1:size(MUA.raw,1)),MUA.raw,[min(MUA.raw(:))/2,maxMUA]);
    xlim([-50 300])
    colormap(ax_mua,'hot')
    ax_mua.XLabel.String='Time (ms)';
    ax_mua.YLabel.String='Channel';
    c_mua = colorbar('Location','eastoutside');
    c_mua.Label.String = 'Spike/s';
    
    sgtitle(sprintf('%s - P%s',ops.recID,int2str(ops.age)))
    
    if ~isfield(ops,'L4')
        %Get L4
        opts.Resize='on';
        opts.WindowStyle='normal';
        opts.Interpreter = 'tex';
        L4channelStr=inputdlg('Enter L4 channel                                           \color{white} .','L4 in CSD',1,{'1-16'},opts);
        L4channelN=str2num(L4channelStr{1});
        ops.L4=L4channelN;

        ops=selectL4channels(ops);

    end
    if strcmp(ops.electrodeType,'Poly2')
        L4channelsToPlot = ops.L4;
    elseif strcmp(ops.electrodeType,'Poly3')&& ~isfield(ops,'L4')
        L4channelsToPlot = L4channelN;
    else 
        L4channelsToPlot =[];
    end   
    
    
    if ops.verbose
        ax_lfp=subplot(1,3,1);
        ax_lfp.XLim=[-50 3000];
        ax_lfp.YLim=[minLFP maxLFP];
        ax_lfp.Box='off';
        ax_lfp.YTick = [];
        ax_lfp.YAxis.Visible='on';
        ax_lfp.YColor=[1 1 1];
        ax_lfp.LineWidth = 1.5;
        ax_lfp.FontSize=20;
        ax_lfp.YLabel.Color=[0 0 0];
        hold on
        plot([0 0],[minLFP maxLFP],'k--','LineWidth',1)
        plot([100 100],[minLFP maxLFP],'k--','LineWidth',1)
        plot([-50 -50],[minLFP minLFP+spacingLFP],'k','LineWidth',1.5)
        text(-200,minLFP+spacingLFP/2,int2str(spacingLFP),'FontSize',20)
        
        % Re-plot interpolated CSD
        ax_csd=subplot(1,3,2);
        imagesc(time,(1:size(CSD.interp,1)),CSD.interp,[min(CSD.interp(:))/5,max(CSD.interp(:))/5]);
        ax_csd.XLim=[-50 250];
        colormap(flipud(jet));    
        ax_csd.XLabel.String='Time (ms)';
        ax_csd.YLabel.String='Depth (\mum)';
        c_csd=colorbar('Ticks',[min(CSD.raw(:))/5 max(CSD.raw(:))/5],'TickLabels',{'Sink','Source'},'Location','eastoutside');
        hold on
        plot([0 0],[0 size(CSD.interp,1)],'k--','LineWidth',1)
        plot([100 100],[0 size(CSD.interp,1)],'k--','LineWidth',1)
        
        if ~isempty(L4channelsToPlot)
            plot([-50 1050],[min(L4channelsToPlot)*size(CSD.interp,1)/ops.NchanTOT min(L4channelsToPlot)*size(CSD.interp,1)/ops.NchanTOT],'k:','LineWidth',1)
            plot([-50 1050],[max(L4channelsToPlot)*size(CSD.interp,1)/ops.NchanTOT max(L4channelsToPlot)*size(CSD.interp,1)/ops.NchanTOT],'k:','LineWidth',1)
            text(-45,max(L4channelsToPlot)*size(CSD.interp,1)/ops.NchanTOT-50,'L4','FontSize',20)
        end


        ax_csd.Box='off';
        ax_csd.LineWidth = 1.5;
        ax_csd.FontSize=20;
        ax_csd.YColor=[0 0 0];


        % Re-plot interpolated MUA
        ax_mua=subplot(1,3,3);
        ax_mua.XLim=[-50 250];
        colormap(ax_mua,'hot')
        ax_mua.XLabel.String='Time (ms)';
        c_mua = colorbar('Location','eastoutside');
        c_mua.Label.String = 'Spike/s';
        hold on
        plot([0 0],[-50 850],'w--','LineWidth',0.5)
        plot([100 100],[-50 850],'w--','LineWidth',1)

        ax_mua.Box='off';
        ax_mua.LineWidth = 1.5;
        ax_mua.FontSize=20;
        ax_mua.YColor=[0 0 0];
        ax_mua.YAxis.Visible='on';
        ax_mua.YTick = [];

        
        %save figure
        figname=strcat('LFP Depth-',ops.condition);
        sg=sgtitle(strcat('P',int2str(ops.age)));
        sg.FontSize=30;
        export_fig(fullfile(ops.dirOUT,figname),'-tiff','-transparent')
    end
    
    close
end

function ops = selectL4channels(ops)

    L4channels = ops.L4;
    if strcmp(ops.electrodeType,'Poly2')
        if size(L4channels,2)==1
            L4channels=[L4channels.*2-1,L4channels.*2];
            L4channels=[min(L4channels)-1 L4channels max(L4channels)+1];
            L4channels=sort(L4channels,'ascend');
        elseif size(L4channels,2)==2 || size(L4channels,2)==3 
            L4channels=[L4channels.*2-1,L4channels.*2];
            L4channels=sort(L4channels,'ascend');
        else
            error('Too many channels selected')
        end
        
    elseif strcmp(ops.electrodeType,'Poly3')
    
        L4channels = ops.idxCentralLine(L4channels);
        if L4channels(1)==1
            L4channels = [1:max(L4channels)+2];   
        else
            L4channels = [min(L4channels)-2:max(L4channels)+2];
        end
    end
    
    ops.L4=L4channels;
    ops.L2=1:min(ops.L4)-1;
    ops.L5=max(ops.L4)+1:ops.NchanTOT;
    
end