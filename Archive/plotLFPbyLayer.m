function plotLFPbyLayer(LFPtrace,MUAtrace,MUA,ops)
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    time=(-ops.fs*ops.LFPwindow(1):ops.fs*ops.LFPwindow(2));
    time=time/ops.fs*1000;
    time=time(1:end-1);
    LFPtrace=mean(LFPtrace,3);
    j=0;
    for i=1:size(LFPtrace,1)
        %% Plot traces
        ax_lfp=subplot(3,2,i+j);
        LFP_plot=detrend(LFPtrace(i,:));
        MUA_plot=(MUAtrace(i,:)+3*min(LFP_plot))./2;
        
        patch([0 100 100 0],[min(MUA_plot)-100 min(MUA_plot)-100 max(LFP_plot)+100 max(LFP_plot)+100],'y','EdgeColor','none') 
        hold on
        plot(time,LFP_plot,'k','LineWidth',2)
        plot(time,MUA_plot,'k')
        

        ax_lfp.YLim=[min(MUA_plot)-100 max(LFP_plot)+100];
        ax_lfp.XLim=[-50 3000];
        ax_lfp.XLabel.String='Time (ms)';
        ax_lfp.YLabel.String='Voltage (\muV)';
        ax_lfp.YAxis.Visible='off';
        ax_lfp.FontSize=20;
        ax_lfp.Box='off';
        ax_lfp.LineWidth = 1.5;
        
        


        
        %% Plot PSTH
        ax_mua=subplot(3,2,i+j+1);
        PSTH=smooth(MUA(i).raw,5);
        patch([0 100 100 0],[min(PSTH)-100 min(PSTH)-100 max(PSTH)+100 max(PSTH)+100],'y','EdgeColor','none') 
        hold on
        plot(MUA(i).bins,PSTH,'k','LineWidth',2)
        if max(PSTH)~=0
            ax_mua.YLim=[0 max(PSTH)];
        end
        ax_mua.XLim=[-50 3000];
        ax_mua.XLabel.String='Time (ms)';
        ax_mua.YLabel.String='Spike/s';
        ax_mua.Box='off';
        ax_mua.FontSize=20;
        ax_mua.LineWidth = 1.5;


    
        
        j=j+1;
    end

    
    %% Save figure
    figname=strcat('MUA PSTH by layer-',ops.condition);
    sgtitle(figname)
    export_fig(fullfile(ops.dirOUT,figname),'-tiff','-transparent')
    close
end