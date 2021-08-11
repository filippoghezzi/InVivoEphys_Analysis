function plotControlvsSalB(c,s)

    MUAcontrol = c.MUA.raw;
    MUAsalB = s.MUA.raw;
    
    %% Fast response
    figure('units','normalized','outerposition',[0 0 1 1]);
    ax_control=subplot(1,3,1);
    imagesc(c.MUA.bins,(1:size(MUAcontrol,1)),MUAcontrol,[0,100]);
    ax_control.XLim=[-50 300];
    colormap(ax_control,'hot')
    ax_control.XLabel.String='Time (ms)';
    ax_control.YLabel.String='Channel';
    c_mua = colorbar('Location','eastoutside');
    c_mua.Label.String = 'Spike/s';
    ax_control.Title.String = 'Control';

    ax_salB=subplot(1,3,2);
    imagesc(c.MUA.bins,(1:size(MUAcontrol,1)),MUAsalB,[0,100]);
    ax_salB.XLim=[-50 300];
    colormap(ax_salB,'hot')
    ax_salB.XLabel.String='Time (ms)';
    ax_salB.YLabel.String='Channel';
    c_mua = colorbar('Location','eastoutside');
    c_mua.Label.String = 'Spike/s';
    ax_salB.Title.String = 'SalB';

    %save figure
    figname='MUAdepth_ControlvsSalB_FirstResponse';
    sg=sgtitle(strcat('P',int2str(c.age)));
    sg.FontSize=30;
    export_fig(fullfile(c.dirOUT,figname),'-tiff','-transparent')
    %% Retinal wave
    ax_control.XLim=[500 4500];
    ax_salB.XLim=[500 4500];
    figname='MUAdepth_ControlvsSalB_SecondResponse';
    export_fig(fullfile(c.dirOUT,figname),'-tiff','-transparent')
    close
    
end