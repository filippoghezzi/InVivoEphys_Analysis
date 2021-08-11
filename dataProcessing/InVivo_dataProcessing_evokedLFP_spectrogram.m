function spectrogram = InVivo_dataProcessing_evokedLFP_spectrogram(data,ops,window)
    
    fprintf(strcat('Obtaining spectrogram...','\n'))

    fs=ops.fs;

    time=(-fs*window(1):fs*window(2));
    time=time/fs*1000;
    time=time(1:end-1);
    
    %Select only portion of data
%     time=time(time>-50 & time<4500);
%     data=data(time>-50 & time<4500,:);
    
    [CFS, f] = getSpectrogram(data,fs);

%     CFS=CFS./mean(CFS(:,1:0.05*fs),2);%min(CFS,[],'all');
    CFS=zscore(CFS,[],2);
    
    %% Plot spectrogram
    ax_spectrogram=subplot(1,1,1);
    matlab.graphics.interaction.disableDefaultAxesInteractions(ax_spectrogram);
    hs = surface('Parent',ax_spectrogram,...
        'XData',[min(time) max(time)],'YData',[max(f(f<200)) min(f(f<200))],...
        'CData',zscore(CFS(f<200,:),[],2), 'ZData', zeros(2,2), ...
         'CDataMapping','scaled', 'FaceColor','texturemap', 'EdgeColor', 'none');
    ax_spectrogram.CLim=[0,5];
    ax_spectrogram.YLim = [1,50];
    ax_spectrogram.YMinorTick = 'on';
    ax_spectrogram.XLim = [-1000 5000];
    ax_spectrogram.Layer = 'top';
    ax_spectrogram.YDir = 'normal';
    ax_spectrogram.YTick = [1,2,3,5,10,20,30,50,100];
    ax_spectrogram.YScale = 'log';
    ax_spectrogram.XLabel.String='Time (ms)'; 
    ax_spectrogram.YLabel.String='Frequency (Hz)'; 
    ax_spectrogram.FontSize=20;
    colormap(ax_spectrogram,'jet');
    hcol = colorbar('peer', ax_spectrogram);
    hcol.Label.String = 'Fold change';
    
    %% Save figure
    figname=strcat('Spectrogram L4',ops.condition);
    sgtitle(figname)
    export_fig(fullfile(ops.dirOUT,figname),'-tiff','-transparent')
    close
    
    %% Make output
    spectrogram=struct;
    spectrogram.CFS=CFS;
    spectrogram.freq=f;
    
end