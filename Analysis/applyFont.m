function ax=applyFont(ax,normalization)
    ax.XLim=[-1000 4000];

    if normalization
        ax.YLabel.String='Z-score';
    else
        ax.YLabel.String='Spike/s'; 
    end
    ax.XAxis.Visible='on';
    ax.FontSize=20;
    ax.Box='off';
    ax.LineWidth = 1.5;
    ax.XLabel.String='Time (ms)';


end