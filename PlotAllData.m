function []= PlotAllData(data)

% plots all channel

if  size(data,2)==1    
    hold on
    plot(data(:,1)); 
    
    title(strcat('Ch',int2str(1)));
    
else

    for i = 1:size(data,1)
        hold on
        subplot(4,round(size(data,1)/4),i);
        title(strcat('Ch',int2str(i)));
        plot(data(i,:));    
    end

    
end
hold off
end