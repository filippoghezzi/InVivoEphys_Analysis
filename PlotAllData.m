function []= PlotAllData(data,stim,color)

% plots all channel
% figure
time=(1:size(data,2));

if  size(data,2)==1    
    hold on
    plot(time-stim,data(:,1),'k','LineWidth',1); 
    
%     title(strcat('Ch',int2str(1)));
    
else

    for i = 1:size(data,1)
        hold on
%         subplot(4,round(size(data,1)/4),i);
%         title(strcat('Ch',int2str(i)));
        if strcmp(color,'Grey')
            plot(time-stim,data(i,:),'Color',[0.5 0.5 0.5],'LineWidth',0.5);    
        else
            plot(time-stim,data(i,:),'Color','k','LineWidth',1);    
        end
    end

    
end
hold off
end