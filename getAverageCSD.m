function [interp_csd,CSDinfo]=getAverageCSD(data,eventTime,preEventTime,electrodeLength,electrodeSpacing,averageChannels,averageEvents,L4,mua,stimulus)
% This function gets the interpolated CSD of average Events
    
    
    %% Determine event start and end Idx
    eventStart=eventTime(:,1);
    eventEnd=eventTime(:,2);
    
    %events must occur after preEventTime ms after recording start
    eventEnd=eventEnd(eventStart>preEventTime); 
    eventStart=eventStart(eventStart>preEventTime);

    maxEventDuration=max(eventEnd-eventStart);

    %Plot eLFP
    meanERP=zeros(size(data,1),maxEventDuration+preEventTime);
        
    for i=1:length(eventStart)
        if ~(eventStart(i)+maxEventDuration-1 > size(data,2))
            ERP=data(:,eventStart(i)-preEventTime:eventStart(i)+maxEventDuration-1);
        end
        meanERP=meanERP+ERP;    
    end
    meanERP=meanERP./length(eventStart);
    
    if L4
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,3,1)
        for channel=1:size(ERP,1)
            plot(ERP(channel,:)-channel*200,'k')
            hold on
        end
        hold off
        xlim([1000 2000])
    end
    
    %% Preprocess data: Z-scoring, average channels 2-by-2, duplicate first and last channels
    data=zscore(data,[],2);

    %Average adiacent channels to have equally spacing channels (A1x32 electrode only)
    if averageChannels
        data=(data(1:2:end,:)+data(2:2:end,:))./2;
    end

    %Duplicate first and last average channels
    data=[data(1,:);data;data(end,:)];
    
    %% Compute CSD from ERP or whole data
    if averageEvents
        %% Average event-related LFP (ERP)
        meanERP=zeros(size(data,1),maxEventDuration+preEventTime);
        
        for i=1:length(eventStart)
            if ~(eventStart(i)+maxEventDuration-1 > size(data,2))
                ERP=data(:,eventStart(i)-preEventTime:eventStart(i)+maxEventDuration-1);
            end
            meanERP=meanERP+ERP;
        end
        meanERP=meanERP./length(eventStart);

        %% Compute CSD
        for i=2:size(data,1)-1
            % spatially filter
            meanERP(i,:)=0.23*meanERP(i-1,:)+0.23*meanERP(i+1,:)+0.54*meanERP(i,:);% sinks are upwards! https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4617414/
        end

        for i=2:size(data,1)-1
            csd(i-1,:)=(meanERP(i-1,:)+meanERP(i+1,:)-2*meanERP(i,:))./electrodeSpacing^2;
        end

    else
        %% Compute CSD
        for i=2:size(data,1)-1
            % spatially filter
            data(i,:)=0.23*data(i-1,:)+0.23*data(i+1,:)+0.54*data(i,:);% sinks are upwards! https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4617414/
        end

        for i=2:size(data,1)-1
            csd(i-1,:)=(data(i-1,:)+data(i+1,:)-2*data(i,:))./electrodeSpacing^2;
        end

        %% Event related CSD
        meanER_CSD=zeros(size(csd,1),maxEventDuration+preEventTime);
        
        for i=1:length(eventStart)
            ER_
            CSD=csd(:,eventStart(i)-preEventTime:eventEnd(i));
            meanER_CSD=meanER_CSD+ER_CSD;    
        end
        csd=meanER_CSD./length(eventStart);

    end

    %% Post-processing CSD
    if averageChannels
        electrodeSpacing=electrodeSpacing.*2;
    end

    %Smooth csd with moving mean window 
    csd=movmean(csd,10,2);
    
    %Source and Sink from each channel
    [CSDinfo] = getLayerData(csd,preEventTime);
    
    % Interpolate CSD on more refined grid. 
    [X,Y]=meshgrid(1:1:size(csd,2),1:electrodeSpacing:electrodeLength);
    [X2,Y2]=meshgrid(1:1:size(csd,2),1:Y(end,1));
    interp_csd=interp2(X,Y,csd,X2,Y2,'linear');
    
    %% Determine Layer 4
    if L4
        
        subplot(1,3,2)
        imagesc(CSDinfo.rawCSD,[min(CSDinfo.rawCSD(:))/5,max(CSDinfo.rawCSD(:))/5])
        colormap('jet')
        hold on
        plot(CSDinfo.Source(:,2)+preEventTime,(1:16),'b')
        plot(CSDinfo.Sink(:,2)+preEventTime,(1:16),'r')
        xlim([1000 2000])
        hold off
        

        subplot(1,3,3)
        k=1;
        for i=1:2:32
            muaTimes=[]; 
            for j=1:size(stimulus,1)
                tmpmua = mua(ismember(mua(:,3),[i,i+1]),:);
                tmpSpikes=tmpmua(:,4)-stimulus(j,1);
                muaTimes=[muaTimes;tmpSpikes(tmpSpikes>0&tmpSpikes<20000)];
            end
            muaTimes=double(muaTimes)./20; % Convert to double and divide by 20 to get time in ms
            [eventRelatedMUA(k,:),~]=histcounts(muaTimes,600);
            k=k+1;
        end
        imagesc(eventRelatedMUA)
        
        
        opts.Resize='on';
        opts.WindowStyle='normal';
        opts.Interpreter = 'tex';
        L4channel=inputdlg('Enter L4 channel                                           \color{white} .','L4 in CSD',1,{'1-16'},opts);
        CSDinfo.L4=[str2num(L4channel{1}).*2-1,str2num(L4channel{1}).*2];
        close
    end
end

function [csdInfo] = getLayerData(csd,pre)
    minV=[];
    MinIdx=[];
    maxV=[];
    maxIdx=[];

    for j = 1:size(csd,1)
        
        %Min value in each channel (Source), from stimulus onset afterwards
        [SourceVal,SourceIdx]=min(csd(j,pre+1:end));
        if isempty(SourceVal)
            SourceVal=NaN;
            SourceIdx=NaN;
        end
        minV=[minV;SourceVal(1)];
        MinIdx=[MinIdx;SourceIdx(1)];
       
        %Max value in each channel(Sink), from stimulus onset afterwards
        [SinkVal,SinkIdx]=max(csd(j,pre+1:end));
        if isempty(SinkVal)
            SinkVal=NaN;
            SinkIdx=NaN;
        end
        maxV=[maxV;SinkVal(1)];
        maxIdx=[maxIdx;SinkIdx(1)];

    end
    
    % Build output structure
    csdInfo.rawCSD=csd;
    csdInfo.Source=[minV,MinIdx];
    csdInfo.Sink=[maxV,maxIdx];

end