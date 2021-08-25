function CSD = getCSD(data,fs,electrodeSpacing,electrodeLength,electrodeType)
% function csdInfo = computeCSD(data,fs,electrodeLength,electrodeSpacing)
% Calculate and average interpolated CSD from data in input.
%
% Inputs: 
%           data -> 3D array [Chan*timestamps*trial];
%           fs -> sampling rate;
%           electrodeLength -> total length of electrode (1 shank) in um;
%           electrodeSpacing -> vertical spacing between contacts in um;
%           electrodeType -> if 'Poly2' some preprocessing is applied.

% Output:
%           CSD -> struct containing rawCSD, interpCSD, Sinks and Sources indices.
    
    fprintf(strcat('Obtaining CSD...','\n'))

    %% Preprocess data
    % Average adiacent channels to have equally spacing channels (A1x32 electrode only)
    if strcmpi(electrodeType,'Poly2')
        data=(data(1:2:end,:,:)+data(2:2:end,:,:))./2;
        electrodeSpacing=electrodeSpacing*2;
    end

    % Duplicate first and last average channels
    data=[data(1,:,:);data;data(end,:,:)];
    
    % Average across trials
    if size(data,3)>1
        meanERP=mean(data,3);
    else
        meanERP=data;
    end
    
    % Detrend: optimise CSD if adjacent channels have different voltage offset baseline
    meanERP=meanERP';
    meanERP=detrend(meanERP,'constant'); 
    meanERP=meanERP';
    
    %% Compute CSD
%     for i=2:size(meanERP,1)-1
        % spatially filter
%         meanERP(i,:)=0.23*meanERP(i-1,:)+0.23*meanERP(i+1,:)+0.54*meanERP(i,:);% sinks are upwards! https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4617414/
%     end

    for i=2:size(data,1)-1
        csd(i-1,:)=-((meanERP(i-1,:)+meanERP(i+1,:)-2*meanERP(i,:))./(electrodeSpacing)^2).*0.30;
    end

    %% Post-processing CSD
    %Smooth csd with moving mean window 
    csd=movmean(csd,0.01*fs,2); %10ms smoothing
    
    %Source and Sink from each channel
    [~,CSD.Source]=max(csd(:,1*fs:end),[],2);
    [~,CSD.Sink]=min(csd(:,1*fs:end),[],2);
    
    % Interpolate CSD 
    [X,Y]=meshgrid(1:1:size(csd,2),1:electrodeSpacing/10^-6:electrodeLength);
    [X2,Y2]=meshgrid(1:1:size(csd,2),1:Y(end,1));
    interpCSD=interp2(X,Y,csd,X2,Y2,'linear');
    
    %% Build output structure
    CSD.raw=csd;
    CSD.Source=CSD.Source+1*fs;
    CSD.Sink=CSD.Sink+1*fs;
    CSD.interp=interpCSD;
    
end
