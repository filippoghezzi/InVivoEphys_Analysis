function [zscoredData, responsive] = zscoreBaseline(data)
    
    %% Z-score for baseline
    baseline=1:300;
    meanData=mean(data(:,baseline),2);
    stdData=std(data(:,baseline),[],2);
    zscoredData=(data-meanData)./stdData;
    zscoredData(isnan(zscoredData(:,1)),:)=zscore(data(isnan(zscoredData(:,1)),:),[],2);
    
    %% Responsive units
    Z=zscore(data,[],'all');
    responsive=any(Z>3,2);


end