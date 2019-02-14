function [avgCsd,CsdInfo] = CalculateAverageCSD(data,eventTimes,pre,electrodeLength,interval,avg,avgEve)
% This function gets the interpolated CSD of average Events

st = eventTimes(:,1);
et = eventTimes(:,2);

et = et(st>pre) ;
st = st(st>pre);

dur = et-st+1;
mDur = max(dur);

data = (data-mean(data,2))./std(data,[],2);

if (avg == true)
    df(1:size(data,1)/2,:) = (data(1:2:end,:)+data(2:2:end,:))./2;
else
    df = data;
end
data = df;
% dn = (data-mean(data,2))./std(data,[],2);
dn = data;
dn = [dn(1,:); dn; dn(end,:) ];
data = dn;
% for i =1 :size(data,1)
%     data(i,:) = smooth(data(i,:),100);
% end

% data = movmean(data,100,2);

if (avgEve == true)
    % find average event over channels
    avgEvent = zeros(size(data,1),mDur+pre+1);
    for i = 1:length(st)
        d = data(:,st(i)-pre:et(i));
        avgEvent(:,1:size(d,2)) = avgEvent(:,1:size(d,2)) + d;    
    end
    avgEvent = avgEvent./length(st);

    % CSD%%%%%%%%%

    for i = 2:size(dn,1)-1
            % spatially filter
            avgEvent(i,:) = 0.23*avgEvent(i-1,:)+0.23*avgEvent(i+1,:)+0.54*avgEvent(i,:);% sinks are upwards! https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4617414/

    end

    j = 1;
    for i = 2:size(dn,1)-1

            csd(j,:) = (avgEvent(i-1,:)+avgEvent(i+1,:)-2*avgEvent(i,:))./interval^2;
    %       csd(j,:) = 0.54.*df(i,:)+.23.*df(i-1,:)+0.23.*df(i+1,:); % sinks are upwards! https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4617414/
        j = j+1;
    end

    %%%%%
    else

        for i = 2:size(dn,1)-1
            % spatially filter
            dn(i,:) = 0.23*dn(i-1,:)+0.23*dn(i+1,:)+0.54*dn(i,:);% sinks are upwards! https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4617414/
        end

       j = 1;
    for i = 2:size(dn,1)-1

            csd(j,:) = (dn(i-1,:)+dn(i+1,:)-2*dn(i,:))./interval^2;
    %       csd(j,:) = 0.54.*df(i,:)+.23.*df(i-1,:)+0.23.*df(i+1,:); % sinks are upwards! https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4617414/
        j = j+1;
    end

    avgEvent = zeros(size(csd,1),mDur+pre);
    for i = 1:length(st)
        d = csd(:,st(i)-pre:et(i));
        avgEvent(:,1:size(d,2)) = avgEvent(:,1:size(d,2)) + d;    
    end
    csd = avgEvent./length(st);

end



if (avg == true)
    interval = interval.*2;
end

csd = movmean(csd,10,2);
[CsdInfo] = GetLayerData(csd,pre);

[X,Y] = meshgrid(1:1:size(csd,2), 1:interval:electrodeLength);
Y = Y(1:end,:);
X = X(1:end,:);
[X2,Y2] = meshgrid(1:1:size(csd,2), 1:Y(end,1));
Y2 = Y2(1:end,:);
X2 = X2(1:end,:);
avgCsd = interp2(X, Y, csd, X2, Y2, 'linear');

end

function [csdInfo] = GetLayerData(csd,pre)

m = [];
i = [];
M = [];
I = [];

csdRaw = csd;

for j = 1:size(csd,1)
%     [mt,it] = findpeaks(-csd(j,:),'MinPeakDistance',25,'MinPeakProminence',0.5*std(csd(j,:)));
    [mt,it] = max(-csd(j,pre+1:end));
    if (isempty(mt))
        mt = NaN;
        it = NaN;
    end
    m = [m;mt(1)];
    i = [i;it(1)];

%     [mt,it] = findpeaks(csd(j,:),'MinPeakDistance',20,'MinPeakProminence',0.5*std(csd(j,:)));
        [mt,it] = max(csd(j,pre+1:end));
       if (isempty(mt))
        mt = NaN;
        it = NaN;
    end
    M = [M;mt(1)];
    I = [I;it(1)];

end
csdInfo.csd = csd;
csdInfo.Source = [m i ];
csdInfo.Sink = [M I ];

end