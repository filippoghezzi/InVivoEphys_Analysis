close all
clear 
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 
load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\NH_LJB_SingleUnitData.mat')
folderFigures='C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\WT Doctoral Programme in Neuroscience\DPhil in Neuroscience\Manuscripts\V1 S1\S1BF_InVivo';

%% Set group logic arrays
% data=data(~isnan(data.endSlope),:);
data=data(data.Tagging~='SST;NrgKO',:);
data=data(data.brainArea=='S1BF',:);
data=data(data.state=='Urethane',:);

%Layer
tmpLayers(data.Layer==1,1)=categorical(cellstr('L2/3'));
tmpLayers(data.Layer==2,1)=categorical(cellstr('L4'));
tmpLayers(data.Layer==3,1)=categorical(cellstr('L5/6'));
data.Layer=tmpLayers;
clear tmpLayers
L23=data.Layer=='L2/3';
L4=data.Layer=='L4';
L56=data.Layer=='L5/6';

%Age
P5P8=data.Age<9;
P9P13=data.Age>=9 & data.Age<14;
P14P18=data.Age>=14;
data.Dev(P5P8)=categorical(cellstr('P5-P8'));
data.Dev(P9P13)=categorical(cellstr('P9-P13'));
data.Dev(P14P18)=categorical(cellstr('P14-P18'));

%Cell Type
SST=(data.Tagging=='SST' & data.responseTag==1);
Nkx=(data.Tagging=='Nkx2-1' & data.responseTag==1);
Untagged=(~SST & ~Nkx);

SSTopto = data.Tagging=='SST';
Nkxopto = data.Tagging=='Nkx2-1';
isLaser = ~(isnan(data.PSTHlaser(:,1)));

%% Chemogenetics V1
isChemo = ~isnan(data.responseVisual_K);

ChemoRealV1 = {'K18','K24','K35','K36','K38','K39','K43','K44','K45','K47','K6'};
ChemoRealS1 = {'K51','K52','K54','K55','K56','K58','K60','K62','K64','K71','K72','K73','K74'};
ChemoReal = [ChemoRealV1(:)',ChemoRealS1(:)'];
ChemoCtrlGFP = {'K82','K83','K84','K85','K86','K87','K88'}; 
ChemoCtrlSaline = {'K75','K76'};
ChemoNoExpression = {'K48','K4','K9','K10','K17','K20','K21','K23','K37','K41','K42','K46','K50','K53','K59'};

KORD=ismember(data.MouseID,ChemoReal);
GFP=ismember(data.MouseID,ChemoCtrlGFP);
Saline=ismember(data.MouseID,ChemoCtrlSaline);

%% Preprocess data
% [Z_PSTHvisual, responsiveVisual]=zscoreBaseline(data.PSTHvisual);
% Z_PSTHvisualOpto=zscoreBaseline(data.PSTHvisualOpto);
% data.Z_PSTHlaser=zscore(data.PSTHlaser,[],'all');
% Z_PSTHoptotagging=zscoreBaseline(data.PSTHoptotagging);
Z_PSTHvisual=zscore(data.PSTHvisual,[],'all');
% responsive=any(Z_PSTHvisual(:,PSTHbins>0 & PSTHbins<200)>=3,2);
responsive=data.responseWhisker==1;


data.peakVisualFast = max(data.PSTHvisual(:,(PSTHbins>0 & PSTHbins<200)),[],2);
data.peakVisualFast_N = data.peakVisualFast-data.baseline_firing;

data.peakVisualFast_K = max(data.PSTHvisual_K(:,(PSTHbins>0 & PSTHbins<200)),[],2);
data.peakVisualFast_K_N = data.peakVisualFast_K-data.rw_baselineFiring;

data.visualOptoBaseline = mean(data.PSTHvisualOpto(:,(PSTHbins>-200 & PSTHbins<0)),2);
data.peakVisualOptoFast = max(data.PSTHvisualOpto(:,(PSTHbins>0 & PSTHbins<200)),[],2);
data.peakVisualOptoFast_N = data.peakVisualOptoFast-data.visualOptoBaseline;

data.optoBaseline = mean(data.PSTHlaser(:,(PSTHbins>-150 & PSTHbins<-50)),2);
data.meanOpto = mean(data.PSTHlaser(:,(PSTHbins>-50 & PSTHbins<50)),2);
data.optoRebound = mean(data.PSTHlaser(:,(PSTHbins>150 & PSTHbins<250)),2);
data.optoChange = log2((data.meanOpto+1)./(data.optoBaseline+1));
data.reboundChange = log2((data.optoRebound+1)./(data.optoBaseline+1));


data.rw_firing_N=data.rw_firing-data.rw_baselineFiring;
data.rw_firing_K_N=data.rw_firing_K-data.rw_baselineFiring_K;
data.rw_firing_N(data.rw_firing_N>10)=NaN;
data.rw_meanPSTH = mean(data.PSTHvisual(:,PSTHbins>400 & PSTHbins < 4000),2);
data.rwOpto_meanPSTH = mean(data.PSTHvisualOpto(:,PSTHbins>400 & PSTHbins < 4000),2);
data.rw_meanPSTH_K = mean(data.PSTHvisual_K(:,PSTHbins>400 & PSTHbins < 4000),2);

SB_entrained=data.sb_entrainedP<0.05;
SB_entrained_K=data.sb_entrainedP_K<0.05;
PPC_entrained=data.pValuePPC(:,2)<0.05;
PPC_entrained_K=data.pValuePPC_K(:,2)<0.05;

LaserInhibited = data.responseLaser(:,1)==2;
LaserExcited = data.responseLaser(:,1)==1;
LaserRebound = data.responseLaser(:,2)==1;

OptoInhibited = data.responseTag(:,1)==2;
OptoExcited = data.responseTag(:,1)==1;
% OptoRebound = data.responseTag(:,2)==1;

data.peakVisualFast_Change = log2((data.peakVisualFast_K)./(data.peakVisualFast));
data.peakVisualSlow_Change = log2((data.rw_firing_K)./(data.rw_firing));
data.baseline_Change = log2((data.baseline_firing_K)./(data.baseline_firing));
data.fano_Change = log2(data.fanoVisual_K./data.fanoVisual);
data.sbSpikeProb_Change = log2((data.sb_spikeProb_K)./(data.sb_spikeProb));
data.PPC_Change = log2((data.PPC_K(:,2)+0.14)./(data.PPC(:,2)+0.14));

data.peakVisualFast_Change(data.peakVisualFast_Change==-Inf | data.peakVisualFast_Change==+Inf)=NaN; 
data.peakVisualSlow_Change(data.peakVisualSlow_Change==-Inf | data.peakVisualSlow_Change==+Inf)=NaN; 
data.sbSpikeProb_Change(data.sbSpikeProb_Change==-Inf | data.sbSpikeProb_Change==+Inf)=NaN; 
data.PPC_Change(data.PPC_Change==-Inf | data.PPC_Change==+Inf)=NaN; 
data.baseline_Change(data.baseline_Change==-Inf | data.baseline_Change==+Inf)=NaN; 
data.vectorAngle = wrapTo360(rad2deg(data.vectorAngle));







%% Single Units identity v1
%P14-P18
figure('units','normalized','outerposition',[0 0 1 1]);
ax=subplot(1,7,1);
[value,edges]=histcounts(data.troughPeakTime(P14P18&Untagged),0:0.1:3);
h1=barh(edges(1:end-1),value,'LineWidth',2,'EdgeColor','none','FaceColor',[217,217,217]/255);
set(get(h1,'Parent'),'xdir','r')
hold on
[value,edges]=histcounts(data.troughPeakTime((P5P8|P9P13)&Untagged),0:0.1:3);
h2=barh(edges(1:end-1),value,'LineWidth',2,'EdgeColor','none','FaceColor',[115,115,115]/255);
set(get(h2,'Parent'),'xdir','r')
% h2.FaceAlpha=0.7;
ylim([0 2.6])
xlim([0,170])

ax.YLabel.String='Trough to Peak Latency (ms)';
ax.XLabel.String='Cell #';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
grid
ax.GridLineStyle=':';
ax.GridColor=[0,0,0]/255;
legend('Pre','Post')



ax=subplot(1,7,(5:7));
hold on
plot(data.halfWidth(P14P18&Untagged),data.troughPeakTime(P14P18&Untagged),'LineStyle','none','Marker','o', 'MarkerSize',7,'MarkerFaceColor',[217,217,217]/255,'MarkerEdgeColor','none')
plot(data.halfWidth(P14P18&Nkx),data.troughPeakTime(P14P18&Nkx),'LineStyle','none','Marker','o','MarkerSize',8, 'MarkerFaceColor',[215,48,39]/255,'MarkerEdgeColor','none')
plot(data.halfWidth(P14P18&SST),data.troughPeakTime(P14P18&SST),'LineStyle','none','Marker','o','MarkerSize',8, 'MarkerFaceColor',[69,117,180]/255,'MarkerEdgeColor','none')
xlim([0 1])
ylim([0 2.6])
ax.XLabel.String='Spike Half-Width (ms)';
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
ax.YAxis.Visible='off';
% ax.Color=[247,247,247]/255;
grid
ax.GridColor=[0,0,0]/255;
ax.GridLineStyle=':';
ax.Title.String='Post eye opening';
legend('Untagged','Nkx2-1','SST')


ax=subplot(1,7,(2:4));
hold on
plot(data.halfWidth((P5P8|P9P13)&Untagged),data.troughPeakTime((P5P8|P9P13)&Untagged),'LineStyle','none','Marker','o', 'MarkerSize',7,'MarkerFaceColor',[217,217,217]/255,'MarkerEdgeColor','none','HandleVisibility','off')
plot(data.halfWidth((P5P8|P9P13)&Nkx),data.troughPeakTime((P5P8|P9P13)&Nkx),'LineStyle','none','Marker','o','MarkerSize',8, 'MarkerFaceColor',[215,48,39]/255,'MarkerEdgeColor','none')
plot(data.halfWidth((P5P8|P9P13)&SST),data.troughPeakTime((P5P8|P9P13)&SST),'LineStyle','none','Marker','o','MarkerSize',8, 'MarkerFaceColor',[69,117,180]/255,'MarkerEdgeColor','none')
xlim([0 1])
ylim([0 2.6])
ax.XLabel.String='Spike Half-Width (ms)';
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
ax.YAxis.Visible='off';
grid
ax.GridColor=[0,0,0]/255;
ax.GridLineStyle=':';
ax.Title.String='Pre eye opening';

% export_fig(fullfile(folderFigures,'WF_scatter_v1'),'-pdf','-transparent','-nocrop')
% close

%% Single Units identity v2
%P14-P18
figure('units','normalized','outerposition',[0 0 1 0.7]);
ax=subplot(1,8,5);
[value,edges]=histcounts(data.troughPeakTime(P14P18&Untagged),0:0.1:3);
h1=barh(edges(1:end-1),value,'LineWidth',1,'EdgeColor','none','FaceColor',[217,217,217]/255);
set(get(h1,'Parent'),'xdir','r')
ylim([0 2.6])
xlim([0,170])
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.XLabel.String='Cell #';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
grid
ax.GridLineStyle=':';
ax.GridColor=[0,0,0]/255;

ax=subplot(1,8,1);
[value,edges]=histcounts(data.troughPeakTime((P5P8|P9P13)&Untagged),0:0.1:3);
h2=barh(edges(1:end-1),value,'LineWidth',1,'EdgeColor','none','FaceColor',[115,115,115]/255);
set(get(h2,'Parent'),'xdir','r')
ylim([0 2.6])
xlim([0,170])
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.XLabel.String='Cell #';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
grid
ax.GridLineStyle=':';
ax.GridColor=[0,0,0]/255;

ax=subplot(1,8,(6:8));
hold on
plot(data.halfWidth(P14P18&Untagged),data.troughPeakTime(P14P18&Untagged),'LineStyle','none','Marker','o', 'MarkerSize',5,'MarkerFaceColor',[217,217,217]/255,'MarkerEdgeColor','none')
plot(data.halfWidth(P14P18&Nkx),data.troughPeakTime(P14P18&Nkx),'LineStyle','none','Marker','o','MarkerSize',6, 'MarkerFaceColor',[215,48,39]/255,'MarkerEdgeColor','none')
plot(data.halfWidth(P14P18&SST),data.troughPeakTime(P14P18&SST),'LineStyle','none','Marker','o','MarkerSize',6, 'MarkerFaceColor',[69,117,180]/255,'MarkerEdgeColor','none')
xlim([0 1])
ylim([0 2.6])
ax.XLabel.String='Spike Half-Width (ms)';
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
ax.YAxis.Visible='off';
% ax.Color=[247,247,247]/255;
grid
ax.GridColor=[0,0,0]/255;
ax.GridLineStyle=':';
ax.Title.String='Post eye opening';
legend('Untagged','Nkx2-1','SST')


ax=subplot(1,8,(2:4));
hold on
plot(data.halfWidth((P5P8|P9P13)&Untagged),data.troughPeakTime((P5P8|P9P13)&Untagged),'LineStyle','none','Marker','o', 'MarkerSize',5,'MarkerFaceColor',[217,217,217]/255,'MarkerEdgeColor','none','HandleVisibility','off')
plot(data.halfWidth((P5P8|P9P13)&Nkx),data.troughPeakTime((P5P8|P9P13)&Nkx),'LineStyle','none','Marker','o','MarkerSize',6, 'MarkerFaceColor',[215,48,39]/255,'MarkerEdgeColor','none')
plot(data.halfWidth((P5P8|P9P13)&SST),data.troughPeakTime((P5P8|P9P13)&SST),'LineStyle','none','Marker','o','MarkerSize',6, 'MarkerFaceColor',[69,117,180]/255,'MarkerEdgeColor','none')
xlim([0 1])
ylim([0 2.6])
ax.XLabel.String='Spike Half-Width (ms)';
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
ax.YAxis.Visible='off';
grid
ax.GridColor=[0,0,0]/255;
ax.GridLineStyle=':';
ax.Title.String='Pre eye opening';

export_fig(fullfile(folderFigures,'4.4','WF_scatter_v2'),'-pdf','-transparent','-nocrop')
close


%% Extract waveforms
age='P8P13';
dimensionalityReductionMethod='tSNE';

clear Waveforms
if strcmp(age,'P14P18')
    Waveform_raw=data.filt_wf(P14P18,:);
elseif strcmp(age,'P8P13')
    Waveform_raw=data.filt_wf(P5P8|P9P13,:);
end
WF_baseline=mean(Waveform_raw(:,1:50),2);
Z_WF=(Waveform_raw-WF_baseline);
[minWFvalue,minWFidx]=min(Z_WF,[],2);
Z_WF=Z_WF./minWFvalue;
for i=1:size(Z_WF,1); Waveforms(i,:)=Z_WF(i,minWFidx(i)-15:minWFidx(i)+60); end

% figure
% plot(data.wf')
% figure
% plot(Waveforms')

% IDs=unique(data.MouseID);
% for i=1:numel(IDs)
%     subset=data(data.MouseID==IDs(i),:);
%     figure
%     plot(subset.wf')
%     title(IDs(i))
% end

%% Single Unit clustering
% %Dimensionality reduction
if strcmp(dimensionalityReductionMethod,'PCA')
    [~,WF_PCAscore,~,~,~]=pca(Waveforms,'Algorithm','svd','Centered',false);
    data_reduced=WF_PCAscore(:,1:10);
    filename='WF_PCA';
elseif strcmp(dimensionalityReductionMethod,'tSNE')
    tsneOptions=struct;
    tsneOptions.MaxIter=5000;
    tsneOptions.OutputFcn=[];
    tsneOptions.TolFun=1e-10;
    WF_tsneScore=tsne(Waveforms,'Algorithm','exact','Distance','cosine','NumDimensions',2,'Perplexity',20,'Options',tsneOptions);
    data_reduced=WF_tsneScore; 
    filename='WF_tSNE';
end
% figure; hold on;
% sc=scatter(data_reduced(:,1),data_reduced(:,2),'filled');
% sc.Marker='o';
% sc.MarkerFaceColor='k';

% K-means clustering
figure('units','normalized','outerposition',[0 0 1 1]);
Klabel = IterateKMeansClustering(data_reduced, 1:5, 'NumIterations', 5000, 'Distance','sqeuclidean','Verbose',true);
export_fig(fullfile(folderFigures,strcat(filename,'_clustering_',age)),'-tiff','-transparent','-nocrop')
close

colours = {'k', 'r', 'g', 'c', 'm'};
figure('units','normalized','outerposition',[0 0 1 1]);
ax=gca;
% subplot(1,2,1)
hold on
for i=1:numel(unique(Klabel))
    sc=scatter(data_reduced(Klabel==i,1),data_reduced(Klabel==i,2),'filled');
    sc.Marker='o';
    sc.MarkerFaceColor=colours{i};
end
% xlim([0 2.6])
% ylim([0 2.6])
if strcmp(dimensionalityReductionMethod,'PCA')
    ax.XLabel.String='PC 1';
    ax.YLabel.String='PC 2';  
elseif strcmp(dimensionalityReductionMethod,'tSNE')
    ax.XLabel.String='tSNE 1';
    ax.YLabel.String='tSNE 2';
end
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
% ax.YAxis.Visible='off';
grid
export_fig(fullfile(folderFigures,strcat(filename,'_embedding_',age)),'-tiff','-transparent','-nocrop')
close


% figure('units','normalized','outerposition',[0 0 1 1]);
figure
subplot(1,2,2)
hold on
ax=gca;
for i=1:numel(unique(Klabel))
    shadedErrorBar(1:size(Waveforms,2),mean(-Waveforms(Klabel==i,:),1),std(Waveforms(Klabel==i,:),1),'lineProps',colours{i})
end
ax.LineWidth = 1.5;
ax.XAxis.Visible='off';
ax.YAxis.Visible='off';
export_fig(fullfile(folderFigures,strcat(filename,'_trace_',age)),'-tiff','-transparent','-nocrop')
close



% Replot scatter with cluster identity
if strcmp(age,'P14P18')
    subset_troughPeakTime=data.troughPeakTime(P14P18);
    subset_halfWidth=data.halfWidth(P14P18);
    subset_Nkx=Nkx(P14P18);
    subset_Untagged=Untagged(P14P18);
    subset_SST=SST(P14P18);
elseif strcmp(age,'P8P13')
    subset_troughPeakTime=data.troughPeakTime(P5P8|P9P13);
    subset_halfWidth=data.halfWidth(P5P8|P9P13);
    subset_Nkx=Nkx(P5P8|P9P13);
    subset_Untagged=Untagged(P5P8|P9P13);
    subset_SST=SST(P5P8|P9P13);
end

%P14-P18
figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle(age, 'FontSize',30)
ax=subplot(1,4,1);
[value,edges]=histcounts(subset_troughPeakTime,0:0.1:3);
h=barh(edges(1:end-1),value,'w','LineWidth',2);
set(get(h,'Parent'),'xdir','r')
ylim([0 2.6])
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.XLabel.String='Cell #';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
grid

ax=subplot(1,4,(2:4));
hold on
for i=1:numel(unique(Klabel))
    scatter(subset_halfWidth(Klabel==i&subset_Untagged),subset_troughPeakTime(Klabel==i&subset_Untagged),100, colours{i}, 'filled')
end
for i=1:numel(unique(Klabel))
    scatter(subset_halfWidth(Klabel==i&subset_Nkx),subset_troughPeakTime(Klabel==i&subset_Nkx),500, colours{i}, '*')
end
for i=1:numel(unique(Klabel))
    scatter(subset_halfWidth(Klabel==i&subset_SST),subset_troughPeakTime(Klabel==i&subset_SST),500, colours{i}, '+')
end
xlim([0 1])
ylim([0 2.6])
ax.XLabel.String='Spike Half-Width (ms)';
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
ax.YAxis.Visible='off';
grid
export_fig(fullfile(folderFigures,strcat(filename,'_scatter_',age)),'-tiff','-transparent','-nocrop')
close

%% Clustering on extracted waveform features 
WF_features=[data.halfWidth,data.troughPeakTime,data.peakTroughRatio,data.endSlope];
if strcmp(age,'P14P18')
    WF_features=WF_features(P14P18,:);
elseif strcmp(age,'P8P13')
    WF_features=WF_features(P5P8|P9P13,:);
end
WF_features=WF_features(~any(isnan(WF_features),2),:);
WF_features=zscore(WF_features);

if strcmp(dimensionalityReductionMethod,'PCA')
    [~,WF_features_PCAscore,~,~,~]=pca(WF_features,'Algorithm','svd','Centered',false);
    features_reduced=WF_features_PCAscore;
    filename='WF_PCA';
elseif strcmp(dimensionalityReductionMethod,'tSNE')
    tsneOptions=struct;
    tsneOptions.MaxIter=5000;
    tsneOptions.OutputFcn=[];
    tsneOptions.TolFun=1e-10;
    WF_features_tsneScore=tsne(WF_features,'Algorithm','exact','Distance','cosine','NumDimensions',2,'Perplexity',20,'Options',tsneOptions);
    features_reduced=WF_features_tsneScore; 
    filename='WF_tSNE';
end
figure; hold on;
sc=scatter(features_reduced(:,1),features_reduced(:,2),'filled');
sc.Marker='o';
sc.MarkerFaceColor='k';

% K-means clustering
figure('units','normalized','outerposition',[0 0 1 1]);
Klabel = IterateKMeansClustering(WF_features, 1:5, 'NumIterations', 5000, 'Distance','sqeuclidean','Verbose',true);
export_fig(fullfile(folderFigures,strcat(filename,'_clustering_features_',age)),'-tiff','-transparent','-nocrop')
close

colours = {'k', 'r', 'g', 'c', 'm'};
figure('units','normalized','outerposition',[0 0 1 1]);
ax=gca;
% subplot(1,2,1)
hold on
for i=1:numel(unique(Klabel))
    sc=scatter(features_reduced(Klabel==i,1),features_reduced(Klabel==i,2),'filled');
    sc.Marker='o';
    sc.MarkerFaceColor=colours{i};
end
% xlim([0 2.6])
% ylim([0 2.6])
if strcmp(dimensionalityReductionMethod,'PCA')
    ax.XLabel.String='PC 1';
    ax.YLabel.String='PC 2';  
elseif strcmp(dimensionalityReductionMethod,'tSNE')
    ax.XLabel.String='tSNE 1';
    ax.YLabel.String='tSNE 2';
end
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
% ax.YAxis.Visible='off';
grid
export_fig(fullfile(folderFigures,strcat(filename,'_embedding_features_',age)),'-tiff','-transparent','-nocrop')
close


% figure('units','normalized','outerposition',[0 0 1 1]);
figure
subplot(1,2,2)
hold on
ax=gca;
for i=1:numel(unique(Klabel))
    shadedErrorBar(1:size(Waveforms,2),mean(-Waveforms(Klabel==i,:),1),std(Waveforms(Klabel==i,:),1),'lineProps',colours{i})
end
ax.LineWidth = 1.5;
ax.XAxis.Visible='off';
ax.YAxis.Visible='off';
export_fig(fullfile(folderFigures,strcat(filename,'_trace_features_',age)),'-tiff','-transparent','-nocrop')
close

% Replot scatter with cluster identity
if strcmp(age,'P14P18')
    subset_troughPeakTime=data.troughPeakTime(P14P18);
    subset_halfWidth=data.halfWidth(P14P18);
    subset_Nkx=Nkx(P14P18);
    subset_Untagged=Untagged(P14P18);
    subset_SST=SST(P14P18);
elseif strcmp(age,'P8P13')
    subset_troughPeakTime=data.troughPeakTime(P5P8|P9P13);
    subset_halfWidth=data.halfWidth(P5P8|P9P13);
    subset_Nkx=Nkx(P5P8|P9P13);
    subset_Untagged=Untagged(P5P8|P9P13);
    subset_SST=SST(P5P8|P9P13);
end

%P14-P18
figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle(age, 'FontSize',30)
ax=subplot(1,4,1);
[value,edges]=histcounts(subset_troughPeakTime,0:0.1:3);
h=barh(edges(1:end-1),value,'w','LineWidth',2);
set(get(h,'Parent'),'xdir','r')
ylim([0 2.6])
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.XLabel.String='Cell #';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
grid

ax=subplot(1,4,(2:4));
hold on
for i=1:numel(unique(Klabel))
    scatter(subset_halfWidth(Klabel==i&subset_Untagged),subset_troughPeakTime(Klabel==i&subset_Untagged),100, colours{i}, 'filled')
end
for i=1:numel(unique(Klabel))
    scatter(subset_halfWidth(Klabel==i&subset_Nkx),subset_troughPeakTime(Klabel==i&subset_Nkx),500, colours{i}, '*')
end
for i=1:numel(unique(Klabel))
    scatter(subset_halfWidth(Klabel==i&subset_SST),subset_troughPeakTime(Klabel==i&subset_SST),500, colours{i}, '+')
end
xlim([0 1])
ylim([0 2.6])
ax.XLabel.String='Spike Half-Width (ms)';
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
ax.YAxis.Visible='off';
grid
export_fig(fullfile(folderFigures,strcat(filename,'_scatter_features_',age)),'-tiff','-transparent','-nocrop')
close

%% Re-set cell identity
FS=(Untagged & P14P18 & data.troughPeakTime<0.75);
RS=(Untagged & ~FS);
Nkx_FS=(Nkx & P14P18 & data.troughPeakTime<0.75);
Nkx_RS=(Nkx & P14P18 & ~Nkx_FS);
Nkx_y=(Nkx & (P5P8|P9P13));
SST=(data.Tagging=='SST' & data.responseTag==1);

data.cellIdentity(RS)=categorical(cellstr('RS'));
data.cellIdentity(FS)=categorical(cellstr('FS'));
data.cellIdentity(Nkx_FS)=categorical(cellstr('Nkx2-1 - FS'));
data.cellIdentity(Nkx_RS)=categorical(cellstr('Nkx2-1 - RS'));
data.cellIdentity(Nkx_y)=categorical(cellstr('Nkx2-1'));
data.cellIdentity(SST)=categorical(cellstr('SST'));

if ~(nnz(FS)+nnz(RS)+nnz(Nkx_FS)+nnz(Nkx_RS)+nnz(SST)+nnz(Nkx_y)==height(data)); error('Problem with cell identity'); end

%% Spindle burst entrainment RS units
clear ax v
figure('units','normalized','outerposition',[0 0 0.2 .7])
ax(1)=subplot(3,2,1);
v(1,:)=violinplot(data.sb_spikeProb(RS&P9P13&SB_entrained),data.Layer(RS&P9P13&SB_entrained));
ax(1).YLim=[0, 1];
ax(1).YLabel.String = 'SB spike probability';

ax(2)=subplot(3,2,2);
v(2,:)=violinplot(data.sb_spikeProb(RS&P14P18&SB_entrained),data.Layer(RS&P14P18&SB_entrained));
ax(2).YLim=[0, 1];
ax(2).YLabel.String = 'SB spike probability';

% 
% swtest(data.sb_spikeProb(RS&SB_entrained&P9P13))
% [p,tbl,stats]=kruskalwallis(data.sb_spikeProb(RS&SB_entrained&P9P13),data.Layer(RS&SB_entrained&P9P13));
% multcompare(stats,'Dimension',[1,2],'CType' ,'dunn-sidak')

ax(3)=subplot(3,2,3);
v(3,:)=violinplot(data.PPC(RS&P9P13&PPC_entrained,2),data.Layer(RS&P9P13&PPC_entrained));
ax(3).YLim=[0, 1];
ax(3).YLabel.String = 'PPC';
ax(3).YScale = 'log';

ax(4)=subplot(3,2,4);
v(4,:)=violinplot(data.PPC(RS&P14P18&PPC_entrained,2),data.Layer(RS&P14P18&PPC_entrained));
ax(4).YLim=[0, .8];
ax(4).YLabel.String = 'PPC';
swtest(data.PPC(RS&PPC_entrained&P9P13,2))
[p,tbl,stats]=kruskalwallis(data.PPC(RS&PPC_entrained&P9P13,2),data.Layer(RS&PPC_entrained&P9P13));
multcompare(stats,'Dimension',[1,2],'CType' ,'dunn-sidak')

ax(5)=subplot(3,2,5);
v(5,:)=violinplot(wrapTo360(rad2deg(data.vectorAngle(RS&P9P13&PPC_entrained,2))),data.Layer(RS&P9P13&PPC_entrained));
ax(5).YLim=[0, 360];
ax(5).YTick = [0,90,180,270,360];
ax(5).YGrid = 'on';
ax(5).YLabel.String = 'Vector Angle';

ax(6)=subplot(3,2,6);
v(6,:)=violinplot(wrapTo360(rad2deg(data.vectorAngle(RS&P14P18&PPC_entrained,2))),data.Layer(RS&P14P18&PPC_entrained));
ax(6).YLim=[0, 360];
ax(6).YTick = [0,90,180,270,360];
ax(6).YGrid = 'on';
ax(6).YLabel.String = 'Vector Angle';

swtest(wrapTo360(rad2deg(data.vectorAngle(RS&PPC_entrained&P9P13,2))))
[p,tbl,stats]=kruskalwallis(wrapTo360(rad2deg(data.vectorAngle(RS&PPC_entrained&P9P13,2))),data.Layer(RS&PPC_entrained&P9P13));
multcompare(stats,'Dimension',[1,2],'CType' ,'dunn-sidak')

for i=1:numel(ax)
    for j=1:size(v,2)
        v(i,j).ViolinColor=[100,100,100]/255;
        v(i,j).ScatterPlot.MarkerFaceColor=[37,37,37]/255;
        v(i,j).ScatterPlot.MarkerFaceAlpha=1;
        v(i,j).ScatterPlot.SizeData=2;
    end

    ax(i).FontSize=12;
    ax(i).LineWidth=1;
    ax(i).XLim=[0.5,3.5];
    
end

print(gcf,'-dpdf','C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\WT Doctoral Programme in Neuroscience\DPhil in Neuroscience\Manuscripts\V1 S1\Figures\Fig. 2\S1BF_PPC_SB_violinRS')
close

clear ax
figure('units','normalized','outerposition',[0 0 0.2 1])

ySB=[nnz(SB_entrained&L23&P9P13)/nnz(L23&P9P13),nnz(SB_entrained&L4&P9P13)/nnz(L4&P9P13),nnz(SB_entrained&L56&P9P13)/nnz(L56&P9P13);...
    nnz(SB_entrained&L23&P14P18)/nnz(L23&P14P18),nnz(SB_entrained&L4&P14P18)/nnz(L4&P14P18),nnz(SB_entrained&L56&P14P18)/nnz(L56&P14P18)];
yPPC=[nnz(PPC_entrained&L23&P9P13)/nnz(L23&P9P13),nnz(PPC_entrained&L4&P9P13)/nnz(L4&P9P13),nnz(PPC_entrained&L56&P9P13)/nnz(L56&P9P13);...
    nnz(PPC_entrained&L23&P14P18)/nnz(L23&P14P18),nnz(PPC_entrained&L4&P14P18)/nnz(L4&P14P18),nnz(PPC_entrained&L56&P14P18)/nnz(L56&P14P18)];
xDev=reordercats(categorical(cellstr({'P9-P13','P14-P18'})),{'P9-P13','P14-P18'});

ax(1)=subplot(3,1,1);
b(1,:)=bar(xDev,ySB);

ax(2)=subplot(3,1,2);
b(2,:)=bar(xDev,yPPC);

for i=1:numel(ax)
    ax(i).YLabel.String='Ratio single units entrained';
    ax(i).Box='off';
    ax(i).LineWidth = 1.5;
    ax(i).FontSize=15;
    legend('L2/3','L4','L5/6','Location','northwest')
    legend('boxoff')
    for j=1:size(b,2)
        b(i,j).LineWidth=1;
    end
end
export_fig(fullfile(folderFigures,'4.6','PPC_SB_bar'),'-pdf','-transparent','-nocrop')
close
    

%% Responsive cell over days

for i=6:18
    fastPerc(i)=nnz(responsive(data.Age==i & RS,:))/numel(responsive(data.Age==i& RS,:));
    slowPerc(i)=nnz(data.rw_responsive(data.Age==i& RS,:))/numel(data.rw_responsive(data.Age==i& RS,:));
    SBPerc(i)=nnz(SB_entrained(data.Age==i& RS,:))/numel(SB_entrained(data.Age==i& RS,:));
    SB_spikeProb(i)=mean(data(data.Age==i&SB_entrained,:).sb_spikeProb);
    reliability(i)=nanmean(data(data.Age==i&RS,:).reliabilityVisual);
    reliability_K(i)=nanmean(data(data.Age==i&RS,:).reliabilityVisual_K);
    
    fano(i)=nanmean(data(data.Age==i&RS&responsive,:).fanoVisual);

end

yHist=[fastPerc(6:18)',slowPerc(6:18)'];
% yHistAlt=[SBPerc(6:18)',SB_spikeProb(6:18)'];
% yHistAlt=[reliability(6:18)',reliability_K(6:18)'];

figure('units','normalized','outerposition',[0 0 0.5 0.5]);
b=bar(6:18,yHist);
for j=1:size(b,2)          
    b(1,j).LineWidth=1;
end
ax=gca;
ax.XLabel.String='Postnatal day';
ax.YLabel.String='Responsive/Entrained ratio';
ax.Box='off';
ax.LineWidth = 1.5;
ax.FontSize=15;
legend('Fast response','Slow response')
% legend('Entrained','SpikeProb')
% legend('Control','SalB')

legend('boxoff')
export_fig(fullfile(folderFigures,'4.7','RatioResponsiveUnits'),'-pdf','-transparent','-nocrop')
close


clear ax
figure('units','normalized','outerposition',[0 0 0.2 1])

ySB=[nnz(responsive&L23&P9P13&RS)/nnz(L23&P9P13&RS),nnz(responsive&L4&P9P13&RS)/nnz(L4&P9P13&RS),nnz(responsive&L56&P9P13&RS)/nnz(L56&P9P13&RS);...
    nnz(responsive&L23&P14P18&RS)/nnz(L23&P14P18&RS),nnz(responsive&L4&P14P18&RS)/nnz(L4&P14P18&RS),nnz(responsive&L56&P14P18&RS)/nnz(L56&P14P18&RS)];
yPPC=[nnz(data.rw_responsive&L23&P9P13&RS)/nnz(L23&P9P13&RS),nnz(data.rw_responsive&L4&P9P13&RS)/nnz(L4&P9P13&RS),nnz(data.rw_responsive&L56&P9P13&RS)/nnz(L56&P9P13&RS);...
    nnz(data.rw_responsive&L23&P14P18&RS)/nnz(L23&P14P18&RS),nnz(data.rw_responsive&L4&P14P18&RS)/nnz(L4&P14P18&RS),nnz(data.rw_responsive&L56&P14P18&RS)/nnz(L56&P14P18&RS)];
xDev=reordercats(categorical(cellstr({'P9-P13','P14-P18'})),{'P9-P13','P14-P18'});

ax(1)=subplot(3,1,1);
b(1,:)=bar(xDev,ySB);

ax(2)=subplot(3,1,2);
b(2,:)=bar(xDev,yPPC);

for i=1:numel(ax)
    ax(i).YLabel.String='Ratio single units responsive';
    ax(i).Box='off';
    ax(i).LineWidth = 1.5;
    ax(i).FontSize=15;
    legend('L2/3','L4','L5/6','Location','northwest')
    legend('boxoff')
    for j=1:size(b,2)
        b(i,j).LineWidth=1;
    end
end
export_fig(fullfile(folderFigures,'4.7','responsiveFastSlow_bar'),'-pdf','-transparent','-nocrop')
close
    

figure('units','normalized','outerposition',[0 0 1 1]);
smoothingVal=5;
examplePSTH=[data.PSTHvisual(910,:);data.PSTHvisual(944,:);data.PSTHvisual(1532,:);data.PSTHvisual(120,:)];
titlesExamplePSTH={'Fast','Both','Slow','Non-responsive'};
for i=1:4
    ax=subplot(2,2,i);
    hold on
%     patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(examplePSTH(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,49];
    ax.Title.String=titlesExamplePSTH{i};
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
% sgtitle('SST - P9-P13 - Responsive Fast','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-RS-Examples'),'-pdf','-transparent','-nocrop')
close

%% RS units firing
clear ax v
figure('units','normalized','outerposition',[0 0 0.2 1])
ax(1)=subplot(3,2,1);
v(1,:)=violinplot(data.peakVisualFast_N(RS&P9P13&responsive),data.Layer(RS&P9P13&responsive));
ax(1).YLim=[0, 150];
ax(1).YAxis.Label.String='Max fast PSTH (spike/s)';

ax(2)=subplot(3,2,2);
v(2,:)=violinplot(data.peakVisualFast_N(RS&P14P18&responsive),data.Layer(RS&P14P18&responsive));
ax(2).YLim=[0, 150];

[p,tbl,stats]=anovan(data.peakVisualFast_N(RS&responsive&(P9P13|P14P18)),{data.Layer(RS&responsive&(P9P13|P14P18)),removecats(data.Dev(RS&responsive&(P9P13|P14P18)))},'model','interaction','varnames',{'Cell Layer','Development'});
multcompare(stats,'Dimension',[1])

ax(3)=subplot(3,2,3);
v(3,:)=violinplot(data.fanoVisual(RS&P9P13&responsive),data.Layer(RS&P9P13&responsive));
ax(3).YLim=[-0.15, 1.2];
ax(3).YAxis.Label.String='Fano factor';

ax(4)=subplot(3,2,4);
v(4,:)=violinplot(data.fanoVisual(RS&P14P18&responsive),data.Layer(RS&P14P18&responsive));
ax(4).YLim=[-0.15, 1.2];

[p,tbl,stats]=anovan(data.fanoVisual(RS&responsive&(P9P13|P14P18)),{data.Layer(RS&responsive&(P9P13|P14P18)),removecats(data.Dev(RS&responsive&(P9P13|P14P18)))},'model','interaction','varnames',{'Cell Layer','Development'});
figure
multcompare(stats,'Dimension',[1,2])

ax(5)=subplot(3,2,5);
v(5,:)=violinplot(data.rw_firing_N(RS&P9P13&data.rw_responsive),data.Layer(RS&P9P13&data.rw_responsive));
ax(5).YLim=[0, 10];
ax(5).YAxis.Label.String='Average slow PSTH (spike/s)';

ax(6)=subplot(3,2,6);
v(6,:)=violinplot(data.rw_firing_N(RS&P14P18&data.rw_responsive),data.Layer(RS&P14P18&data.rw_responsive));
ax(6).YLim=[0, 10];

[p,tbl,stats]=anovan(data.rw_firing_N(RS&data.rw_responsive&(P9P13|P14P18)),{data.Layer(RS&data.rw_responsive&(P9P13|P14P18)),removecats(data.Dev(RS&data.rw_responsive&(P9P13|P14P18)))},'model','interaction','varnames',{'Cell Layer','Development'});
figure
multcompare(stats,'Dimension',[1,2])

for i=1:numel(ax)
    for j=1:size(v,2)
        v(i,j).ViolinColor=[100,100,100]/255;
        v(i,j).ScatterPlot.MarkerFaceColor=[37,37,37]/255;
        v(i,j).ScatterPlot.MarkerFaceAlpha=1;
        v(i,j).ScatterPlot.SizeData=2;
    end

    ax(i).FontSize=12;
    ax(i).LineWidth=1;
    ax(i).XLim=[0.5,3.5];
    
end
print(gcf,'-dpdf',fullfile(folderFigures,'4.7','ViolinPlots'))
close 

% figure('units','normalized','outerposition',[0 0 0.2 1])
% ax(1)=subplot(3,1,1);
% hold on
% errorbar(5:5:80,mean(data.rw_rho(RS&P9P13&data.rw_responsive,:),1),std(data.rw_rho(RS&P9P13&data.rw_responsive,:),1))
% errorbar(5:5:80,mean(data.rw_rho(RS&P14P18&data.rw_responsive,:),1),std(data.rw_rho(RS&P14P18&data.rw_responsive,:),1))

%% Responsive units stats
% resp=data.rw_responsive;
resp=responsive;
resp2=data.rw_responsive;

colours={[255,227,145]/255,[0,0,179]/255,[203,24,29]/255,[107,174,214]/255,[253,141,60]/255};

figure('units','normalized','outerposition',[0 0 1 1]);
ax=subplot(3,4,[2,3]);
%All Ages
R1=nnz(resp&P5P8)/nnz(P5P8)*100;
R2=nnz(resp&P9P13)/nnz(P9P13)*100;
R3=nnz(resp&P14P18)/nnz(P14P18)*100;

Y_bar=[R1,R2,R3];
X_bar=categorical(cellstr({'P5-P8','P9-P13','P14-P18'}));
X_bar = reordercats(X_bar,{'P5-P8','P9-P13','P14-P18'});
b=bar(X_bar,Y_bar);
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='Percentage responsive cells';
ax.Title.FontSize=30;
ax.YLim=[0,100];
b.LineWidth=1;
b.FaceColor='w';
text(0.9,R1+8,strcat(int2str(nnz(resp&P5P8)),'/',int2str(nnz(P5P8))),'FontSize',20)
text(1.80,R2+8,strcat(int2str(nnz(resp&P9P13)),'/',int2str(nnz(P9P13))),'FontSize',20)
text(2.80,R3+8,strcat(int2str(nnz(resp&P14P18)),'/',int2str(nnz(P14P18))),'FontSize',20)
box('off')

% P9-P13
ax=subplot(3,4,[5,6,9,10]);
R2_RS=nnz(resp&P9P13&RS)/nnz(P9P13&RS)*100;
R2_SST=nnz(resp&P9P13&SST)/nnz(P9P13&SST)*100;
R2_Nkx=nnz(resp&P9P13&Nkx_y)/nnz(P9P13&Nkx_y)*100;

Y_bar=[R2_RS,R2_SST,R2_Nkx];
X_bar=categorical(cellstr({'RS','SST','Nkx2-1'}));
X_bar = reordercats(X_bar,{'RS','SST','Nkx2-1'});
for i=1:numel(Y_bar)
    hold on
    b=bar(X_bar(i),Y_bar(i));
    b.LineWidth=1;
    b.FaceColor=colours{i};
end
hold off
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P9-P13';
ax.Title.FontSize=25;
ax.YLim=[0,109];
text(0.80,R2_RS+4,strcat(int2str(nnz(resp&P9P13&RS)),'/',int2str(nnz(P9P13&RS))),'FontSize',20)
text(1.90,R2_SST+4,strcat(int2str(nnz(resp&P9P13&SST)),'/',int2str(nnz(P9P13&SST))),'FontSize',20)
text(2.88,R2_Nkx+4,strcat(int2str(nnz(resp&P9P13&Nkx_y)),'/',int2str(nnz(P9P13&Nkx_y))),'FontSize',20)

% P14 - P18
ax=subplot(3,4,[7,8,11,12]);
R3_RS=nnz(resp&P14P18&RS)/nnz(P14P18&RS)*100;
R3_FS=nnz(resp&P14P18&FS)/nnz(P14P18&FS)*100;
R3_Nkx_FS=nnz(resp&P14P18&Nkx_FS)/nnz(P14P18&Nkx_FS)*100;
R3_SST=nnz(resp&P14P18&SST)/nnz(P14P18&SST)*100;
R3_Nkx_RS=nnz(resp&P14P18&Nkx_RS)/nnz(P14P18&Nkx_RS)*100;

Y_bar=[R3_RS,R3_SST,R3_FS,R3_Nkx_RS,R3_Nkx_FS];
X_bar=categorical(cellstr({'RS','SST','FS','Nkx2-1 - RS','Nkx2-1 - FS'}));
X_bar = reordercats(X_bar,{'RS','SST','FS','Nkx2-1 - RS','Nkx2-1 - FS'});
for i=1:numel(Y_bar)
    hold on
    b=bar(X_bar(i),Y_bar(i));
    b.LineWidth=1;
    b.FaceColor=colours{i};
end
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P14-P18';
ax.Title.FontSize=25;
ax.YLim=[0,109];
text(0.65,R3_RS+4,strcat(int2str(nnz(resp&P14P18&RS)),'/',int2str(nnz(P14P18&RS))),'FontSize',20)
text(1.85,R3_SST+4,strcat(int2str(nnz(resp&P14P18&SST)),'/',int2str(nnz(P14P18&SST))),'FontSize',20)
text(2.78,R3_FS+4,strcat(int2str(nnz(resp&P14P18&FS)),'/',int2str(nnz(P14P18&FS))),'FontSize',20)
text(3.83,R3_Nkx_RS+4,strcat(int2str(nnz(resp&P14P18&Nkx_RS)),'/',int2str(nnz(P14P18&Nkx_RS))),'FontSize',20)
text(4.83,R3_Nkx_FS+4,strcat(int2str(nnz(resp&P14P18&Nkx_FS)),'/',int2str(nnz(P14P18&Nkx_FS))),'FontSize',20)

export_fig(fullfile(folderFigures,'PercentageResponsiveCells'),'-tiff','-transparent','-nocrop')
close

%% Visual response stats
[MaxPSTH_Val,MaxPSTH_Idx]=max(data.PSTHvisual(:,PSTHbins>0&PSTHbins<200),[],2);
MaxPSTH_Lat=PSTHbins(MaxPSTH_Idx+find(PSTHbins>0,1)-1);
[data.MaxPSTH_Val,data.MaxPSTH_Idx]=max(data.PSTHvisual(:,PSTHbins>0&PSTHbins<3000),[],2);

cats_Y=categories(data.cellIdentity(resp&P9P13));

figure('units','normalized','outerposition',[0 0 1 1]);
ax1=subplot(2,2,1);
v=violinplot(MaxPSTH_Val(resp&P9P13), removecats(data.cellIdentity(resp&P9P13)));
for i=1:numel(v)
    v(i).ViolinColor=colours{i};
    v(i).ScatterPlot.MarkerFaceColor=[189,189,189]/255;
    v(i).ScatterPlot.MarkerFaceAlpha=1;
    v(i).ScatterPlot.SizeData=20;
end
ax1.FontSize=18;
ax1.LineWidth=2;
ax1.YLim=[0,100];
ax1.YLabel.String='Max PSTH rate (Hz)';
ax1.Title.String='P9-P13';
ax1.Title.FontSize=23;
ax1.XTickLabelRotation=20;

ax3=subplot(2,2,3);
v=violinplot(MaxPSTH_Lat(resp&P9P13), removecats(data.cellIdentity(resp&P9P13)));
for i=1:numel(v)
    v(i).ViolinColor=colours{i};
    v(i).ScatterPlot.MarkerFaceColor=[189,189,189]/255;
    v(i).ScatterPlot.MarkerFaceAlpha=1;
    v(i).ScatterPlot.SizeData=20;
end
ax3.FontSize=18;
ax3.LineWidth=2;
ax3.YLabel.String='Max PSTH latency (Hz)';
ax3.XTickLabelRotation=20;

ax2=subplot(2,2,2);
v=violinplot(MaxPSTH_Val(resp&P14P18), removecats(data.cellIdentity(resp&P14P18)), 'Width', 0.5);
for i=1:numel(v)
    v(i).ViolinColor=colours{i};
     v(i).ScatterPlot.MarkerFaceColor=[189,189,189]/255;
    v(i).ScatterPlot.MarkerFaceAlpha=1;
    v(i).ScatterPlot.SizeData=20;
end
ax2.FontSize=18;
ax2.LineWidth=2;
ax2.YLim=[0,200];
ax2.YLabel.String='Max PSTH rate (Hz)';
ax2.Title.String='P14-P18';
ax2.Title.FontSize=23;
ax2.XTickLabelRotation=20;

ax4=subplot(2,2,4);
v=violinplot(MaxPSTH_Lat(resp&P14P18), removecats(data.cellIdentity(resp&P14P18)), 'Width', 0.5);
for i=1:numel(v)
    v(i).ViolinColor=colours{i};
    v(i).ScatterPlot.MarkerFaceColor=[189,189,189]/255;
    v(i).ScatterPlot.MarkerFaceAlpha=1;
    v(i).ScatterPlot.SizeData=20;
end
ax4.FontSize=18;
ax4.LineWidth=2;
ax4.YLabel.String='Max PSTH latency (Hz)';
ax4.XTickLabelRotation=20;

% export_fig(fullfile(folderFigures,'ViolinPlotPSTHAmplitudeLatency'),'-tiff','-transparent','-nocrop')
% ax3.YLim=[0,400];
% ax4.YLim=[0,400];
% export_fig(fullfile(folderFigures,'ViolinPlotPSTHAmplitudeLatency_Zoom'),'-tiff','-transparent','-nocrop')
% close

%% PSTH single units
smoothingVal=5;

figure('units','normalized','outerposition',[0 0 1 1]);
ax1=subplot(2,1,1);
patch([0 100 100 0],[-0.8 -0.8 -0.3 -0.3],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot(PSTHbins,smooth(mean(Z_PSTHvisual(resp & P9P13 & RS,:),1),smoothingVal),'Color',[255,255,191]/255,'LineWidth',2)
plot(PSTHbins,smooth(mean(Z_PSTHvisual(resp & P9P13 & SST,:),1),smoothingVal),'Color',colours{2},'LineWidth',2)
plot(PSTHbins,smooth(mean(Z_PSTHvisual(resp & P9P13 & Nkx_y,:),1),smoothingVal),'Color',colours{3},'LineWidth',2)
ax1=applyFont(ax1,1);
ax1.Color=[237,237,237]/255;
ax1.YLim=[-1 25];
ax1.Title.String='P9-P13';
ax1.Title.FontSize=23;
legend('RS','SST','Nkx')

ax2=subplot(2,1,2);
patch([0 100 100 0],[-0.8 -0.8 -0.3 -0.3],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot(PSTHbins,smooth(mean(Z_PSTHvisual(resp & P14P18 & RS,:),1),smoothingVal),'Color',[255,255,191]/255,'LineWidth',2)
plot(PSTHbins,smooth(mean(Z_PSTHvisual(resp & P14P18 & FS,:),1),smoothingVal),'Color',colours{3},'LineWidth',2)
plot(PSTHbins,smooth(mean(Z_PSTHvisual(resp & P14P18 & SST,:),1),smoothingVal),'Color',colours{2},'LineWidth',2)
plot(PSTHbins,smooth(mean(Z_PSTHvisual(resp & P14P18 & Nkx_FS,:),1),smoothingVal),'Color',colours{5},'LineWidth',2)
plot(PSTHbins,smooth(mean(Z_PSTHvisual(resp & P14P18 & Nkx_RS,:),1),smoothingVal),'Color',colours{4},'LineWidth',2)
ax2=applyFont(ax2,1);
ax2.YLim=[-1 50];
ax2.Color=[237,237,237]/255;
ax2.Title.String='P14-P18';
ax1.Title.FontSize=23;
legend('RS','FS','SST','Nkx-FS','Nkx-RS')

export_fig(fullfile(folderFigures,'AveragePSTH'),'-tiff','-transparent','-nocrop')
ax1.XLim=[-50,1000];
ax2.XLim=[-50,1000];
export_fig(fullfile(folderFigures,'AveragePSTH_Zoom'),'-tiff','-transparent','-nocrop')
close

%% PSTH plot all single units
smoothingVal=5;
%SST young
subset=data(resp & SST & P9P13,:);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:height(subset)
    ax=subplot(ceil(sqrt(height(subset))),ceil(sqrt(height(subset))),i);
    hold on
%     patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,30];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('SST - P9-P13 - Responsive Fast','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-SST-P9P13-ResponsiveFast'),'-tiff','-transparent','-nocrop')
close

subset=data(resp2 & SST & P9P13,:);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:height(subset)
    ax=subplot(ceil(sqrt(height(subset))),ceil(sqrt(height(subset))),i);
    hold on
%     patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,30];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('SST - P9-P13 - Responsive Slow','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-SST-P9P13-ResponsiveSlow'),'-tiff','-transparent','-nocrop')
close

subset=data(~resp & ~resp2 & SST & P9P13,:);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:height(subset)
    ax=subplot(ceil(sqrt(height(subset))),ceil(sqrt(height(subset))),i);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,25];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('SST - P9-P13 - Non-responsive','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-SST-P9P13-NonResponsive'),'-tiff','-transparent','-nocrop')
close

%Nkx2-1 young
subset=data(resp & Nkx_y & P9P13,:);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:height(subset)
    ax=subplot(ceil(sqrt(height(subset))),ceil(sqrt(height(subset))),i);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,35];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('Nkx2-1 - P9-P13 - Responsive Fast','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-Nkx21-P9P13-Responsive Fast'),'-tiff','-transparent','-nocrop')
close

subset=data(resp2 & Nkx_y & P9P13,:);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:height(subset)
    ax=subplot(ceil(sqrt(height(subset))),ceil(sqrt(height(subset))),i);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,35];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('Nkx2-1 - P9-P13 - Responsive Slow','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-Nkx21-P9P13-Responsive Slow'),'-tiff','-transparent','-nocrop')
close

subset=data(~resp &~resp2 & Nkx_y & P9P13,:);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:height(subset)
    ax=subplot(ceil(sqrt(height(subset))),ceil(sqrt(height(subset))),i);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,30];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('Nkx2-1 - P9-P13 - Non-responsive','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-Nkx21-P9P13-NonResponsive'),'-tiff','-transparent','-nocrop')
close

%SST old
subset=data(resp & SST & P14P18,:);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:height(subset)
    ax=subplot(ceil(sqrt(height(subset))),ceil(sqrt(height(subset))),i);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,44];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('SST - P14-P18 - Responsive Fast','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-SST-P14P18-Responsive Fast'),'-tiff','-transparent','-nocrop')
close

subset=data(resp2 & SST & P14P18,:);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:height(subset)
    ax=subplot(ceil(sqrt(height(subset))),ceil(sqrt(height(subset))),i);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,44];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('SST - P14-P18 - Responsive Slow','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-SST-P14P18-Responsive Slow'),'-tiff','-transparent','-nocrop')
close

subset=data(~resp & ~resp2 & SST & P14P18,:);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:height(subset)
    ax=subplot(ceil(sqrt(height(subset))),ceil(sqrt(height(subset))),i);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,60];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('SST - P14-P18 - Non-responsive','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-SST-P14P18-NonResponsive'),'-tiff','-transparent','-nocrop')
close

%Nkx2-1-FS old
subset=data(resp & Nkx_FS & P14P18,:);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:height(subset)
    ax=subplot(ceil(sqrt(height(subset))),ceil(sqrt(height(subset))),i);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,80];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('Nkx-FS - P14-P18 - Responsive Fast','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-Nkx_FS-P14P18-ResponsiveFast'),'-tiff','-transparent','-nocrop')
close

subset=data(resp2 & Nkx_FS & P14P18,:);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:height(subset)
    ax=subplot(ceil(sqrt(height(subset))),ceil(sqrt(height(subset))),i);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,80];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('Nkx-FS - P14-P18 - Responsive Slow','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-Nkx_FS-P14P18-ResponsiveSlow'),'-tiff','-transparent','-nocrop')
close

subset=data(~resp & ~resp2 & Nkx_FS & P14P18,:);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:height(subset)
    ax=subplot(ceil(sqrt(height(subset))),ceil(sqrt(height(subset))),i);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,55];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('Nkx-FS - P14-P18 - Non-responsive','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-Nkx_FS-P14P18-NonResponsive'),'-tiff','-transparent','-nocrop')
close

%Nkx2-1-RS old
subset=data(resp & Nkx_RS & P14P18,:);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:height(subset)
    ax=subplot(ceil(sqrt(height(subset))),ceil(sqrt(height(subset))),i);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,80];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('Nkx-RS - P14-P18 - Responsive Fast','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-Nkx_RS-P14P18-ResponsiveFast'),'-tiff','-transparent','-nocrop')
close

subset=data(resp2 & Nkx_RS & P14P18,:);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:height(subset)
    ax=subplot(ceil(sqrt(height(subset))),ceil(sqrt(height(subset))),i);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,80];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('Nkx-RS - P14-P18 - Responsive Slow','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-Nkx_RS-P14P18-ResponsiveSlow'),'-tiff','-transparent','-nocrop')
close

subset=data(~resp &~resp2 & Nkx_RS & P14P18,:);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:height(subset)
    ax=subplot(ceil(sqrt(height(subset))),ceil(sqrt(height(subset))),i);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,55];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('Nkx-RS - P14-P18 - Non-responsive','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-Nkx_RS-P14P18-NonResponsive'),'-tiff','-transparent','-nocrop')
close

%FS old
subset=data(resp & FS & P14P18,:);
c=0;
f=0;
for i=1:height(subset)
    c=c+1;
    if c==1
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    ax=subplot(4,4,c);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,100];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
    if c==16
        f=f+1;
        sgtitle('FS - P14-P18 - Responsive Fast','FontSize',25)
        export_fig(fullfile(folderFigures,strcat('PSTH-FS-P14P18-ResponsiveFast-',int2str(f))),'-tiff','-transparent','-nocrop')
        close
        c=0; 
    end
end
if c~=0
    export_fig(fullfile(folderFigures,strcat('PSTH-FS-P14P18-ResponsiveFast-',int2str(f+1))),'-tiff','-transparent','-nocrop')
    close
end

subset=data(resp2 & FS & P14P18,:);
c=0;
f=0;
for i=1:height(subset)
    c=c+1;
    if c==1
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    ax=subplot(4,4,c);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,100];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
    if c==16
        f=f+1;
        sgtitle('FS - P14-P18 - Responsive Slow','FontSize',25)
        export_fig(fullfile(folderFigures,strcat('PSTH-FS-P14P18-ResponsiveSlow-',int2str(f))),'-tiff','-transparent','-nocrop')
        close
        c=0; 
    end
end
if c~=0
    export_fig(fullfile(folderFigures,strcat('PSTH-FS-P14P18-ResponsiveSlow-',int2str(f+1))),'-tiff','-transparent','-nocrop')
    close
end

subset=data(~resp & ~resp2 & FS & P14P18,:);
c=0;
f=0;
for i=1:height(subset)
    c=c+1;
    if c==1
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    ax=subplot(4,4,c);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,100];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
    if c==16
        f=f+1;
        sgtitle('FS - P14-P18 - Non-Responsive','FontSize',25)
        export_fig(fullfile(folderFigures,strcat('PSTH-FS-P14P18-NonResponsive-',int2str(f))),'-tiff','-transparent','-nocrop')
        close
        c=0; 
    end
end
if c~=0
    export_fig(fullfile(folderFigures,strcat('PSTH-FS-P14P18-NonResponsive-',int2str(f+1))),'-tiff','-transparent','-nocrop')
    close
end

%RS young
subset=data(resp & RS & P9P13,:);
c=0;
f=0;
for i=1:height(subset)
    c=c+1;
    if c==1
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    ax=subplot(4,4,c);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,50];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
    if c==16
        f=f+1;
        sgtitle('RS - P9-P13 - Responsive Fast','FontSize',25)
        export_fig(fullfile(folderFigures,strcat('PSTH-RS-P9P13-ResponsiveFast-',int2str(f))),'-tiff','-transparent','-nocrop')
        close
        c=0; 
    end
end
if c~=0
    sgtitle('RS - P9-P13 - Responsive Fast','FontSize',25)
    export_fig(fullfile(folderFigures,strcat('PSTH-RS-P9P13-ResponsiveFast-',int2str(f+1))),'-tiff','-transparent','-nocrop')
    close
end

subset=data(resp2 & RS & P9P13,:);
c=0;
f=0;
for i=1:height(subset)
    c=c+1;
    if c==1
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    ax=subplot(4,4,c);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,50];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
    if c==16
        f=f+1;
        sgtitle('RS - P9-P13 - Responsive Slow','FontSize',25)
        export_fig(fullfile(folderFigures,strcat('PSTH-RS-P9P13-ResponsiveSlow-',int2str(f))),'-tiff','-transparent','-nocrop')
        close
        c=0; 
    end
end
if c~=0
    sgtitle('RS - P9-P13 - Responsive Slow','FontSize',25)
    export_fig(fullfile(folderFigures,strcat('PSTH-RS-P9P13-ResponsiveSlow-',int2str(f+1))),'-tiff','-transparent','-nocrop')
    close
end

subset=data(~resp & ~resp2 & RS & P9P13,:);
c=0;
f=0;
for i=1:height(subset)
    c=c+1;
    if c==1
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    ax=subplot(4,4,c);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,50];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
    if c==16
        f=f+1;
        sgtitle('RS - P9-P13 - Non-Responsive','FontSize',25)
        export_fig(fullfile(folderFigures,strcat('PSTH-RS-P9P13-NonResponsive-',int2str(f+1))),'-tiff','-transparent','-nocrop')
        close
        c=0; 
    end
end
if c~=0
    sgtitle('RS - P9-P13 - Non-Responsive','FontSize',25)
    export_fig(fullfile(folderFigures,strcat('PSTH-RS-P9P13-NonResponsive-',int2str(f+1))),'-tiff','-transparent','-nocrop')
    close
end

%RS old
subset=data(resp & RS & P14P18,:);
c=0;
f=0;
for i=1:height(subset)
    c=c+1;
    if c==1
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    ax=subplot(4,4,c);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,50];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
    if c==16
        f=f+1;
        sgtitle('RS - P14-P18 - Responsive','FontSize',25)
        export_fig(fullfile(folderFigures,strcat('PSTH-RS-P14P18-Responsive-',int2str(f))),'-tiff','-transparent','-nocrop')
        close
        c=0; 
    end
end
if c~=0
    sgtitle('RS - P14-P18 - Responsive Fast','FontSize',25)
    export_fig(fullfile(folderFigures,strcat('PSTH-RS-P14P18-ResponsiveFast-',int2str(f+1))),'-tiff','-transparent','-nocrop')
    close
end

subset=data(~resp & ~resp2 & RS & P14P18,:);
c=0;
f=0;
for i=1:height(subset)
    c=c+1;
    if c==1
        figure('units','normalized','outerposition',[0 0 1 1]);
    end
    ax=subplot(4,4,c);
    hold on
    patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(subset.PSTHvisual(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,50];
    ax.Title.String=strcat(char(subset.MouseID(i)),' - U',int2str(subset.suid(i)),' - P',int2str(subset.Age(i)),' -',char(subset.Layer(i)));
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
    if c==16
        f=f+1;
        sgtitle('RS - P14-P18 - Non-Responsive','FontSize',25)
        export_fig(fullfile(folderFigures,strcat('PSTH-RS-P14P18-NonResponsive-',int2str(f+1))),'-tiff','-transparent','-nocrop')
        close
        c=0; 
    end
end
if c~=0
    sgtitle('RS - P14-P18 - Non-Responsive','FontSize',25)
    export_fig(fullfile(folderFigures,strcat('PSTH-RS-P14P18-NonResponsive-',int2str(f+1))),'-tiff','-transparent','-nocrop')
    close
end

%% Imagesc plots - responsive only
data.Z_PSTHvisual=movmean(Z_PSTHvisual,10,2);
resp1=resp;
resp=(resp1|resp2);
data.Depth=-data.Depth;
%%
%SST fast
figure('units','normalized','outerposition',[0 0 1 0.5]);
ax=subplot(1,2,1);
imagesc(PSTHbins,1:height(data(resp&SST&P9P13,:)),sortrows(data(resp&SST&P9P13,:),{'Layer','Depth','MaxPSTH_Idx'}).Z_PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp&SST&P9P13,:).Layer=='L2/3'),nnz(data(resp&SST&P9P13,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp&SST&P9P13,:).Layer=='L5/6')),nnz(not(data(resp&SST&P9P13,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-200,400];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P9-P13';
ax.CLim=[-0.5,6];

ax=subplot(1,2,2);
imagesc(PSTHbins,1:height(data(resp&SST&P14P18,:)),sortrows(data(resp&SST&P14P18,:),{'Layer','Depth','MaxPSTH_Idx'}).Z_PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp&SST&P14P18,:).Layer=='L2/3'),nnz(data(resp&SST&P14P18,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp&SST&P14P18,:).Layer=='L5/6')),nnz(not(data(resp&SST&P14P18,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-200,400];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P14-P18';
ax.CLim=[-0.5,6];

sgtitle('SST','FontSize',25)
export_fig(fullfile(folderFigures,'4.9','PSTH_img_SST_fast'),'-tiff','-pdf','-transparent','-nocrop')
close

%SST slow
figure('units','normalized','outerposition',[0 0 1 0.5]);
ax=subplot(1,2,1);
imagesc(PSTHbins,1:height(data(resp&SST&P9P13,:)),sortrows(data(resp&SST&P9P13,:),{'Layer','Depth','MaxPSTH_Idx'}).Z_PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp&SST&P9P13,:).Layer=='L2/3'),nnz(data(resp&SST&P9P13,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp&SST&P9P13,:).Layer=='L5/6')),nnz(not(data(resp&SST&P9P13,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[400,4000];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P9-P13';
ax.CLim=[-0.5,6];

ax=subplot(1,2,2);
imagesc(PSTHbins,1:height(data(resp&SST&P14P18,:)),sortrows(data(resp&SST&P14P18,:),{'Layer','Depth','MaxPSTH_Idx'}).Z_PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp&SST&P14P18,:).Layer=='L2/3'),nnz(data(resp&SST&P14P18,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp&SST&P14P18,:).Layer=='L5/6')),nnz(not(data(resp&SST&P14P18,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[400,4000];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P14-P18';
ax.CLim=[-0.5,6];

sgtitle('SST','FontSize',25)
export_fig(fullfile(folderFigures,'4.9','PSTH_img_SST_slow'),'-tiff','-pdf','-transparent','-nocrop')
close
    
%Nkx fast
figure('units','normalized','outerposition',[0 0 1 0.5]);
ax=subplot(1,2,1);
imagesc(PSTHbins,1:height(data(resp&Nkx&P9P13,:)),sortrows(data(resp&Nkx&P9P13,:),{'Layer','Depth','MaxPSTH_Idx'}).Z_PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp&Nkx&P9P13,:).Layer=='L2/3'),nnz(data(resp&Nkx&P9P13,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp&Nkx&P9P13,:).Layer=='L5/6')),nnz(not(data(resp&Nkx&P9P13,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-200,400];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P9-P13';
ax.CLim=[-0.5,6];

ax=subplot(1,2,2);
imagesc(PSTHbins,1:height(data(resp&Nkx&P14P18,:)),sortrows(data(resp&Nkx&P14P18,:),{'Layer','Depth','MaxPSTH_Idx'}).Z_PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp&Nkx&P14P18,:).Layer=='L2/3'),nnz(data(resp&Nkx&P14P18,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp&Nkx&P14P18,:).Layer=='L5/6')),nnz(not(data(resp&Nkx&P14P18,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-200,400];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P14-P18';
ax.CLim=[-0.5,6];

sgtitle('Nkx2-1','FontSize',25)
export_fig(fullfile(folderFigures,'4.10','PSTH_img_Nkx_fast'),'-tiff','-pdf','-transparent','-nocrop')
close

%Nkx slow
figure('units','normalized','outerposition',[0 0 1 0.5]);
ax=subplot(1,2,1);
imagesc(PSTHbins,1:height(data(resp&Nkx&P9P13,:)),sortrows(data(resp&Nkx&P9P13,:),{'Layer','Depth','MaxPSTH_Idx'}).Z_PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp&Nkx&P9P13,:).Layer=='L2/3'),nnz(data(resp&Nkx&P9P13,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp&Nkx&P9P13,:).Layer=='L5/6')),nnz(not(data(resp&Nkx&P9P13,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[400,4000];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P9-P13';
ax.CLim=[-0.5,6];

ax=subplot(1,2,2);
imagesc(PSTHbins,1:height(data(resp&Nkx&P14P18,:)),sortrows(data(resp&Nkx&P14P18,:),{'Layer','Depth','MaxPSTH_Idx'}).Z_PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp&Nkx&P14P18,:).Layer=='L2/3'),nnz(data(resp&Nkx&P14P18,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp&Nkx&P14P18,:).Layer=='L5/6')),nnz(not(data(resp&Nkx&P14P18,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[400,4000];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P14-P18';
ax.CLim=[-0.5,6];

sgtitle('Nkx2-1','FontSize',25)
export_fig(fullfile(folderFigures,'4.10','PSTH_img_Nkx_slow'),'-tiff','-pdf','-transparent','-nocrop')
close




figure('units','normalized','outerposition',[0 0 1 0.5]);
ax=subplot(1,2,1);
imagesc(PSTHbins,1:height(data(resp&Nkx&P9P13,:)),sortrows(data(resp&Nkx&P9P13,:),{'Layer','Depth','MaxPSTH_Idx'}).PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp&Nkx&P9P13,:).Layer=='L2/3'),nnz(data(resp&Nkx&P9P13,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp&Nkx&P9P13,:).Layer=='L5/6')),nnz(not(data(resp&Nkx&P9P13,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-100,400];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P9-P13';
ax.CLim=[0,20];

ax=subplot(1,2,2);
imagesc(PSTHbins,1:height(data(resp&Nkx&P14P18,:)),sortrows(data(resp&Nkx&P14P18,:),{'Layer','Depth','MaxPSTH_Idx'}).PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp&Nkx&P14P18,:).Layer=='L2/3'),nnz(data(resp&Nkx&P14P18,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp&Nkx&P14P18,:).Layer=='L5/6')),nnz(not(data(resp&Nkx&P14P18,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-100,400];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P14-P18';

sgtitle('Nkx2-1','FontSize',25)
export_fig(fullfile(folderFigures,'4.10','PSTH_img_Nkx2-1'),'-tiff','-pdf','-transparent','-nocrop')
close

%RS
figure('units','normalized','outerposition',[0 0 1 1]);
ax=subplot(1,2,1);
imagesc(PSTHbins,1:height(data(resp1&RS&P9P13,:)),sortrows(data(resp1&RS&P9P13,:),{'Layer','Depth','MaxPSTH_Idx'}).PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp1&RS&P9P13,:).Layer=='L2/3'),nnz(data(resp1&RS&P9P13,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp1&RS&P9P13,:).Layer=='L5/6')),nnz(not(data(resp1&RS&P9P13,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-100,400];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P9-P13';

ax=subplot(1,2,2);
imagesc(PSTHbins,1:height(data(resp1&RS&P14P18,:)),sortrows(data(resp1&RS&P14P18,:),{'Layer','Depth','MaxPSTH_Idx'}).PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp1&RS&P14P18,:).Layer=='L2/3'),nnz(data(resp1&RS&P14P18,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp1&RS&P14P18,:).Layer=='L5/6')),nnz(not(data(resp1&RS&P14P18,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-100,400];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P14-P18';

sgtitle('RS','FontSize',25)
export_fig(fullfile(folderFigures,'4.9','PSTH_img_RS'),'-tiff','-pdf','-transparent','-nocrop')
close

%FS
figure('units','normalized','outerposition',[0 0 1 1]);
ax=subplot(1,2,2);
imagesc(PSTHbins,1:height(data(resp&FS&P14P18,:)),sortrows(data(resp&FS&P14P18,:),{'Layer','Depth','MaxPSTH_Idx'}).PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp&FS&P14P18,:).Layer=='L2/3'),nnz(data(resp&FS&P14P18,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp&FS&P14P18,:).Layer=='L5/6')),nnz(not(data(resp&FS&P14P18,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-100,400];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P14-P18';

sgtitle('FS','FontSize',25)
export_fig(fullfile(folderFigures,'4.9','PSTH_img_FS'),'-tiff','-pdf','-transparent','-nocrop')
close

%% Bar plot percentage IN respnsive
% SST
clear ax
yHistSB=[nnz(SST&resp1&P9P13)/nnz(SST&P9P13),nnz(SST&resp2&P9P13)/nnz(SST&P9P13);nnz(SST&resp1&P14P18)/nnz(SST&P14P18),nnz(SST&resp2&P14P18)/nnz(SST&P14P18)];
xHistDev=reordercats(categorical(cellstr({'P9-P13','P14-P18'})),{'P9-P13','P14-P18'});
figure('units','normalized','outerposition',[0 0 0.2 1]);
ax=subplot(3,1,1);
b=bar(xHistDev,yHistSB);
b(1,1).LineWidth=1;
b(1,2).LineWidth=1;
ax=gca;
% ax.XLabel.String='Postnatal day';
ax.YLabel.String='Ratio responsive single units';
ax.Box='off';
ax.LineWidth = 1.5;
ax.FontSize=10;
ax.YLim=[0,1];
legend('Fast response','Slow response','Location','northwest')
legend('boxoff')
text(0.76,yHistSB(1,1)+0.05,strcat(int2str(nnz(SST&resp1&P9P13)),'/',int2str(nnz(SST&P9P13))),'FontSize',10)
text(1.03,yHistSB(1,2)+0.05,strcat(int2str(nnz(SST&resp2&P9P13)),'/',int2str(nnz(SST&P9P13))),'FontSize',10)
text(1.73,yHistSB(2,1)+0.05,strcat(int2str(nnz(SST&resp1&P14P18)),'/',int2str(nnz(SST&P14P18))),'FontSize',10)
text(2.02,yHistSB(2,2)+0.05,strcat(int2str(nnz(SST&resp2&P14P18)),'/',int2str(nnz(SST&P14P18))),'FontSize',10)

% export_fig(fullfile(folderFigures,'RatioResponsiveUnits_SST'),'-pdf','-transparent','-nocrop')
% close

% Nkx2-1
yHistNkx=[nnz(Nkx&resp1&P9P13)/nnz(Nkx&P9P13),nnz(Nkx&resp2&P9P13)/nnz(Nkx&P9P13);nnz(Nkx&resp1&P14P18)/nnz(Nkx&P14P18),nnz(Nkx&resp2&P14P18)/nnz(Nkx&P14P18)];
xHistDev=reordercats(categorical(cellstr({'P9-P13','P14-P18'})),{'P9-P13','P14-P18'});
ax=subplot(3,1,2);
b=bar(xHistDev,yHistNkx);
b(1,1).LineWidth=1;
b(1,2).LineWidth=1;
ax=gca;
% ax.XLabel.String='Postnatal day';
ax.YLabel.String='Ratio responsive single units';
ax.Box='off';
ax.LineWidth = 1.5;
ax.FontSize=10;
ax.YLim=[0,1];
% legend('Fast response','Slow response','Location','northwest')
% legend('boxoff')
text(0.76,yHistNkx(1,1)+0.05,strcat(int2str(nnz(Nkx&resp1&P9P13)),'/',int2str(nnz(Nkx&P9P13))),'FontSize',10)
text(1.03,yHistNkx(1,2)+0.05,strcat(int2str(nnz(Nkx&resp2&P9P13)),'/',int2str(nnz(Nkx&P9P13))),'FontSize',10)
text(1.73,yHistNkx(2,1)+0.05,strcat(int2str(nnz(Nkx&resp1&P14P18)),'/',int2str(nnz(Nkx&P14P18))),'FontSize',10)
text(2.02,yHistNkx(2,2)+0.05,strcat(int2str(nnz(Nkx&resp2&P14P18)),'/',int2str(nnz(Nkx&P14P18))),'FontSize',10)

export_fig(fullfile(folderFigures,'4.9','RatioResponsiveUnits_SST&Nkx2-1'),'-pdf','-transparent','-nocrop')
close

%% Response SST and Nkx cells to visual stim
figure('units','normalized','outerposition',[0 0 0.2 1])
clear ax v
ax(1)=subplot(3,2,1);
v(1,:)=violinplot(data.peakVisualFast_N(SST&responsive),data.Dev(SST&responsive));
% ax(1).YLim=[0, 150];
ax(1).YAxis.Label.String='Max fast PSTH (spike/s)';
[h,p,ci,stats] = ttest2(data.peakVisualFast_N(SST&responsive&P9P13),data.peakVisualFast_N(SST&responsive&P14P18));

ax(2)=subplot(3,2,2);
v(2,:)=violinplot(data.rw_firing_N(SST&data.rw_responsive),data.Dev(SST&data.rw_responsive));
% ax(2).YLim=[0, 150];
[h,p,ci,stats] = ttest2(data.rw_firing_N(SST&data.rw_responsive&P9P13),data.rw_firing_N(SST&data.rw_responsive&P14P18));

ax(3)=subplot(3,2,3);
v(3,:)=violinplot(data.fanoVisual(SST&responsive),data.Dev(SST&responsive));
ax(3).YLim=[-0.15, 1.2];
ax(3).YAxis.Label.String='Fano factor';
[h,p,ci,stats] = ttest2(data.fanoVisual(SST&responsive&P9P13),data.peakVisualFast_N(SST&responsive&P14P18));

for i=1:numel(ax)
    for j=1:size(v,2)
        v(i,j).ViolinColor=[100,100,100]/255;
        v(i,j).ScatterPlot.MarkerFaceColor=[37,37,37]/255;
        v(i,j).ScatterPlot.MarkerFaceAlpha=1;
        v(i,j).ScatterPlot.SizeData=5;
    end

    ax(i).FontSize=12;
    ax(i).LineWidth=1;
    ax(i).XLim=[0.5,2.5];
    
end
print(gcf,'-dpdf',fullfile(folderFigures,'4.9','ViolinPlots'))

figure('units','normalized','outerposition',[0 0 0.2 1])
clear ax v
ax(1)=subplot(3,2,1);
v(1,:)=violinplot(data.peakVisualFast_N(Nkx&responsive),data.Dev(Nkx&responsive));
% ax(1).YLim=[0, 150];
ax(1).YAxis.Label.String='Max fast PSTH (spike/s)';
[h,p,ci,stats] = ttest2(data.peakVisualFast_N(Nkx&responsive&P9P13),data.peakVisualFast_N(Nkx&responsive&P14P18));

ax(2)=subplot(3,2,2);
v(2,:)=violinplot(data.rw_firing_N(Nkx&data.rw_responsive),data.Dev(Nkx&data.rw_responsive));
% ax(2).YLim=[0, 150];
[h,p,ci,stats] = ttest2(data.rw_firing_N(Nkx&data.rw_responsive&P9P13),data.rw_firing_N(Nkx&data.rw_responsive&P14P18));

ax(3)=subplot(3,2,3);
v(3,:)=violinplot(data.fanoVisual(Nkx&responsive),data.Dev(Nkx&responsive));
ax(3).YLim=[-0.15, 1.2];
ax(3).YAxis.Label.String='Fano factor';
[h,p,ci,stats] = ttest2(data.fanoVisual(Nkx&responsive&P9P13),data.peakVisualFast_N(Nkx&responsive&P14P18));

for i=1:numel(ax)
    for j=1:size(v,2)
        v(i,j).ViolinColor=[100,100,100]/255;
        v(i,j).ScatterPlot.MarkerFaceColor=[37,37,37]/255;
        v(i,j).ScatterPlot.MarkerFaceAlpha=1;
        v(i,j).ScatterPlot.SizeData=5;
    end

    ax(i).FontSize=12;
    ax(i).LineWidth=1;
    ax(i).XLim=[0.5,2.5];
    
end
print(gcf,'-dpdf',fullfile(folderFigures,'4.10','ViolinPlots'))

%% Laser only SST
% data.PSTHlaser = movmean(data.PSTHlaser,10,2);
% data.PSTHlaser =data.PSTHoptotagging;

%Imagesc plot 
figure('units','normalized','outerposition',[0 0 1 0.5]);
ax=subplot(1,2,1);
imagesc(PSTHbins+50,1:height(data(RS&SSTopto&P9P13&isLaser,:)),zscore(sortrows(data(RS&SSTopto&P9P13&isLaser,:),{'Layer','Depth','meanOpto'}).PSTHlaser,[],'all'))
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(RS&SSTopto&P9P13&isLaser,:).Layer=='L2/3'),nnz(data(RS&SSTopto&P9P13&isLaser,:).Layer=='L2/3')]+0.5,'w--','LineWidth',1)
plot([-1000,5000],[nnz(not(data(RS&SSTopto&P9P13&isLaser,:).Layer=='L5/6')),nnz(not(data(RS&SSTopto&P9P13&isLaser,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',1)
plot([0,0],[1,1000],'w-')
plot([200,200],[1,1000],'w-')
ax.XLim=[-450,+550];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=12;
ax.Title.String='P9-P13';
ax.CLim=[-0.5,8];

ax=subplot(1,2,2);
imagesc(PSTHbins+50,1:height(data(RS&SSTopto&P14P18&isLaser,:)),zscore(sortrows(data(RS&SSTopto&P14P18&isLaser,:),{'Layer','Depth','meanOpto'}).PSTHlaser,[],'all'))
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(RS&SSTopto&P14P18&isLaser,:).Layer=='L2/3'),nnz(data(RS&SSTopto&P14P18&isLaser,:).Layer=='L2/3')]+0.5,'w--','LineWidth',1)
plot([-1000,5000],[nnz(not(data(RS&SSTopto&P14P18&isLaser,:).Layer=='L5/6')),nnz(not(data(RS&SSTopto&P14P18&isLaser,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',1)
plot([0,0],[1,1000],'w-')
plot([200,200],[1,1000],'w-')
ax.XLim=[-450,+550];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=12;
ax.Title.String='P14-P18';
ax.CLim=[-0.5,8];
export_fig(fullfile(folderFigures,'4.11','Heatmaps'),'-pdf','-transparent','-nocrop')
close

%Bar plot
figure('units','normalized','outerposition',[0 0 0.2 1]);
ax=subplot(3,1,1);
barLaserOnlySST = [nnz(RS&P9P13&L23&SSTopto&LaserInhibited),nnz(RS&P9P13&L23&SSTopto&~(LaserInhibited|LaserExcited)),nnz(RS&P9P13&L23&SSTopto&LaserExcited);...
    nnz(RS&P9P13&L4&SSTopto&LaserInhibited),nnz(RS&P9P13&L4&SSTopto&~(LaserInhibited|LaserExcited)),nnz(RS&P9P13&L4&SSTopto&LaserExcited);...
    nnz(RS&P9P13&L56&SSTopto&LaserInhibited),nnz(RS&P9P13&L56&SSTopto&~(LaserInhibited|LaserExcited)),nnz(RS&P9P13&L56&SSTopto&LaserExcited);...
    nnz(RS&P14P18&L23&SSTopto&LaserInhibited),nnz(RS&P14P18&L23&SSTopto&~(LaserInhibited|LaserExcited)),nnz(RS&P14P18&L23&SSTopto&LaserExcited);...
    nnz(RS&P14P18&L4&SSTopto&LaserInhibited),nnz(RS&P14P18&L4&SSTopto&~(LaserInhibited|LaserExcited)),nnz(RS&P14P18&L4&SSTopto&LaserExcited);...
    nnz(RS&P14P18&L56&SSTopto&LaserInhibited),nnz(RS&P14P18&L56&SSTopto&~(LaserInhibited|LaserExcited)),nnz(RS&P14P18&L56&SSTopto&LaserExcited);...
    nnz(FS&P14P18&SSTopto&LaserInhibited),nnz(FS&P14P18&SSTopto&~(LaserInhibited|LaserExcited)),nnz(FS&P14P18&SSTopto&LaserExcited)];
nBarLaserOnlySST = [nnz(RS&P9P13&L23&SSTopto);nnz(RS&P9P13&L4&SSTopto);nnz(RS&P9P13&L56&SSTopto);nnz(RS&P14P18&L23&SSTopto);nnz(RS&P14P18&L4&SSTopto);nnz(RS&P14P18&L56&SSTopto);nnz(FS&P14P18&SSTopto)];
b=bar([1,2,3,5,6,7,9],barLaserOnlySST./nBarLaserOnlySST,'stacked');
for i=1:size(b,2)
    b(1,i).LineWidth=1;
end
ax.YLabel.String='Ratio single units';
ax.Box='off';
ax.LineWidth = 1;
ax.FontSize=10;
ax.YLim=[0,1];
ax.Title.String = '200ms';
legend('Inhibited','Non-responsive','Excited')
legend('boxoff')


ax=subplot(3,1,2);
barLaserOnlySST = [nnz(RS&P9P13&L23&SSTopto&OptoInhibited),nnz(RS&P9P13&L23&SSTopto&~(OptoInhibited|OptoExcited)),nnz(RS&P9P13&L23&SSTopto&OptoExcited);...
    nnz(RS&P9P13&L4&SSTopto&OptoInhibited),nnz(RS&P9P13&L4&SSTopto&~(OptoInhibited|OptoExcited)),nnz(RS&P9P13&L4&SSTopto&OptoExcited);...
    nnz(RS&P9P13&L56&SSTopto&OptoInhibited),nnz(RS&P9P13&L56&SSTopto&~(OptoInhibited|OptoExcited)),nnz(RS&P9P13&L56&SSTopto&OptoExcited);...
    nnz(RS&P14P18&L23&SSTopto&OptoInhibited),nnz(RS&P14P18&L23&SSTopto&~(OptoInhibited|OptoExcited)),nnz(RS&P14P18&L23&SSTopto&OptoExcited);...
    nnz(RS&P14P18&L4&SSTopto&OptoInhibited),nnz(RS&P14P18&L4&SSTopto&~(OptoInhibited|OptoExcited)),nnz(RS&P14P18&L4&SSTopto&OptoExcited);...
    nnz(RS&P14P18&L56&SSTopto&OptoInhibited),nnz(RS&P14P18&L56&SSTopto&~(OptoInhibited|OptoExcited)),nnz(RS&P14P18&L56&SSTopto&OptoExcited);...
    nnz(FS&P14P18&SSTopto&OptoInhibited),nnz(FS&P14P18&SSTopto&~(OptoInhibited|OptoExcited)),nnz(FS&P14P18&SSTopto&OptoExcited)];
nBarLaserOnlySST = [nnz(RS&P9P13&L23&SSTopto);nnz(RS&P9P13&L4&SSTopto);nnz(RS&P9P13&L56&SSTopto);nnz(RS&P14P18&L23&SSTopto);nnz(RS&P14P18&L4&SSTopto);nnz(RS&P14P18&L56&SSTopto);nnz(FS&P14P18&SSTopto)];
b=bar([1,2,3,5,6,7,9],barLaserOnlySST./nBarLaserOnlySST,'stacked');
for i=1:size(b,2)
    b(1,i).LineWidth=1;
end
ax.Title.String = '50ms';
ax.YLabel.String='Ratio single units';
ax.Box='off';
ax.LineWidth = 1;
ax.FontSize=10;
ax.YLim=[0,1];
legend('Inhibited','Non-responsive','Excited')
legend('boxoff')
export_fig(fullfile(folderFigures,'4.11','BarPlot'),'-pdf','-transparent','-nocrop')
close




figure('units','normalized','outerposition',[0 0 0.2 1]);
ax(1)=subplot(3,2,1);
hold on
plot([-3,-2],[data.optoBaseline(RS&P9P13&L23&SSTopto),data.meanOpto(RS&P9P13&L23&SSTopto)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.optoBaseline(RS&P9P13&L23&SSTopto),'omitnan'),mean(data.meanOpto(RS&P9P13&L23&SSTopto),'omitnan')],[sem(data.optoBaseline(RS&P9P13&L23&SSTopto)),sem(data.meanOpto(RS&P9P13&L23&SSTopto))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.optoBaseline(RS&P14P18&L23&SSTopto),data.meanOpto(RS&P14P18&L23&SSTopto)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.optoBaseline(RS&P14P18&L23&SSTopto),'omitnan'),mean(data.meanOpto(RS&P14P18&L23&SSTopto),'omitnan')],[sem(data.optoBaseline(RS&P14P18&L23&SSTopto)),sem(data.meanOpto(RS&P14P18&L23&SSTopto))],'k-o','Linewidth',1,'CapSize',10)

ax(2)=subplot(3,2,2);
hold on
plot([-3,-2],[data.optoBaseline(RS&P9P13&L4&SSTopto),data.meanOpto(RS&P9P13&L4&SSTopto)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.optoBaseline(RS&P9P13&L4&SSTopto),'omitnan'),mean(data.meanOpto(RS&P9P13&L4&SSTopto),'omitnan')],[sem(data.optoBaseline(RS&P9P13&L4&SSTopto)),sem(data.meanOpto(RS&P9P13&L4&SSTopto))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.optoBaseline(RS&P14P18&L4&SSTopto),data.meanOpto(RS&P14P18&L4&SSTopto)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.optoBaseline(RS&P14P18&L4&SSTopto),'omitnan'),mean(data.meanOpto(RS&P14P18&L4&SSTopto),'omitnan')],[sem(data.optoBaseline(RS&P14P18&L4&SSTopto)),sem(data.meanOpto(RS&P14P18&L4&SSTopto))],'k-o','Linewidth',1,'CapSize',10)

ax(3)=subplot(3,2,3);
hold on
plot([-3,-2],[data.optoBaseline(RS&P9P13&L56&SSTopto),data.meanOpto(RS&P9P13&L56&SSTopto)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.optoBaseline(RS&P9P13&L56&SSTopto),'omitnan'),mean(data.meanOpto(RS&P9P13&L56&SSTopto),'omitnan')],[sem(data.optoBaseline(RS&P9P13&L56&SSTopto)),sem(data.meanOpto(RS&P9P13&L56&SSTopto))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.optoBaseline(RS&P14P18&L56&SSTopto),data.meanOpto(RS&P14P18&L56&SSTopto)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.optoBaseline(RS&P14P18&L56&SSTopto),'omitnan'),mean(data.meanOpto(RS&P14P18&L56&SSTopto),'omitnan')],[sem(data.optoBaseline(RS&P14P18&L56&SSTopto)),sem(data.meanOpto(RS&P14P18&L56&SSTopto))],'k-o','Linewidth',1,'CapSize',10)

ax(4)=subplot(3,2,4);
hold on
plot([-1,0],[data.optoBaseline(FS&P14P18&SSTopto),data.meanOpto(FS&P14P18&SSTopto)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.optoBaseline(FS&P14P18&SSTopto),'omitnan'),mean(data.meanOpto(FS&P14P18&SSTopto),'omitnan')],[sem(data.optoBaseline(FS&P14P18&SSTopto)),sem(data.meanOpto(FS&P14P18&SSTopto))],'k-o','Linewidth',1,'CapSize',10)

[h,p,stats] = my_ttest(data.optoBaseline(RS&P9P13&L23&SSTopto),data.meanOpto(RS&P9P13&L23&SSTopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P14P18&L23&SSTopto),data.meanOpto(RS&P14P18&L23&SSTopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P9P13&L4&SSTopto),data.meanOpto(RS&P9P13&L4&SSTopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P14P18&L4&SSTopto),data.meanOpto(RS&P14P18&L4&SSTopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P9P13&L56&SSTopto),data.meanOpto(RS&P9P13&L56&SSTopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P14P18&L56&SSTopto),data.meanOpto(RS&P14P18&L56&SSTopto))
[h,p,stats] = my_ttest(data.optoBaseline(FS&SSTopto),data.meanOpto(FS&SSTopto))

ylims=[0,3;0,5;0,10;0,15];
for i=1:numel(ax)
    ax(i).XLim=[-3.5,.5];
    ax(i).YAxis.Label.String='Spike/s';
    ax(i).FontSize=10;
    ax(i).YLim=ylims(i,:);
end

export_fig(fullfile(folderFigures,'4.11','ErrorBars'),'-pdf','-transparent','-nocrop')
close

%% Laser only Nkx
%Imagesc plot 
figure('units','normalized','outerposition',[0 0 1 0.5]);
ax=subplot(1,2,1);
imagesc(PSTHbins+50,1:height(data(RS&Nkxopto&P9P13&isLaser,:)),zscore(sortrows(data(RS&Nkxopto&P9P13&isLaser,:),{'Layer','Depth','meanOpto'}).PSTHlaser,[],'all'))
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(RS&Nkxopto&P9P13&isLaser,:).Layer=='L2/3'),nnz(data(RS&Nkxopto&P9P13&isLaser,:).Layer=='L2/3')]+0.5,'w--','LineWidth',1)
plot([-1000,5000],[nnz(not(data(RS&Nkxopto&P9P13&isLaser,:).Layer=='L5/6')),nnz(not(data(RS&Nkxopto&P9P13&isLaser,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',1)
plot([0,0],[1,1000],'w-')
plot([200,200],[1,1000],'w-')
ax.XLim=[-450,+550];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=12;
ax.Title.String='P9-P13';
ax.CLim=[-0.5,8];

ax=subplot(1,2,2);
imagesc(PSTHbins+50,1:height(data(RS&Nkxopto&P14P18&isLaser,:)),zscore(sortrows(data(RS&Nkxopto&P14P18&isLaser,:),{'Layer','Depth','meanOpto'}).PSTHlaser,[],'all'))
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(RS&Nkxopto&P14P18&isLaser,:).Layer=='L2/3'),nnz(data(RS&Nkxopto&P14P18&isLaser,:).Layer=='L2/3')]+0.5,'w--','LineWidth',1)
plot([-1000,5000],[nnz(not(data(RS&Nkxopto&P14P18&isLaser,:).Layer=='L5/6')),nnz(not(data(RS&Nkxopto&P14P18&isLaser,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',1)
plot([0,0],[1,1000],'w-')
plot([200,200],[1,1000],'w-')
ax.XLim=[-450,+550];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=12;
ax.Title.String='P14-P18';
ax.CLim=[-0.5,8];
export_fig(fullfile(folderFigures,'4.12','Heatmaps'),'-pdf','-transparent','-nocrop')
close

%Bar plot
figure('units','normalized','outerposition',[0 0 0.2 1]);
ax=subplot(3,1,1);
barLaserOnlySST = [nnz(RS&P9P13&L23&Nkxopto&LaserInhibited),nnz(RS&P9P13&L23&Nkxopto&~(LaserInhibited|LaserExcited)),nnz(RS&P9P13&L23&Nkxopto&LaserExcited);...
    nnz(RS&P9P13&L4&Nkxopto&LaserInhibited),nnz(RS&P9P13&L4&Nkxopto&~(LaserInhibited|LaserExcited)),nnz(RS&P9P13&L4&Nkxopto&LaserExcited);...
    nnz(RS&P9P13&L56&Nkxopto&LaserInhibited),nnz(RS&P9P13&L56&Nkxopto&~(LaserInhibited|LaserExcited)),nnz(RS&P9P13&L56&Nkxopto&LaserExcited);...
    nnz(RS&P14P18&L23&Nkxopto&LaserInhibited),nnz(RS&P14P18&L23&Nkxopto&~(LaserInhibited|LaserExcited)),nnz(RS&P14P18&L23&Nkxopto&LaserExcited);...
    nnz(RS&P14P18&L4&Nkxopto&LaserInhibited),nnz(RS&P14P18&L4&Nkxopto&~(LaserInhibited|LaserExcited)),nnz(RS&P14P18&L4&Nkxopto&LaserExcited);...
    nnz(RS&P14P18&L56&Nkxopto&LaserInhibited),nnz(RS&P14P18&L56&Nkxopto&~(LaserInhibited|LaserExcited)),nnz(RS&P14P18&L56&Nkxopto&LaserExcited);...
    nnz(FS&P14P18&Nkxopto&LaserInhibited),nnz(FS&P14P18&Nkxopto&~(LaserInhibited|LaserExcited)),nnz(FS&P14P18&Nkxopto&LaserExcited)];
nBarLaserOnlySST = [nnz(RS&P9P13&L23&Nkxopto);nnz(RS&P9P13&L4&Nkxopto);nnz(RS&P9P13&L56&Nkxopto);nnz(RS&P14P18&L23&Nkxopto);nnz(RS&P14P18&L4&Nkxopto);nnz(RS&P14P18&L56&Nkxopto);nnz(FS&P14P18&Nkxopto)];
b=bar([1,2,3,5,6,7,9],barLaserOnlySST./nBarLaserOnlySST,'stacked');
for i=1:size(b,2)
    b(1,i).LineWidth=1;
end
ax.YLabel.String='Ratio single units';
ax.Box='off';
ax.LineWidth = 1;
ax.FontSize=10;
ax.YLim=[0,1];
ax.Title.String = '200 ms';
legend('Inhibited','Non-responsive','Excited')
legend('boxoff')

ax=subplot(3,1,2);
barLaserOnlySST = [nnz(RS&P9P13&L23&Nkxopto&OptoInhibited),nnz(RS&P9P13&L23&Nkxopto&~(OptoInhibited|OptoExcited)),nnz(RS&P9P13&L23&Nkxopto&OptoExcited);...
    nnz(RS&P9P13&L4&Nkxopto&OptoInhibited),nnz(RS&P9P13&L4&Nkxopto&~(OptoInhibited|OptoExcited)),nnz(RS&P9P13&L4&Nkxopto&OptoExcited);...
    nnz(RS&P9P13&L56&Nkxopto&OptoInhibited),nnz(RS&P9P13&L56&Nkxopto&~(OptoInhibited|OptoExcited)),nnz(RS&P9P13&L56&Nkxopto&OptoExcited);...
    nnz(RS&P14P18&L23&Nkxopto&OptoInhibited),nnz(RS&P14P18&L23&Nkxopto&~(OptoInhibited|OptoExcited)),nnz(RS&P14P18&L23&Nkxopto&OptoExcited);...
    nnz(RS&P14P18&L4&Nkxopto&OptoInhibited),nnz(RS&P14P18&L4&Nkxopto&~(OptoInhibited|OptoExcited)),nnz(RS&P14P18&L4&Nkxopto&OptoExcited);...
    nnz(RS&P14P18&L56&Nkxopto&OptoInhibited),nnz(RS&P14P18&L56&Nkxopto&~(OptoInhibited|OptoExcited)),nnz(RS&P14P18&L56&Nkxopto&OptoExcited);...
    nnz(FS&P14P18&Nkxopto&OptoInhibited),nnz(FS&P14P18&Nkxopto&~(OptoInhibited|OptoExcited)),nnz(FS&P14P18&Nkxopto&OptoExcited)];
nBarLaserOnlySST = [nnz(RS&P9P13&L23&Nkxopto);nnz(RS&P9P13&L4&Nkxopto);nnz(RS&P9P13&L56&Nkxopto);nnz(RS&P14P18&L23&Nkxopto);nnz(RS&P14P18&L4&Nkxopto);nnz(RS&P14P18&L56&Nkxopto);nnz(FS&P14P18&Nkxopto)];
b=bar([1,2,3,5,6,7,9],barLaserOnlySST./nBarLaserOnlySST,'stacked');
for i=1:size(b,2)
    b(1,i).LineWidth=1;
end
ax.YLabel.String='Ratio single units';
ax.Box='off';
ax.LineWidth = 1;
ax.FontSize=10;
ax.YLim=[0,1];
ax.Title.String = '50 ms';
legend('Inhibited','Non-responsive','Excited')
legend('boxoff')
export_fig(fullfile(folderFigures,'4.12','BarPlot'),'-pdf','-transparent','-nocrop')
close





figure('units','normalized','outerposition',[0 0 0.2 1]);
ax(1)=subplot(3,2,1);
hold on
plot([-3,-2],[data.optoBaseline(RS&P9P13&L23&Nkxopto),data.meanOpto(RS&P9P13&L23&Nkxopto)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.optoBaseline(RS&P9P13&L23&Nkxopto),'omitnan'),mean(data.meanOpto(RS&P9P13&L23&Nkxopto),'omitnan')],[sem(data.optoBaseline(RS&P9P13&L23&Nkxopto)),sem(data.meanOpto(RS&P9P13&L23&Nkxopto))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.optoBaseline(RS&P14P18&L23&Nkxopto),data.meanOpto(RS&P14P18&L23&Nkxopto)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.optoBaseline(RS&P14P18&L23&Nkxopto),'omitnan'),mean(data.meanOpto(RS&P14P18&L23&Nkxopto),'omitnan')],[sem(data.optoBaseline(RS&P14P18&L23&Nkxopto)),sem(data.meanOpto(RS&P14P18&L23&Nkxopto))],'k-o','Linewidth',1,'CapSize',10)

ax(2)=subplot(3,2,2);
hold on
plot([-3,-2],[data.optoBaseline(RS&P9P13&L4&Nkxopto),data.meanOpto(RS&P9P13&L4&Nkxopto)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.optoBaseline(RS&P9P13&L4&Nkxopto),'omitnan'),mean(data.meanOpto(RS&P9P13&L4&Nkxopto),'omitnan')],[sem(data.optoBaseline(RS&P9P13&L4&Nkxopto)),sem(data.meanOpto(RS&P9P13&L4&Nkxopto))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.optoBaseline(RS&P14P18&L4&Nkxopto),data.meanOpto(RS&P14P18&L4&Nkxopto)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.optoBaseline(RS&P14P18&L4&Nkxopto),'omitnan'),mean(data.meanOpto(RS&P14P18&L4&Nkxopto),'omitnan')],[sem(data.optoBaseline(RS&P14P18&L4&Nkxopto)),sem(data.meanOpto(RS&P14P18&L4&Nkxopto))],'k-o','Linewidth',1,'CapSize',10)

ax(3)=subplot(3,2,3);
hold on
plot([-3,-2],[data.optoBaseline(RS&P9P13&L56&Nkxopto),data.meanOpto(RS&P9P13&L56&Nkxopto)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.optoBaseline(RS&P9P13&L56&Nkxopto),'omitnan'),mean(data.meanOpto(RS&P9P13&L56&Nkxopto),'omitnan')],[sem(data.optoBaseline(RS&P9P13&L56&Nkxopto)),sem(data.meanOpto(RS&P9P13&L56&Nkxopto))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.optoBaseline(RS&P14P18&L56&Nkxopto),data.meanOpto(RS&P14P18&L56&Nkxopto)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.optoBaseline(RS&P14P18&L56&Nkxopto),'omitnan'),mean(data.meanOpto(RS&P14P18&L56&Nkxopto),'omitnan')],[sem(data.optoBaseline(RS&P14P18&L56&Nkxopto)),sem(data.meanOpto(RS&P14P18&L56&Nkxopto))],'k-o','Linewidth',1,'CapSize',10)

ax(4)=subplot(3,2,4);
hold on
plot([-1,0],[data.optoBaseline(FS&P14P18&Nkxopto),data.meanOpto(FS&P14P18&Nkxopto)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.optoBaseline(FS&P14P18&Nkxopto),'omitnan'),mean(data.meanOpto(FS&P14P18&Nkxopto),'omitnan')],[sem(data.optoBaseline(FS&P14P18&Nkxopto)),sem(data.meanOpto(FS&P14P18&Nkxopto))],'k-o','Linewidth',1,'CapSize',10)

ylims=[0,4;0,10;0,10;0,15];
for i=1:numel(ax)
    ax(i).XLim=[-3.5,.5];
    ax(i).YAxis.Label.String='Spike/s';
    ax(i).FontSize=10;
    ax(i).YLim=ylims(i,:);
end

export_fig(fullfile(folderFigures,'4.12','ErrorBars'),'-pdf','-transparent','-nocrop')
close

[h,p,stats] = my_ttest(data.optoBaseline(RS&P9P13&L23&Nkxopto),data.meanOpto(RS&P9P13&L23&Nkxopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P14P18&L23&Nkxopto),data.meanOpto(RS&P14P18&L23&Nkxopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P9P13&L4&Nkxopto),data.meanOpto(RS&P9P13&L4&Nkxopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P14P18&L4&Nkxopto),data.meanOpto(RS&P14P18&L4&Nkxopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P9P13&L56&Nkxopto),data.meanOpto(RS&P9P13&L56&Nkxopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P14P18&L56&Nkxopto),data.meanOpto(RS&P14P18&L56&Nkxopto))
[h,p,stats] = my_ttest(data.optoBaseline(FS&Nkxopto),data.meanOpto(FS&Nkxopto))

%% Rebound excited laser only 
%Bar plot
figure('units','normalized','outerposition',[0 0 0.2 1]);
ax=subplot(3,1,1);
barLaserOnlySST = [nnz(RS&P9P13&L23&SSTopto&LaserRebound),nnz(RS&P9P13&L23&SSTopto&~LaserRebound);...
    nnz(RS&P9P13&L4&SSTopto&LaserRebound),nnz(RS&P9P13&L4&SSTopto&~LaserRebound);...
    nnz(RS&P9P13&L56&SSTopto&LaserRebound),nnz(RS&P9P13&L56&SSTopto&~LaserRebound);...
    nnz(RS&P14P18&L23&SSTopto&LaserRebound),nnz(RS&P14P18&L23&SSTopto&~LaserRebound);...
    nnz(RS&P14P18&L4&SSTopto&LaserRebound),nnz(RS&P14P18&L4&SSTopto&~LaserRebound);...
    nnz(RS&P14P18&L56&SSTopto&LaserRebound),nnz(RS&P14P18&L56&SSTopto&~LaserRebound);...
    nnz(FS&P14P18&SSTopto&LaserRebound),nnz(FS&P14P18&SSTopto&~LaserRebound)];
nBarLaserOnlySST = [nnz(RS&P9P13&L23&SSTopto);nnz(RS&P9P13&L4&SSTopto);nnz(RS&P9P13&L56&SSTopto);nnz(RS&P14P18&L23&SSTopto);nnz(RS&P14P18&L4&SSTopto);nnz(RS&P14P18&L56&SSTopto);nnz(FS&P14P18&SSTopto)];
b=bar([1,2,3,5,6,7,9],barLaserOnlySST./nBarLaserOnlySST,'stacked');
for i=1:size(b,2)
    b(1,i).LineWidth=1;
end
ax.YLabel.String='Ratio single units';
ax.Box='off';
ax.LineWidth = 1;
ax.FontSize=10;
ax.YLim=[0,1];
legend('Rebound excited','Non-responsive')
legend('boxoff')


%Bar plot
ax=subplot(3,1,2);
barLaserOnlySST = [nnz(RS&P9P13&L23&Nkxopto&LaserRebound),nnz(RS&P9P13&L23&Nkxopto&~LaserRebound);...
    nnz(RS&P9P13&L4&Nkxopto&LaserRebound),nnz(RS&P9P13&L4&Nkxopto&~LaserRebound);...
    nnz(RS&P9P13&L56&Nkxopto&LaserRebound),nnz(RS&P9P13&L56&Nkxopto&~LaserRebound);...
    nnz(RS&P14P18&L23&Nkxopto&LaserRebound),nnz(RS&P14P18&L23&Nkxopto&~LaserRebound);...
    nnz(RS&P14P18&L4&Nkxopto&LaserRebound),nnz(RS&P14P18&L4&Nkxopto&~LaserRebound);...
    nnz(RS&P14P18&L56&Nkxopto&LaserRebound),nnz(RS&P14P18&L56&Nkxopto&~LaserRebound);...
    nnz(FS&P14P18&Nkxopto&LaserRebound),nnz(FS&P14P18&Nkxopto&~LaserRebound)];
nBarLaserOnlySST = [nnz(RS&P9P13&L23&Nkxopto);nnz(RS&P9P13&L4&Nkxopto);nnz(RS&P9P13&L56&Nkxopto);nnz(RS&P14P18&L23&Nkxopto);nnz(RS&P14P18&L4&Nkxopto);nnz(RS&P14P18&L56&Nkxopto);nnz(FS&P14P18&Nkxopto)];
b=bar([1,2,3,5,6,7,9],barLaserOnlySST./nBarLaserOnlySST,'stacked');
for i=1:size(b,2)
    b(1,i).LineWidth=1;
end
ax.YLabel.String='Ratio single units';
ax.Box='off';
ax.LineWidth = 1;
ax.FontSize=10;
ax.YLim=[0,1];
legend('Rebound excited','Non-responsive')
legend('boxoff')

export_fig(fullfile(folderFigures,'4.13','Rebound_BarPlot'),'-pdf','-transparent','-nocrop')
close


figure('units','normalized','outerposition',[0 0 0.2 1]);
ax(1)=subplot(3,2,1);
hold on
plot([-3,-2],[data.optoBaseline(RS&P9P13&L23&SSTopto),data.optoRebound(RS&P9P13&L23&SSTopto)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.optoBaseline(RS&P9P13&L23&SSTopto),'omitnan'),mean(data.optoRebound(RS&P9P13&L23&SSTopto),'omitnan')],[sem(data.optoBaseline(RS&P9P13&L23&SSTopto)),sem(data.optoRebound(RS&P9P13&L23&SSTopto))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.optoBaseline(RS&P14P18&L23&SSTopto),data.optoRebound(RS&P14P18&L23&SSTopto)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.optoBaseline(RS&P14P18&L23&SSTopto),'omitnan'),mean(data.optoRebound(RS&P14P18&L23&SSTopto),'omitnan')],[sem(data.optoBaseline(RS&P14P18&L23&SSTopto)),sem(data.optoRebound(RS&P14P18&L23&SSTopto))],'k-o','Linewidth',1,'CapSize',10)

ax(2)=subplot(3,2,2);
hold on
plot([-3,-2],[data.optoBaseline(RS&P9P13&L4&SSTopto),data.optoRebound(RS&P9P13&L4&SSTopto)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.optoBaseline(RS&P9P13&L4&SSTopto),'omitnan'),mean(data.optoRebound(RS&P9P13&L4&SSTopto),'omitnan')],[sem(data.optoBaseline(RS&P9P13&L4&SSTopto)),sem(data.optoRebound(RS&P9P13&L4&SSTopto))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.optoBaseline(RS&P14P18&L4&SSTopto),data.optoRebound(RS&P14P18&L4&SSTopto)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.optoBaseline(RS&P14P18&L4&SSTopto),'omitnan'),mean(data.optoRebound(RS&P14P18&L4&SSTopto),'omitnan')],[sem(data.optoBaseline(RS&P14P18&L4&SSTopto)),sem(data.optoRebound(RS&P14P18&L4&SSTopto))],'k-o','Linewidth',1,'CapSize',10)

ax(3)=subplot(3,2,3);
hold on
plot([-3,-2],[data.optoBaseline(RS&P9P13&L56&SSTopto),data.optoRebound(RS&P9P13&L56&SSTopto)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.optoBaseline(RS&P9P13&L56&SSTopto),'omitnan'),mean(data.optoRebound(RS&P9P13&L56&SSTopto),'omitnan')],[sem(data.optoBaseline(RS&P9P13&L56&SSTopto)),sem(data.optoRebound(RS&P9P13&L56&SSTopto))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.optoBaseline(RS&P14P18&L56&SSTopto),data.optoRebound(RS&P14P18&L56&SSTopto)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.optoBaseline(RS&P14P18&L56&SSTopto),'omitnan'),mean(data.optoRebound(RS&P14P18&L56&SSTopto),'omitnan')],[sem(data.optoBaseline(RS&P14P18&L56&SSTopto)),sem(data.optoRebound(RS&P14P18&L56&SSTopto))],'k-o','Linewidth',1,'CapSize',10)

ax(4)=subplot(3,2,4);
hold on
plot([-1,0],[data.optoBaseline(FS&P14P18&SSTopto),data.optoRebound(FS&P14P18&SSTopto)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.optoBaseline(FS&P14P18&SSTopto),'omitnan'),mean(data.optoRebound(FS&P14P18&SSTopto),'omitnan')],[sem(data.optoBaseline(FS&P14P18&SSTopto)),sem(data.optoRebound(FS&P14P18&SSTopto))],'k-o','Linewidth',1,'CapSize',10)

ylims=[0,15;0,15;0,15;0,40];
for i=1:numel(ax)
    ax(i).XLim=[-3.5,.5];
    ax(i).YAxis.Label.String='Spike/s';
    ax(i).FontSize=10;
    ax(i).YLim=ylims(i,:);
end

export_fig(fullfile(folderFigures,'4.13','Rebound_ErrorPlotSST'),'-pdf','-transparent','-nocrop')
close
[h,p,stats] = my_ttest(data.optoBaseline(RS&P9P13&L23&SSTopto),data.optoRebound(RS&P9P13&L23&SSTopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P14P18&L23&SSTopto),data.optoRebound(RS&P14P18&L23&SSTopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P9P13&L4&SSTopto),data.optoRebound(RS&P9P13&L4&SSTopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P14P18&L4&SSTopto),data.optoRebound(RS&P14P18&L4&SSTopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P9P13&L56&SSTopto),data.optoRebound(RS&P9P13&L56&SSTopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P14P18&L56&SSTopto),data.optoRebound(RS&P14P18&L56&SSTopto))
[h,p,stats] = my_ttest(data.optoBaseline(FS&SSTopto),data.optoRebound(FS&SSTopto))

figure('units','normalized','outerposition',[0 0 0.2 1]);
ax(1)=subplot(3,2,1);
hold on
plot([-3,-2],[data.optoBaseline(RS&P9P13&L23&Nkxopto),data.optoRebound(RS&P9P13&L23&Nkxopto)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.optoBaseline(RS&P9P13&L23&Nkxopto),'omitnan'),mean(data.optoRebound(RS&P9P13&L23&Nkxopto),'omitnan')],[sem(data.optoBaseline(RS&P9P13&L23&Nkxopto)),sem(data.optoRebound(RS&P9P13&L23&Nkxopto))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.optoBaseline(RS&P14P18&L23&Nkxopto),data.optoRebound(RS&P14P18&L23&Nkxopto)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.optoBaseline(RS&P14P18&L23&Nkxopto),'omitnan'),mean(data.optoRebound(RS&P14P18&L23&Nkxopto),'omitnan')],[sem(data.optoBaseline(RS&P14P18&L23&Nkxopto)),sem(data.optoRebound(RS&P14P18&L23&Nkxopto))],'k-o','Linewidth',1,'CapSize',10)

ax(2)=subplot(3,2,2);
hold on
plot([-3,-2],[data.optoBaseline(RS&P9P13&L4&Nkxopto),data.optoRebound(RS&P9P13&L4&Nkxopto)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.optoBaseline(RS&P9P13&L4&Nkxopto),'omitnan'),mean(data.optoRebound(RS&P9P13&L4&Nkxopto),'omitnan')],[sem(data.optoBaseline(RS&P9P13&L4&Nkxopto)),sem(data.optoRebound(RS&P9P13&L4&Nkxopto))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.optoBaseline(RS&P14P18&L4&Nkxopto),data.optoRebound(RS&P14P18&L4&Nkxopto)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.optoBaseline(RS&P14P18&L4&Nkxopto),'omitnan'),mean(data.optoRebound(RS&P14P18&L4&Nkxopto),'omitnan')],[sem(data.optoBaseline(RS&P14P18&L4&Nkxopto)),sem(data.optoRebound(RS&P14P18&L4&Nkxopto))],'k-o','Linewidth',1,'CapSize',10)

ax(3)=subplot(3,2,3);
hold on
plot([-3,-2],[data.optoBaseline(RS&P9P13&L56&Nkxopto),data.optoRebound(RS&P9P13&L56&Nkxopto)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.optoBaseline(RS&P9P13&L56&Nkxopto),'omitnan'),mean(data.optoRebound(RS&P9P13&L56&Nkxopto),'omitnan')],[sem(data.optoBaseline(RS&P9P13&L56&Nkxopto)),sem(data.optoRebound(RS&P9P13&L56&Nkxopto))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.optoBaseline(RS&P14P18&L56&Nkxopto),data.optoRebound(RS&P14P18&L56&Nkxopto)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.optoBaseline(RS&P14P18&L56&Nkxopto),'omitnan'),mean(data.optoRebound(RS&P14P18&L56&Nkxopto),'omitnan')],[sem(data.optoBaseline(RS&P14P18&L56&Nkxopto)),sem(data.optoRebound(RS&P14P18&L56&Nkxopto))],'k-o','Linewidth',1,'CapSize',10)

ax(4)=subplot(3,2,4);
hold on
plot([-1,0],[data.optoBaseline(FS&P14P18&Nkxopto),data.optoRebound(FS&P14P18&Nkxopto)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.optoBaseline(FS&P14P18&Nkxopto),'omitnan'),mean(data.optoRebound(FS&P14P18&Nkxopto),'omitnan')],[sem(data.optoBaseline(FS&P14P18&Nkxopto)),sem(data.optoRebound(FS&P14P18&Nkxopto))],'k-o','Linewidth',1,'CapSize',10)

ylims=[0,15;0,15;0,15;0,40];
for i=1:numel(ax)
    ax(i).XLim=[-3.5,.5];
    ax(i).YAxis.Label.String='Spike/s';
    ax(i).FontSize=10;
    ax(i).YLim=ylims(i,:);
end
export_fig(fullfile(folderFigures,'4.13','Rebound_ErrorPlotNkx'),'-pdf','-transparent','-nocrop')
close

[h,p,stats] = my_ttest(data.optoBaseline(RS&P9P13&L23&Nkxopto),data.optoRebound(RS&P9P13&L23&Nkxopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P14P18&L23&Nkxopto),data.optoRebound(RS&P14P18&L23&Nkxopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P9P13&L4&Nkxopto),data.optoRebound(RS&P9P13&L4&Nkxopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P14P18&L4&Nkxopto),data.optoRebound(RS&P14P18&L4&Nkxopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P9P13&L56&Nkxopto),data.optoRebound(RS&P9P13&L56&Nkxopto))
[h,p,stats] = my_ttest(data.optoBaseline(RS&P14P18&L56&Nkxopto),data.optoRebound(RS&P14P18&L56&Nkxopto))
[h,p,stats] = my_ttest(data.optoBaseline(FS&Nkxopto),data.optoRebound(FS&Nkxopto))


%% Plot all optotagged units
figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle('SST tagged - P9-P13')
for i=1:nnz(SST&P9P13)
    subplot(5,6,i)
    plot(PSTHbins,data(SST&P9P13,:).PSTHoptotagging(i,:),'k-','LineWidth',0.75)
    xlim([-50,150])
    patch([0 50 50 0],[0 0 max(data(SST&P9P13,:).PSTHoptotagging(i,:)) max(data(SST&P9P13,:).PSTHoptotagging(i,:))],'c','FaceAlpha',.1,'EdgeColor','none')
    title(data(SST&P9P13,:).Layer(i))
end
export_fig(fullfile(folderFigures,'4.11','AllSSTUnits_P9P13'),'-tiff','-transparent','-nocrop')
close

figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle('SST tagged - P14-P18')
for i=1:nnz(SST&P14P18)
    subplot(5,6,i)
    plot(PSTHbins,data(SST&P14P18,:).PSTHoptotagging(i,:),'k-','LineWidth',0.75)
    xlim([-50,150])
    patch([0 50 50 0],[0 0 max(data(SST&P14P18,:).PSTHoptotagging(i,:)) max(data(SST&P14P18,:).PSTHoptotagging(i,:))],'c','FaceAlpha',.1,'EdgeColor','none')
    title(data(SST&P14P18,:).Layer(i))
end
export_fig(fullfile(folderFigures,'4.11','AllSSTUnits_P14P18'),'-tiff','-transparent','-nocrop')
close

figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle('Nkx2-1 tagged - P9-P13')
for i=1:nnz(Nkx&P9P13)
    subplot(5,6,i)
    plot(PSTHbins,data(Nkx&P9P13,:).PSTHoptotagging(i,:),'k-','LineWidth',0.75)
    xlim([-50,150])
    patch([0 50 50 0],[0 0 max(data(Nkx&P9P13,:).PSTHoptotagging(i,:)) max(data(Nkx&P9P13,:).PSTHoptotagging(i,:))],'c','FaceAlpha',.1,'EdgeColor','none')
    title(data(Nkx&P9P13,:).Layer(i))
end
export_fig(fullfile(folderFigures,'4.11','AllNkxUnits_P9P13'),'-tiff','-transparent','-nocrop')
close

figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle('Nkx2-1 tagged - P14-P18')
for i=1:nnz(Nkx&P14P18)
    subplot(5,6,i)
    plot(PSTHbins,data(Nkx&P14P18,:).PSTHoptotagging(i,:),'k-','LineWidth',0.75)
    xlim([-50,150])
    patch([0 50 50 0],[0 0 max(data(Nkx&P14P18,:).PSTHoptotagging(i,:)) max(data(Nkx&P14P18,:).PSTHoptotagging(i,:))],'c','FaceAlpha',.1,'EdgeColor','none')
    title(data(Nkx&P14P18,:).Layer(i))
end
export_fig(fullfile(folderFigures,'4.11','AllNkxUnits_P14P18'),'-tiff','-transparent','-nocrop')
close

%%
clear ax
figure('units','normalized','outerposition',[0 0 0.2 1]);
ax(1)=subplot(3,1,1);
hold on
plot(-3,data.optoChange(RS&P9P13&L23&SSTopto),'o','Color',[150,150,150]/255)
errorbar(-3+0.2,mean(data.optoChange(RS&P9P13&L23&SSTopto),'omitnan'),sem(data.optoChange(RS&P9P13&L23&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P9P13&L23&SSTopto)) 

plot(-2,data.optoChange(RS&P9P13&L4&SSTopto),'o','Color',[150,150,150]/255)
errorbar(-2+0.2,mean(data.optoChange(RS&P9P13&L4&SSTopto),'omitnan'),sem(data.optoChange(RS&P9P13&L4&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P9P13&L4&SSTopto))

plot(-1,data.optoChange(RS&P9P13&L56&SSTopto),'o','Color',[150,150,150]/255)
errorbar(-1+0.2,mean(data.optoChange(RS&P9P13&L56&SSTopto),'omitnan'),sem(data.optoChange(RS&P9P13&L56&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P9P13&L56&SSTopto))

plot(1,data.optoChange(RS&P14P18&L23&SSTopto),'o','Color',[150,150,150]/255)
errorbar(1+0.2,mean(data.optoChange(RS&P14P18&L23&SSTopto),'omitnan'),sem(data.optoChange(RS&P14P18&L23&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P14P18&L23&SSTopto))

plot(2,data.optoChange(RS&P14P18&L4&SSTopto),'o','Color',[150,150,150]/255)
errorbar(2+0.2,mean(data.optoChange(RS&P14P18&L4&SSTopto),'omitnan'),sem(data.optoChange(RS&P14P18&L4&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P14P18&L4&SSTopto))

plot(3,data.optoChange(RS&P14P18&L56&SSTopto),'o','Color',[150,150,150]/255)
errorbar(3+0.2,mean(data.optoChange(RS&P14P18&L56&SSTopto),'omitnan'),sem(data.optoChange(RS&P14P18&L56&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P14P18&L56&SSTopto))

plot(5,data.optoChange(L23&FS&SSTopto),'o','Color',[150,150,150]/255)
errorbar(5+0.2,mean(data.optoChange(L23&FS&SSTopto),'omitnan'),sem(data.optoChange(L23&FS&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(L23&FS&SSTopto))

plot(6,data.optoChange(L4&FS&SSTopto),'o','Color',[150,150,150]/255)
errorbar(6+0.2,mean(data.optoChange(L4&FS&SSTopto),'omitnan'),sem(data.optoChange(L4&FS&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(L4&FS&SSTopto))

plot(7,data.optoChange(L56&FS&SSTopto),'o','Color',[150,150,150]/255)
errorbar(7+0.2,mean(data.optoChange(L56&FS&SSTopto),'omitnan'),sem(data.optoChange(L56&FS&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(L56&FS&SSTopto))


ax(2)=subplot(3,1,2);
hold on
plot(-3,data.reboundChange(RS&P9P13&L23&SSTopto),'o','Color',[150,150,150]/255)
errorbar(-3+0.2,mean(data.reboundChange(RS&P9P13&L23&SSTopto),'omitnan'),sem(data.reboundChange(RS&P9P13&L23&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(RS&P9P13&L23&SSTopto))     

plot(-2,data.reboundChange(RS&P9P13&L4&SSTopto),'o','Color',[150,150,150]/255)
errorbar(-2+0.2,mean(data.reboundChange(RS&P9P13&L4&SSTopto),'omitnan'),sem(data.reboundChange(RS&P9P13&L4&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(RS&P9P13&L4&SSTopto))

plot(-1,data.reboundChange(RS&P9P13&L56&SSTopto),'o','Color',[150,150,150]/255)
errorbar(-1+0.2,mean(data.reboundChange(RS&P9P13&L56&SSTopto),'omitnan'),sem(data.reboundChange(RS&P9P13&L56&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(RS&P9P13&L56&SSTopto))

plot(1,data.reboundChange(RS&P14P18&L23&SSTopto),'o','Color',[150,150,150]/255)
errorbar(1+0.2,mean(data.reboundChange(RS&P14P18&L23&SSTopto),'omitnan'),sem(data.reboundChange(RS&P14P18&L23&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(RS&P14P18&L23&SSTopto))

plot(2,data.reboundChange(RS&P14P18&L4&SSTopto),'o','Color',[150,150,150]/255)
errorbar(2+0.2,mean(data.reboundChange(RS&P14P18&L4&SSTopto),'omitnan'),sem(data.reboundChange(RS&P14P18&L4&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(RS&P14P18&L4&SSTopto))

plot(3,data.reboundChange(RS&P14P18&L56&SSTopto),'o','Color',[150,150,150]/255)
errorbar(3+0.2,mean(data.reboundChange(RS&P14P18&L56&SSTopto),'omitnan'),sem(data.reboundChange(RS&P14P18&L56&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(RS&P14P18&L56&SSTopto))

plot(5,data.reboundChange(L23&FS&SSTopto),'o','Color',[150,150,150]/255)
errorbar(5+0.2,mean(data.reboundChange(L23&FS&SSTopto),'omitnan'),sem(data.reboundChange(L23&FS&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(L23&FS&SSTopto))

plot(6,data.reboundChange(L4&FS&SSTopto),'o','Color',[150,150,150]/255)
errorbar(6+0.2,mean(data.reboundChange(L4&FS&SSTopto),'omitnan'),sem(data.reboundChange(L4&FS&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(L4&FS&SSTopto))

plot(7,data.reboundChange(L56&FS&SSTopto),'o','Color',[150,150,150]/255)
errorbar(7+0.2,mean(data.reboundChange(L56&FS&SSTopto),'omitnan'),sem(data.reboundChange(L56&FS&SSTopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(L56&FS&SSTopto))


for i=1:numel(ax)
    ax(i).XLim=[-3.5,7.5];
    ax(i).YAxis.Label.String='Mean firing change (ratio log2)';
    ax(i).FontSize=10;
    ax(i).YLim=[-3,3];
    ax(i).XTick=[-3,-2,-1,1,2,3,5,6,7];
    ax(i).XTickLabel={'L2/3','L4','L5/6','L2/3','L4','L5/6','L2/3','L4','L5/6'};
end
export_fig(fullfile(folderFigures,'4.11','ChangePlot'),'-pdf','-transparent','-nocrop')
close

figure('units','normalized','outerposition',[0 0 0.2 1]);
ax(1)=subplot(3,1,1);
hold on
plot(-3,data.optoChange(RS&P9P13&L23&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(-3+0.2,mean(data.optoChange(RS&P9P13&L23&Nkxopto),'omitnan'),sem(data.optoChange(RS&P9P13&L23&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P9P13&L23&Nkxopto))   

plot(-2,data.optoChange(RS&P9P13&L4&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(-2+0.2,mean(data.optoChange(RS&P9P13&L4&Nkxopto),'omitnan'),sem(data.optoChange(RS&P9P13&L4&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P9P13&L4&Nkxopto))

plot(-1,data.optoChange(RS&P9P13&L56&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(-1+0.2,mean(data.optoChange(RS&P9P13&L56&Nkxopto),'omitnan'),sem(data.optoChange(RS&P9P13&L56&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P9P13&L56&Nkxopto))

plot(1,data.optoChange(RS&P14P18&L23&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(1+0.2,mean(data.optoChange(RS&P14P18&L23&Nkxopto),'omitnan'),sem(data.optoChange(RS&P14P18&L23&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P14P18&L23&Nkxopto))

plot(2,data.optoChange(RS&P14P18&L4&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(2+0.2,mean(data.optoChange(RS&P14P18&L4&Nkxopto),'omitnan'),sem(data.optoChange(RS&P14P18&L4&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P14P18&L4&Nkxopto))

plot(3,data.optoChange(RS&P14P18&L56&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(3+0.2,mean(data.optoChange(RS&P14P18&L56&Nkxopto),'omitnan'),sem(data.optoChange(RS&P14P18&L56&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P14P18&L56&Nkxopto))

plot(5,data.optoChange(L23&FS&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(5+0.2,mean(data.optoChange(L23&FS&Nkxopto),'omitnan'),sem(data.optoChange(L23&FS&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(L23&FS&Nkxopto))

plot(6,data.optoChange(L4&FS&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(6+0.2,mean(data.optoChange(L4&FS&Nkxopto),'omitnan'),sem(data.optoChange(L4&FS&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(L4&FS&Nkxopto))

plot(7,data.optoChange(L56&FS&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(7+0.2,mean(data.optoChange(L56&FS&Nkxopto),'omitnan'),sem(data.optoChange(L56&FS&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(L56&FS&Nkxopto))


ax(2)=subplot(3,1,2);
hold on
plot(-3,data.reboundChange(RS&P9P13&L23&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(-3+0.2,mean(data.reboundChange(RS&P9P13&L23&Nkxopto),'omitnan'),sem(data.reboundChange(RS&P9P13&L23&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(RS&P9P13&L23&Nkxopto))     

plot(-2,data.reboundChange(RS&P9P13&L4&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(-2+0.2,mean(data.reboundChange(RS&P9P13&L4&Nkxopto),'omitnan'),sem(data.reboundChange(RS&P9P13&L4&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(RS&P9P13&L4&Nkxopto))

plot(-1,data.reboundChange(RS&P9P13&L56&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(-1+0.2,mean(data.reboundChange(RS&P9P13&L56&Nkxopto),'omitnan'),sem(data.reboundChange(RS&P9P13&L56&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(RS&P9P13&L56&Nkxopto))

plot(1,data.reboundChange(RS&P14P18&L23&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(1+0.2,mean(data.reboundChange(RS&P14P18&L23&Nkxopto),'omitnan'),sem(data.reboundChange(RS&P14P18&L23&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(RS&P14P18&L23&Nkxopto))

plot(2,data.reboundChange(RS&P14P18&L4&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(2+0.2,mean(data.reboundChange(RS&P14P18&L4&Nkxopto),'omitnan'),sem(data.reboundChange(RS&P14P18&L4&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(RS&P14P18&L4&Nkxopto))

plot(3,data.reboundChange(RS&P14P18&L56&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(3+0.2,mean(data.reboundChange(RS&P14P18&L56&Nkxopto),'omitnan'),sem(data.reboundChange(RS&P14P18&L56&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(RS&P14P18&L56&Nkxopto))

plot(5,data.reboundChange(L23&FS&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(5+0.2,mean(data.reboundChange(L23&FS&Nkxopto),'omitnan'),sem(data.reboundChange(L23&FS&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(L23&FS&Nkxopto))

plot(6,data.reboundChange(L4&FS&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(6+0.2,mean(data.reboundChange(L4&FS&Nkxopto),'omitnan'),sem(data.reboundChange(L4&FS&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(L4&FS&Nkxopto))

plot(7,data.reboundChange(L56&FS&Nkxopto),'o','Color',[150,150,150]/255)
errorbar(7+0.2,mean(data.reboundChange(L56&FS&Nkxopto),'omitnan'),sem(data.reboundChange(L56&FS&Nkxopto)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.reboundChange(L56&FS&Nkxopto))



for i=1:numel(ax)
    ax(i).XLim=[-3.5,7.5];
    ax(i).YAxis.Label.String='Mean firing change (ratio log2)';
    ax(i).FontSize=10;
    ax(i).YLim=[-3,3];
    ax(i).XTick=[-3,-2,-1,1,2,3,5,6,7];
    ax(i).XTickLabel={'L2/3','L4','L5/6','L2/3','L4','L5/6','L2/3','L4','L5/6'};
end
export_fig(fullfile(folderFigures,'4.12','ChangePlot'),'-pdf','-transparent','-nocrop')
close

%% SST opto on visual stim
figure('units','normalized','outerposition',[0 0 0.2 1]);
ax(1)=subplot(4,2,1);
hold on
plot([-3,-2],[data.peakVisualFast_N(RS&P9P13&L23&SSTopto&resp1),data.peakVisualOptoFast_N(RS&P9P13&L23&SSTopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.peakVisualFast_N(RS&P9P13&L23&SSTopto&resp1),'omitnan'),mean(data.peakVisualOptoFast_N(RS&P9P13&L23&SSTopto&resp1),'omitnan')],[sem(data.peakVisualFast_N(RS&P9P13&L23&SSTopto&resp1)),sem(data.peakVisualOptoFast_N(RS&P9P13&L23&SSTopto&resp1))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.peakVisualFast_N(RS&P14P18&L23&SSTopto&resp1),data.peakVisualOptoFast_N(RS&P14P18&L23&SSTopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.peakVisualFast_N(RS&P14P18&L23&SSTopto&resp1),'omitnan'),mean(data.peakVisualOptoFast_N(RS&P14P18&L23&SSTopto&resp1),'omitnan')],[sem(data.peakVisualFast_N(RS&P14P18&L23&SSTopto&resp1)),sem(data.peakVisualOptoFast_N(RS&P14P18&L23&SSTopto&resp1))],'k-o','Linewidth',1,'CapSize',10)

ax(2)=subplot(4,2,2);
hold on
plot([-3,-2],[data.peakVisualFast_N(RS&P9P13&L4&SSTopto&resp1),data.peakVisualOptoFast_N(RS&P9P13&L4&SSTopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.peakVisualFast_N(RS&P9P13&L4&SSTopto&resp1),'omitnan'),mean(data.peakVisualOptoFast_N(RS&P9P13&L4&SSTopto&resp1),'omitnan')],[sem(data.peakVisualFast_N(RS&P9P13&L4&SSTopto&resp1)),sem(data.peakVisualOptoFast_N(RS&P9P13&L4&SSTopto&resp1))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.peakVisualFast_N(RS&P14P18&L4&SSTopto&resp1),data.peakVisualOptoFast_N(RS&P14P18&L4&SSTopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.peakVisualFast_N(RS&P14P18&L4&SSTopto&resp1),'omitnan'),mean(data.peakVisualOptoFast_N(RS&P14P18&L4&SSTopto&resp1),'omitnan')],[sem(data.peakVisualFast_N(RS&P14P18&L4&SSTopto&resp1)),sem(data.peakVisualOptoFast_N(RS&P14P18&L4&SSTopto&resp1))],'k-o','Linewidth',1,'CapSize',10)

ax(3)=subplot(4,2,3);
hold on
plot([-3,-2],[data.peakVisualFast_N(RS&P9P13&L56&SSTopto&resp1),data.peakVisualOptoFast_N(RS&P9P13&L56&SSTopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.peakVisualFast_N(RS&P9P13&L56&SSTopto&resp1),'omitnan'),mean(data.peakVisualOptoFast_N(RS&P9P13&L56&SSTopto&resp1),'omitnan')],[sem(data.peakVisualFast_N(RS&P9P13&L56&SSTopto&resp1)),sem(data.peakVisualOptoFast_N(RS&P9P13&L56&SSTopto&resp1))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.peakVisualFast_N(RS&P14P18&L56&SSTopto&resp1),data.peakVisualOptoFast_N(RS&P14P18&L56&SSTopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.peakVisualFast_N(RS&P14P18&L56&SSTopto&resp1),'omitnan'),mean(data.peakVisualOptoFast_N(RS&P14P18&L56&SSTopto&resp1),'omitnan')],[sem(data.peakVisualFast_N(RS&P14P18&L56&SSTopto&resp1)),sem(data.peakVisualOptoFast_N(RS&P14P18&L56&SSTopto&resp1))],'k-o','Linewidth',1,'CapSize',10)

ax(4)=subplot(4,2,4);
hold on
plot([-1,0],[data.peakVisualFast_N(FS&P14P18&SSTopto&resp1),data.peakVisualOptoFast_N(FS&P14P18&SSTopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.peakVisualFast_N(FS&P14P18&SSTopto&resp1),'omitnan'),mean(data.peakVisualOptoFast_N(FS&P14P18&SSTopto&resp1),'omitnan')],[sem(data.peakVisualFast_N(FS&P14P18&SSTopto&resp1)),sem(data.peakVisualOptoFast_N(FS&P14P18&SSTopto&resp1))],'k-o','Linewidth',1,'CapSize',10)

ylims=[0,100;0,100;0,100;0,150];
for i=1:numel(ax)
    ax(i).XLim=[-3.5,.5];
    ax(i).YAxis.Label.String='Spike/s';
    ax(i).FontSize=10;
    ax(i).YLim=ylims(i,:);
end

export_fig(fullfile(folderFigures,'4.13','ErrorBars_FastRaw'),'-pdf','-transparent','-nocrop')
close

[h,p,stats]=my_ttest(data.peakVisualFast_N(RS&P9P13&L23&SSTopto&resp1),data.peakVisualOptoFast_N(RS&P9P13&L23&SSTopto&resp1))
[h,p,stats]=my_ttest(data.peakVisualFast_N(RS&P14P18&L23&SSTopto&resp1),data.peakVisualOptoFast_N(RS&P14P18&L23&SSTopto&resp1))
[h,p,stats]=my_ttest(data.peakVisualFast_N(RS&P9P13&L4&SSTopto&resp1),data.peakVisualOptoFast_N(RS&P9P13&L4&SSTopto&resp1))
[h,p,stats]=my_ttest(data.peakVisualFast_N(RS&P14P18&L4&SSTopto&resp1),data.peakVisualOptoFast_N(RS&P14P18&L4&SSTopto&resp1))
[h,p,stats]=my_ttest(data.peakVisualFast_N(RS&P9P13&L56&SSTopto&resp1),data.peakVisualOptoFast_N(RS&P9P13&L56&SSTopto&resp1))
[h,p,stats]=my_ttest(data.peakVisualFast_N(RS&P14P18&L56&SSTopto&resp1),data.peakVisualOptoFast_N(RS&P14P18&L56&SSTopto&resp1))
[h,p,stats]=my_ttest(data.peakVisualFast_N(FS&SSTopto&resp1),data.peakVisualOptoFast_N(FS&SSTopto&resp1))

figure('units','normalized','outerposition',[0 0 0.2 1]);
ax(1)=subplot(4,2,1);
hold on
plot([-3,-2],[data.rw_meanPSTH(RS&P9P13&L23&SSTopto&resp2),data.rwOpto_meanPSTH(RS&P9P13&L23&SSTopto&resp2)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.rw_meanPSTH(RS&P9P13&L23&SSTopto&resp2),'omitnan'),mean(data.rwOpto_meanPSTH(RS&P9P13&L23&SSTopto&resp2),'omitnan')],[sem(data.rw_meanPSTH(RS&P9P13&L23&SSTopto&resp2)),sem(data.rwOpto_meanPSTH(RS&P9P13&L23&SSTopto&resp2))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.rw_meanPSTH(RS&P14P18&L23&SSTopto&resp2),data.rwOpto_meanPSTH(RS&P14P18&L23&SSTopto&resp2)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.rw_meanPSTH(RS&P14P18&L23&SSTopto&resp2),'omitnan'),mean(data.rwOpto_meanPSTH(RS&P14P18&L23&SSTopto&resp2),'omitnan')],[sem(data.rw_meanPSTH(RS&P14P18&L23&SSTopto&resp2)),sem(data.rwOpto_meanPSTH(RS&P14P18&L23&SSTopto&resp2))],'k-o','Linewidth',1,'CapSize',10)

ax(2)=subplot(4,2,2);
hold on
plot([-3,-2],[data.rw_meanPSTH(RS&P9P13&L4&SSTopto&resp2),data.rwOpto_meanPSTH(RS&P9P13&L4&SSTopto&resp2)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.rw_meanPSTH(RS&P9P13&L4&SSTopto&resp2),'omitnan'),mean(data.rwOpto_meanPSTH(RS&P9P13&L4&SSTopto&resp2),'omitnan')],[sem(data.rw_meanPSTH(RS&P9P13&L4&SSTopto&resp2)),sem(data.rwOpto_meanPSTH(RS&P9P13&L4&SSTopto&resp2))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.rw_meanPSTH(RS&P14P18&L4&SSTopto&resp2),data.rwOpto_meanPSTH(RS&P14P18&L4&SSTopto&resp2)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.rw_meanPSTH(RS&P14P18&L4&SSTopto&resp2),'omitnan'),mean(data.rwOpto_meanPSTH(RS&P14P18&L4&SSTopto&resp2),'omitnan')],[sem(data.rw_meanPSTH(RS&P14P18&L4&SSTopto&resp2)),sem(data.rwOpto_meanPSTH(RS&P14P18&L4&SSTopto&resp2))],'k-o','Linewidth',1,'CapSize',10)

ax(3)=subplot(4,2,3);
hold on
plot([-3,-2],[data.rw_meanPSTH(RS&P9P13&L56&SSTopto&resp2),data.rwOpto_meanPSTH(RS&P9P13&L56&SSTopto&resp2)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.rw_meanPSTH(RS&P9P13&L56&SSTopto&resp2),'omitnan'),mean(data.rwOpto_meanPSTH(RS&P9P13&L56&SSTopto&resp2),'omitnan')],[sem(data.rw_meanPSTH(RS&P9P13&L56&SSTopto&resp2)),sem(data.rwOpto_meanPSTH(RS&P9P13&L56&SSTopto&resp2))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.rw_meanPSTH(RS&P14P18&L56&SSTopto&resp2),data.rwOpto_meanPSTH(RS&P14P18&L56&SSTopto&resp2)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.rw_meanPSTH(RS&P14P18&L56&SSTopto&resp2),'omitnan'),mean(data.rwOpto_meanPSTH(RS&P14P18&L56&SSTopto&resp2),'omitnan')],[sem(data.rw_meanPSTH(RS&P14P18&L56&SSTopto&resp2)),sem(data.rwOpto_meanPSTH(RS&P14P18&L56&SSTopto&resp2))],'k-o','Linewidth',1,'CapSize',10)

ax(4)=subplot(4,2,4);
hold on
plot([-1,0],[data.rw_meanPSTH(FS&P14P18&SSTopto&resp2),data.rwOpto_meanPSTH(FS&P14P18&SSTopto&resp2)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.rw_meanPSTH(FS&P14P18&SSTopto&resp2),'omitnan'),mean(data.rwOpto_meanPSTH(FS&P14P18&SSTopto&resp2),'omitnan')],[sem(data.rw_meanPSTH(FS&P14P18&SSTopto&resp2)),sem(data.rwOpto_meanPSTH(FS&P14P18&SSTopto&resp2))],'k-o','Linewidth',1,'CapSize',10)

ylims=[0,15;0,15;0,15;0,25];
for i=1:numel(ax)
    ax(i).XLim=[-3.5,.5];
    ax(i).YAxis.Label.String='Spike/s';
    ax(i).FontSize=10;
    ax(i).YLim=ylims(i,:);
end

export_fig(fullfile(folderFigures,'4.13','ErrorBars_SlowRaw'),'-pdf','-transparent','-nocrop')
close

[h,p,stats]=my_ttest(data.rw_meanPSTH(RS&P9P13&L23&SSTopto&resp2),data.rwOpto_meanPSTH(RS&P9P13&L23&SSTopto&resp2))
[h,p,stats]=my_ttest(data.rw_meanPSTH(RS&P14P18&L23&SSTopto&resp2),data.rwOpto_meanPSTH(RS&P14P18&L23&SSTopto&resp2))
[h,p,stats]=my_ttest(data.rw_meanPSTH(RS&P9P13&L4&SSTopto&resp2),data.rwOpto_meanPSTH(RS&P9P13&L4&SSTopto&resp2))
[h,p,stats]=my_ttest(data.rw_meanPSTH(RS&P14P18&L4&SSTopto&resp2),data.rwOpto_meanPSTH(RS&P14P18&L4&SSTopto&resp2))
[h,p,stats]=my_ttest(data.rw_meanPSTH(RS&P9P13&L56&SSTopto&resp2),data.rwOpto_meanPSTH(RS&P9P13&L56&SSTopto&resp2))
[h,p,stats]=my_ttest(data.rw_meanPSTH(RS&P14P18&L56&SSTopto&resp2),data.rwOpto_meanPSTH(RS&P14P18&L56&SSTopto&resp2))
[h,p,stats]=my_ttest(data.rw_meanPSTH(FS&SSTopto&resp2),data.rwOpto_meanPSTH(FS&SSTopto&resp2))

figure('units','normalized','outerposition',[0 0 0.2 1]);
ax(1)=subplot(4,2,1);
hold on
plot([-3,-2],[data.fanoVisual(RS&P9P13&L23&SSTopto&resp1),data.fanoVisualOpto(RS&P9P13&L23&SSTopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.fanoVisual(RS&P9P13&L23&SSTopto&resp1),'omitnan'),mean(data.fanoVisualOpto(RS&P9P13&L23&SSTopto&resp1),'omitnan')],[sem(data.fanoVisual(RS&P9P13&L23&SSTopto&resp1)),sem(data.fanoVisualOpto(RS&P9P13&L23&SSTopto&resp1))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.fanoVisual(RS&P14P18&L23&SSTopto&resp1),data.fanoVisualOpto(RS&P14P18&L23&SSTopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.fanoVisual(RS&P14P18&L23&SSTopto&resp1),'omitnan'),mean(data.fanoVisualOpto(RS&P14P18&L23&SSTopto&resp1),'omitnan')],[sem(data.fanoVisual(RS&P14P18&L23&SSTopto&resp1)),sem(data.fanoVisualOpto(RS&P14P18&L23&SSTopto&resp1))],'k-o','Linewidth',1,'CapSize',10)

ax(2)=subplot(4,2,2);
hold on
plot([-3,-2],[data.fanoVisual(RS&P9P13&L4&SSTopto&resp1),data.fanoVisualOpto(RS&P9P13&L4&SSTopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.fanoVisual(RS&P9P13&L4&SSTopto&resp1),'omitnan'),mean(data.fanoVisualOpto(RS&P9P13&L4&SSTopto&resp1),'omitnan')],[sem(data.fanoVisual(RS&P9P13&L4&SSTopto&resp1)),sem(data.fanoVisualOpto(RS&P9P13&L4&SSTopto&resp1))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.fanoVisual(RS&P14P18&L4&SSTopto&resp1),data.fanoVisualOpto(RS&P14P18&L4&SSTopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.fanoVisual(RS&P14P18&L4&SSTopto&resp1),'omitnan'),mean(data.fanoVisualOpto(RS&P14P18&L4&SSTopto&resp1),'omitnan')],[sem(data.fanoVisual(RS&P14P18&L4&SSTopto&resp1)),sem(data.fanoVisualOpto(RS&P14P18&L4&SSTopto&resp1))],'k-o','Linewidth',1,'CapSize',10)

ax(3)=subplot(4,2,3);
hold on
plot([-3,-2],[data.fanoVisual(RS&P9P13&L56&SSTopto&resp1),data.fanoVisualOpto(RS&P9P13&L56&SSTopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.fanoVisual(RS&P9P13&L56&SSTopto&resp1),'omitnan'),mean(data.fanoVisualOpto(RS&P9P13&L56&SSTopto&resp1),'omitnan')],[sem(data.fanoVisual(RS&P9P13&L56&SSTopto&resp1)),sem(data.fanoVisualOpto(RS&P9P13&L56&SSTopto&resp1))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.fanoVisual(RS&P14P18&L56&SSTopto&resp1),data.fanoVisualOpto(RS&P14P18&L56&SSTopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.fanoVisual(RS&P14P18&L56&SSTopto&resp1),'omitnan'),mean(data.fanoVisualOpto(RS&P14P18&L56&SSTopto&resp1),'omitnan')],[sem(data.fanoVisual(RS&P14P18&L56&SSTopto&resp1)),sem(data.fanoVisualOpto(RS&P14P18&L56&SSTopto&resp1))],'k-o','Linewidth',1,'CapSize',10)

ax(4)=subplot(4,2,4);
hold on
plot([-1,0],[data.fanoVisual(FS&P14P18&SSTopto&resp1),data.fanoVisualOpto(FS&P14P18&SSTopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.fanoVisual(FS&P14P18&SSTopto&resp1),'omitnan'),mean(data.fanoVisualOpto(FS&P14P18&SSTopto&resp1),'omitnan')],[sem(data.fanoVisual(FS&P14P18&SSTopto&resp1)),sem(data.fanoVisualOpto(FS&P14P18&SSTopto&resp1))],'k-o','Linewidth',1,'CapSize',10)

% ylims=[0,1.5;0,1.5;0,1.5;0,1.5];
for i=1:numel(ax)
    ax(i).XLim=[-3.5,.5];
    ax(i).YAxis.Label.String='Fano factor';
    ax(i).FontSize=10;
    ax(i).YLim=[0, 2.5];
end

export_fig(fullfile(folderFigures,'4.13','FanoFactor'),'-pdf','-transparent','-nocrop')
close

[h,p,stats]=my_ttest(data.fanoVisual(RS&P9P13&L23&SSTopto&resp1),data.fanoVisualOpto(RS&P9P13&L23&SSTopto&resp1))
[h,p,stats]=my_ttest(data.fanoVisual(RS&P14P18&L23&SSTopto&resp1),data.fanoVisualOpto(RS&P14P18&L23&SSTopto&resp1))
[h,p,stats]=my_ttest(data.fanoVisual(RS&P9P13&L4&SSTopto&resp1),data.fanoVisualOpto(RS&P9P13&L4&SSTopto&resp1))
[h,p,stats]=my_ttest(data.fanoVisual(RS&P14P18&L4&SSTopto&resp1),data.fanoVisualOpto(RS&P14P18&L4&SSTopto&resp1))
[h,p,stats]=my_ttest(data.fanoVisual(RS&P9P13&L56&SSTopto&resp1),data.fanoVisualOpto(RS&P9P13&L56&SSTopto&resp1))
[h,p,stats]=my_ttest(data.fanoVisual(RS&P14P18&L56&SSTopto&resp1),data.fanoVisualOpto(RS&P14P18&L56&SSTopto&resp1))
[h,p,stats]=my_ttest(data.fanoVisual(FS&SSTopto&resp1),data.fanoVisualOpto(FS&SSTopto&resp1))

%% Nkx opto on visual stim
figure('units','normalized','outerposition',[0 0 0.2 1]);
ax(1)=subplot(4,2,1);
hold on
plot([-3,-2],[data.peakVisualFast_N(RS&P9P13&L23&Nkxopto&resp1),data.peakVisualOptoFast_N(RS&P9P13&L23&Nkxopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.peakVisualFast_N(RS&P9P13&L23&Nkxopto&resp1),'omitnan'),mean(data.peakVisualOptoFast_N(RS&P9P13&L23&Nkxopto&resp1),'omitnan')],[sem(data.peakVisualFast_N(RS&P9P13&L23&Nkxopto&resp1)),sem(data.peakVisualOptoFast_N(RS&P9P13&L23&Nkxopto&resp1))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.peakVisualFast_N(RS&P14P18&L23&Nkxopto&resp1),data.peakVisualOptoFast_N(RS&P14P18&L23&Nkxopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.peakVisualFast_N(RS&P14P18&L23&Nkxopto&resp1),'omitnan'),mean(data.peakVisualOptoFast_N(RS&P14P18&L23&Nkxopto&resp1),'omitnan')],[sem(data.peakVisualFast_N(RS&P14P18&L23&Nkxopto&resp1)),sem(data.peakVisualOptoFast_N(RS&P14P18&L23&Nkxopto&resp1))],'k-o','Linewidth',1,'CapSize',10)

ax(2)=subplot(4,2,2);
hold on
plot([-3,-2],[data.peakVisualFast_N(RS&P9P13&L4&Nkxopto&resp1),data.peakVisualOptoFast_N(RS&P9P13&L4&Nkxopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.peakVisualFast_N(RS&P9P13&L4&Nkxopto&resp1),'omitnan'),mean(data.peakVisualOptoFast_N(RS&P9P13&L4&Nkxopto&resp1),'omitnan')],[sem(data.peakVisualFast_N(RS&P9P13&L4&Nkxopto&resp1)),sem(data.peakVisualOptoFast_N(RS&P9P13&L4&Nkxopto&resp1))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.peakVisualFast_N(RS&P14P18&L4&Nkxopto&resp1),data.peakVisualOptoFast_N(RS&P14P18&L4&Nkxopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.peakVisualFast_N(RS&P14P18&L4&Nkxopto&resp1),'omitnan'),mean(data.peakVisualOptoFast_N(RS&P14P18&L4&Nkxopto&resp1),'omitnan')],[sem(data.peakVisualFast_N(RS&P14P18&L4&Nkxopto&resp1)),sem(data.peakVisualOptoFast_N(RS&P14P18&L4&Nkxopto&resp1))],'k-o','Linewidth',1,'CapSize',10)

ax(3)=subplot(4,2,3);
hold on
plot([-3,-2],[data.peakVisualFast_N(RS&P9P13&L56&Nkxopto&resp1),data.peakVisualOptoFast_N(RS&P9P13&L56&Nkxopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.peakVisualFast_N(RS&P9P13&L56&Nkxopto&resp1),'omitnan'),mean(data.peakVisualOptoFast_N(RS&P9P13&L56&Nkxopto&resp1),'omitnan')],[sem(data.peakVisualFast_N(RS&P9P13&L56&Nkxopto&resp1)),sem(data.peakVisualOptoFast_N(RS&P9P13&L56&Nkxopto&resp1))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.peakVisualFast_N(RS&P14P18&L56&Nkxopto&resp1),data.peakVisualOptoFast_N(RS&P14P18&L56&Nkxopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.peakVisualFast_N(RS&P14P18&L56&Nkxopto&resp1),'omitnan'),mean(data.peakVisualOptoFast_N(RS&P14P18&L56&Nkxopto&resp1),'omitnan')],[sem(data.peakVisualFast_N(RS&P14P18&L56&Nkxopto&resp1)),sem(data.peakVisualOptoFast_N(RS&P14P18&L56&Nkxopto&resp1))],'k-o','Linewidth',1,'CapSize',10)

ax(4)=subplot(4,2,4);
hold on
plot([-1,0],[data.peakVisualFast_N(FS&P14P18&Nkxopto&resp1),data.peakVisualOptoFast_N(FS&P14P18&Nkxopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.peakVisualFast_N(FS&P14P18&Nkxopto&resp1),'omitnan'),mean(data.peakVisualOptoFast_N(FS&P14P18&Nkxopto&resp1),'omitnan')],[sem(data.peakVisualFast_N(FS&P14P18&Nkxopto&resp1)),sem(data.peakVisualOptoFast_N(FS&P14P18&Nkxopto&resp1))],'k-o','Linewidth',1,'CapSize',10)

ylims=[0,100;0,100;0,100;0,150];
for i=1:numel(ax)
    ax(i).XLim=[-3.5,.5];
    ax(i).YAxis.Label.String='Spike/s';
    ax(i).FontSize=10;
    ax(i).YLim=ylims(i,:);
end

export_fig(fullfile(folderFigures,'4.14','ErrorBars_Raw'),'-pdf','-transparent','-nocrop')
close

[h,p,stats]=my_ttest(data.peakVisualFast_N(RS&P9P13&L23&Nkxopto&resp1),data.peakVisualOptoFast_N(RS&P9P13&L23&Nkxopto&resp1))
[h,p,stats]=my_ttest(data.peakVisualFast_N(RS&P14P18&L23&Nkxopto&resp1),data.peakVisualOptoFast_N(RS&P14P18&L23&Nkxopto&resp1))
[h,p,stats]=my_ttest(data.peakVisualFast_N(RS&P14P18&L4&Nkxopto&resp1),data.peakVisualOptoFast_N(RS&P14P18&L4&Nkxopto&resp1))
[h,p,stats]=my_ttest(data.peakVisualFast_N(RS&P9P13&L56&Nkxopto&resp1),data.peakVisualOptoFast_N(RS&P9P13&L56&Nkxopto&resp1))
[h,p,stats]=my_ttest(data.peakVisualFast_N(RS&P14P18&L56&Nkxopto&resp1),data.peakVisualOptoFast_N(RS&P14P18&L56&Nkxopto&resp1))
[h,p,stats]=my_ttest(data.peakVisualFast_N(FS&Nkxopto&resp1),data.peakVisualOptoFast_N(FS&Nkxopto&resp1))


figure('units','normalized','outerposition',[0 0 0.2 1]);
ax(1)=subplot(4,2,1);
hold on
plot([-3,-2],[data.rw_meanPSTH(RS&P9P13&L23&Nkxopto&resp2),data.rwOpto_meanPSTH(RS&P9P13&L23&Nkxopto&resp2)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.rw_meanPSTH(RS&P9P13&L23&Nkxopto&resp2),'omitnan'),mean(data.rwOpto_meanPSTH(RS&P9P13&L23&Nkxopto&resp2),'omitnan')],[sem(data.rw_meanPSTH(RS&P9P13&L23&Nkxopto&resp2)),sem(data.rwOpto_meanPSTH(RS&P9P13&L23&Nkxopto&resp2))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.rw_meanPSTH(RS&P14P18&L23&Nkxopto&resp2),data.rwOpto_meanPSTH(RS&P14P18&L23&Nkxopto&resp2)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.rw_meanPSTH(RS&P14P18&L23&Nkxopto&resp2),'omitnan'),mean(data.rwOpto_meanPSTH(RS&P14P18&L23&Nkxopto&resp2),'omitnan')],[sem(data.rw_meanPSTH(RS&P14P18&L23&Nkxopto&resp2)),sem(data.rwOpto_meanPSTH(RS&P14P18&L23&Nkxopto&resp2))],'k-o','Linewidth',1,'CapSize',10)

ax(2)=subplot(4,2,2);
hold on
plot([-3,-2],[data.rw_meanPSTH(RS&P9P13&L4&Nkxopto&resp2),data.rwOpto_meanPSTH(RS&P9P13&L4&Nkxopto&resp2)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.rw_meanPSTH(RS&P9P13&L4&Nkxopto&resp2),'omitnan'),mean(data.rwOpto_meanPSTH(RS&P9P13&L4&Nkxopto&resp2),'omitnan')],[sem(data.rw_meanPSTH(RS&P9P13&L4&Nkxopto&resp2)),sem(data.rwOpto_meanPSTH(RS&P9P13&L4&Nkxopto&resp2))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.rw_meanPSTH(RS&P14P18&L4&Nkxopto&resp2),data.rwOpto_meanPSTH(RS&P14P18&L4&Nkxopto&resp2)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.rw_meanPSTH(RS&P14P18&L4&Nkxopto&resp2),'omitnan'),mean(data.rwOpto_meanPSTH(RS&P14P18&L4&Nkxopto&resp2),'omitnan')],[sem(data.rw_meanPSTH(RS&P14P18&L4&Nkxopto&resp2)),sem(data.rwOpto_meanPSTH(RS&P14P18&L4&Nkxopto&resp2))],'k-o','Linewidth',1,'CapSize',10)

ax(3)=subplot(4,2,3);
hold on
plot([-3,-2],[data.rw_meanPSTH(RS&P9P13&L56&Nkxopto&resp2),data.rwOpto_meanPSTH(RS&P9P13&L56&Nkxopto&resp2)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.rw_meanPSTH(RS&P9P13&L56&Nkxopto&resp2),'omitnan'),mean(data.rwOpto_meanPSTH(RS&P9P13&L56&Nkxopto&resp2),'omitnan')],[sem(data.rw_meanPSTH(RS&P9P13&L56&Nkxopto&resp2)),sem(data.rwOpto_meanPSTH(RS&P9P13&L56&Nkxopto&resp2))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.rw_meanPSTH(RS&P14P18&L56&Nkxopto&resp2),data.rwOpto_meanPSTH(RS&P14P18&L56&Nkxopto&resp2)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.rw_meanPSTH(RS&P14P18&L56&Nkxopto&resp2),'omitnan'),mean(data.rwOpto_meanPSTH(RS&P14P18&L56&Nkxopto&resp2),'omitnan')],[sem(data.rw_meanPSTH(RS&P14P18&L56&Nkxopto&resp2)),sem(data.rwOpto_meanPSTH(RS&P14P18&L56&Nkxopto&resp2))],'k-o','Linewidth',1,'CapSize',10)

ax(4)=subplot(4,2,4);
hold on
plot([-1,0],[data.rw_meanPSTH(FS&P14P18&Nkxopto&resp2),data.rwOpto_meanPSTH(FS&P14P18&Nkxopto&resp2)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.rw_meanPSTH(FS&P14P18&Nkxopto&resp2),'omitnan'),mean(data.rwOpto_meanPSTH(FS&P14P18&Nkxopto&resp2),'omitnan')],[sem(data.rw_meanPSTH(FS&P14P18&Nkxopto&resp2)),sem(data.rwOpto_meanPSTH(FS&P14P18&Nkxopto&resp2))],'k-o','Linewidth',1,'CapSize',10)

ylims=[0,15;0,15;0,15;0,25];
for i=1:numel(ax)
    ax(i).XLim=[-3.5,.5];
    ax(i).YAxis.Label.String='Spike/s';
    ax(i).FontSize=10;
    ax(i).YLim=ylims(i,:);
end

export_fig(fullfile(folderFigures,'4.14','ErrorBars_SlowRaw'),'-pdf','-transparent','-nocrop')
close

[h,p,stats]=my_ttest(data.rw_meanPSTH(RS&P9P13&L23&Nkxopto&resp2),data.rwOpto_meanPSTH(RS&P9P13&L23&Nkxopto&resp2))
[h,p,stats]=my_ttest(data.rw_meanPSTH(RS&P14P18&L23&Nkxopto&resp2),data.rwOpto_meanPSTH(RS&P14P18&L23&Nkxopto&resp2))
[h,p,stats]=my_ttest(data.rw_meanPSTH(RS&P9P13&L4&Nkxopto&resp2),data.rwOpto_meanPSTH(RS&P9P13&L4&Nkxopto&resp2))
[h,p,stats]=my_ttest(data.rw_meanPSTH(RS&P14P18&L4&Nkxopto&resp2),data.rwOpto_meanPSTH(RS&P14P18&L4&Nkxopto&resp2))
[h,p,stats]=my_ttest(data.rw_meanPSTH(RS&P9P13&L56&Nkxopto&resp2),data.rwOpto_meanPSTH(RS&P9P13&L56&Nkxopto&resp2))
[h,p,stats]=my_ttest(data.rw_meanPSTH(RS&P14P18&L56&Nkxopto&resp2),data.rwOpto_meanPSTH(RS&P14P18&L56&Nkxopto&resp2))
[h,p,stats]=my_ttest(data.rw_meanPSTH(FS&Nkxopto&resp2),data.rwOpto_meanPSTH(FS&Nkxopto&resp2))


figure('units','normalized','outerposition',[0 0 0.2 1]);
ax(1)=subplot(4,2,1);
hold on
plot([-3,-2],[data.fanoVisual(RS&P9P13&L23&Nkxopto&resp1),data.fanoVisualOpto(RS&P9P13&L23&Nkxopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.fanoVisual(RS&P9P13&L23&Nkxopto&resp1),'omitnan'),mean(data.fanoVisualOpto(RS&P9P13&L23&Nkxopto&resp1),'omitnan')],[sem(data.fanoVisual(RS&P9P13&L23&Nkxopto&resp1)),sem(data.fanoVisualOpto(RS&P9P13&L23&Nkxopto&resp1))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.fanoVisual(RS&P14P18&L23&Nkxopto&resp1),data.fanoVisualOpto(RS&P14P18&L23&Nkxopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.fanoVisual(RS&P14P18&L23&Nkxopto&resp1),'omitnan'),mean(data.fanoVisualOpto(RS&P14P18&L23&Nkxopto&resp1),'omitnan')],[sem(data.fanoVisual(RS&P14P18&L23&Nkxopto&resp1)),sem(data.fanoVisualOpto(RS&P14P18&L23&Nkxopto&resp1))],'k-o','Linewidth',1,'CapSize',10)

ax(2)=subplot(4,2,2);
hold on
plot([-3,-2],[data.fanoVisual(RS&P9P13&L4&Nkxopto&resp1),data.fanoVisualOpto(RS&P9P13&L4&Nkxopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.fanoVisual(RS&P9P13&L4&Nkxopto&resp1),'omitnan'),mean(data.fanoVisualOpto(RS&P9P13&L4&Nkxopto&resp1),'omitnan')],[sem(data.fanoVisual(RS&P9P13&L4&Nkxopto&resp1)),sem(data.fanoVisualOpto(RS&P9P13&L4&Nkxopto&resp1))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.fanoVisual(RS&P14P18&L4&Nkxopto&resp1),data.fanoVisualOpto(RS&P14P18&L4&Nkxopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.fanoVisual(RS&P14P18&L4&Nkxopto&resp1),'omitnan'),mean(data.fanoVisualOpto(RS&P14P18&L4&Nkxopto&resp1),'omitnan')],[sem(data.fanoVisual(RS&P14P18&L4&Nkxopto&resp1)),sem(data.fanoVisualOpto(RS&P14P18&L4&Nkxopto&resp1))],'k-o','Linewidth',1,'CapSize',10)

ax(3)=subplot(4,2,3);
hold on
plot([-3,-2],[data.fanoVisual(RS&P9P13&L56&Nkxopto&resp1),data.fanoVisualOpto(RS&P9P13&L56&Nkxopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-3,-2],[mean(data.fanoVisual(RS&P9P13&L56&Nkxopto&resp1),'omitnan'),mean(data.fanoVisualOpto(RS&P9P13&L56&Nkxopto&resp1),'omitnan')],[sem(data.fanoVisual(RS&P9P13&L56&Nkxopto&resp1)),sem(data.fanoVisualOpto(RS&P9P13&L56&Nkxopto&resp1))],'k-o','Linewidth',1,'CapSize',10)
plot([-1,0],[data.fanoVisual(RS&P14P18&L56&Nkxopto&resp1),data.fanoVisualOpto(RS&P14P18&L56&Nkxopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.fanoVisual(RS&P14P18&L56&Nkxopto&resp1),'omitnan'),mean(data.fanoVisualOpto(RS&P14P18&L56&Nkxopto&resp1),'omitnan')],[sem(data.fanoVisual(RS&P14P18&L56&Nkxopto&resp1)),sem(data.fanoVisualOpto(RS&P14P18&L56&Nkxopto&resp1))],'k-o','Linewidth',1,'CapSize',10)

ax(4)=subplot(4,2,4);
hold on
plot([-1,0],[data.fanoVisual(FS&P14P18&Nkxopto&resp1),data.fanoVisualOpto(FS&P14P18&Nkxopto&resp1)],'-','Color',[150,150,150]/255)
errorbar([-1,0],[mean(data.fanoVisual(FS&P14P18&Nkxopto&resp1),'omitnan'),mean(data.fanoVisualOpto(FS&P14P18&Nkxopto&resp1),'omitnan')],[sem(data.fanoVisual(FS&P14P18&Nkxopto&resp1)),sem(data.fanoVisualOpto(FS&P14P18&Nkxopto&resp1))],'k-o','Linewidth',1,'CapSize',10)

% ylims=[0,1.5;0,1.5;0,1.5;0,1.5];
for i=1:numel(ax)
    ax(i).XLim=[-3.5,.5];
    ax(i).YAxis.Label.String='Fano factor';
    ax(i).FontSize=10;
    ax(i).YLim=[0, 2.5];
end

export_fig(fullfile(folderFigures,'4.14','FanoFactor'),'-pdf','-transparent','-nocrop')
close

[h,p,stats]=my_ttest(data.fanoVisual(RS&P9P13&L23&Nkxopto&resp1),data.fanoVisualOpto(RS&P9P13&L23&Nkxopto&resp1))
[h,p,stats]=my_ttest(data.fanoVisual(RS&P14P18&L23&Nkxopto&resp1),data.fanoVisualOpto(RS&P14P18&L23&Nkxopto&resp1))
[h,p,stats]=my_ttest(data.fanoVisual(RS&P14P18&L4&Nkxopto&resp1),data.fanoVisualOpto(RS&P14P18&L4&Nkxopto&resp1))
[h,p,stats]=my_ttest(data.fanoVisual(RS&P9P13&L56&Nkxopto&resp1),data.fanoVisualOpto(RS&P9P13&L56&Nkxopto&resp1))
[h,p,stats]=my_ttest(data.fanoVisual(RS&P14P18&L56&Nkxopto&resp1),data.fanoVisualOpto(RS&P14P18&L56&Nkxopto&resp1))
[h,p,stats]=my_ttest(data.fanoVisual(FS&Nkxopto&resp1),data.fanoVisualOpto(FS&Nkxopto&resp1))

%% Chemogenetic visual
chemoYBar=[nnz(data.responseVisual(KORD&P9P13))/nnz(KORD&P9P13),nnz(data.responseVisual_K(KORD&P9P13))/nnz(KORD&P9P13);...
    nnz(data.rw_responsive(KORD&P9P13))/nnz(KORD&P9P13),nnz(data.rw_responsive_K(KORD&P9P13))/nnz(KORD&P9P13);...
    nnz(data.responseVisual(KORD&P14P18))/nnz(KORD&P14P18),nnz(data.responseVisual_K(KORD&P14P18))/nnz(KORD&P14P18);...
    nnz(data.rw_responsive(KORD&P14P18))/nnz(KORD&P14P18),nnz(data.rw_responsive(KORD&P14P18))/nnz(KORD&P14P18);];
chemoYBar2=[nnz(SB_entrained(KORD&P9P13))/nnz(KORD&P9P13),nnz(SB_entrained_K(KORD&P9P13))/nnz(KORD&P9P13);...
            nnz(PPC_entrained(KORD&P9P13))/nnz(KORD&P9P13),nnz(PPC_entrained_K(KORD&P9P13))/nnz(KORD&P9P13);...
            nnz(SB_entrained(KORD&P14P18))/nnz(KORD&P14P18),nnz(SB_entrained_K(KORD&P14P18))/nnz(KORD&P14P18);...
            nnz(PPC_entrained(KORD&P14P18))/nnz(KORD&P14P18),nnz(PPC_entrained_K(KORD&P14P18))/nnz(KORD&P14P18);];

figure('units','normalized','outerposition',[0 0 0.15 1]);
ax=subplot(3,1,1);
b = bar([1,2,4,5],chemoYBar);
for i=1:size(b,2)
    b(1,i).LineWidth=1;
end
ax.YLabel.String='Ratio single units';
ax.Box='off';
ax.LineWidth = 1;
ax.FontSize=10;
ax.YLim=[0,1];
legend('Control','KORD-SalB')
legend('boxoff')

ax=subplot(3,1,2);
b = bar([1,2,4,5],chemoYBar2);
for i=1:size(b,2)
    b(1,i).LineWidth=1;
end
ax.YLabel.String='Ratio single units';
ax.Box='off';
ax.LineWidth = 1;
ax.FontSize=10;
ax.YLim=[0,1];
legend('Control','KORD-SalB')
legend('boxoff')
export_fig(fullfile(folderFigures,'4.18','Barplots'),'-pdf','-transparent','-nocrop')
close

figure('units','normalized','outerposition',[0 0 0.15 1]);
ax(1)=subplot(3,1,1);
hold on
plot(ones(height(data(KORD&P9P13&resp1,:)),1),data.peakVisualFast_Change(KORD&P9P13&resp1),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.peakVisualFast_Change(KORD&P9P13&resp1),'omitnan'),sem(data.peakVisualFast_Change(KORD&P9P13&resp1)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18&resp1,:)),1)*4,data.peakVisualFast_Change(KORD&P14P18&resp1),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.peakVisualFast_Change(KORD&P14P18&resp1),'omitnan'),sem(data.peakVisualFast_Change(KORD&P14P18&resp1)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P9P13&resp1,:)),1)*2,data.peakVisualFast_Change(GFP&P9P13&resp1),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.peakVisualFast_Change(GFP&P9P13&resp1),'omitnan'),sem(data.peakVisualFast_Change(GFP&P9P13&resp1)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18&resp1,:)),1)*5,data.peakVisualFast_Change(GFP&P14P18&resp1),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.peakVisualFast_Change(GFP&P14P18&resp1),'omitnan'),sem(data.peakVisualFast_Change(GFP&P14P18&resp1)),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(1).YAxis.Label.String='Max PSTH fast change (ratio log2)';

[h,p, stats] = my_ttest2(data.peakVisualFast_Change(KORD&P14P18&resp1),data.peakVisualFast_Change(GFP&P14P18&resp1))

ax(2)=subplot(3,1,2);
hold on
plot(ones(height(data(KORD&P9P13&resp2,:)),1),data.peakVisualSlow_Change(KORD&P9P13&resp2),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.peakVisualSlow_Change(KORD&P9P13&resp2),'omitnan'),sem(data.peakVisualSlow_Change(KORD&P9P13&resp2)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18&resp2,:)),1)*4,data.peakVisualSlow_Change(KORD&P14P18&resp2),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.peakVisualSlow_Change(KORD&P14P18&resp2),'omitnan'),sem(data.peakVisualSlow_Change(KORD&P14P18&resp2)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P9P13&resp2,:)),1)*2,data.peakVisualSlow_Change(GFP&P9P13&resp2),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.peakVisualSlow_Change(GFP&P9P13&resp2),'omitnan'),sem(data.peakVisualSlow_Change(GFP&P9P13&resp2)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18&resp2,:)),1)*5,data.peakVisualSlow_Change(GFP&P14P18&resp2),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.peakVisualSlow_Change(GFP&P14P18&resp2),'omitnan'),sem(data.peakVisualSlow_Change(GFP&P14P18&resp2)),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(2).YAxis.Label.String='Max PSTH slow change (ratio log2)';

[h,p, stats] = my_ttest2(data.peakVisualSlow_Change(KORD&P9P13&resp2),data.peakVisualSlow_Change(GFP&P9P13&resp2))
[h,p,~, stats] = ttest2(data.peakVisualSlow_Change(KORD&P14P18&resp2),data.peakVisualSlow_Change(GFP&P14P18&resp2))


ax(3)=subplot(3,1,3);
hold on
plot(ones(height(data(KORD&P9P13,:)),1),data.baseline_Change(KORD&P9P13),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.baseline_Change(KORD&P9P13),'omitnan'),sem(data.baseline_Change(KORD&P9P13)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*4,data.baseline_Change(KORD&P14P18),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.baseline_Change(KORD&P14P18),'omitnan'),sem(data.baseline_Change(KORD&P14P18)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P9P13,:)),1)*2,data.baseline_Change(GFP&P9P13),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.baseline_Change(GFP&P9P13),'omitnan'),sem(data.baseline_Change(GFP&P9P13)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18,:)),1)*5,data.baseline_Change(GFP&P14P18),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.baseline_Change(GFP&P14P18),'omitnan'),sem(data.baseline_Change(GFP&P14P18)),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(3).YAxis.Label.String='Baseline change (ratio log2)';

[h,p,stats] = my_ttest2(data.baseline_Change(KORD&P9P13),data.baseline_Change(GFP&P9P13))
[h,p,stats] = my_ttest2(data.baseline_Change(KORD&P14P18),data.baseline_Change(GFP&P14P18))

for i=1:numel(ax)
    ax(i).XLim=[0.5,5.5];
    ax(i).FontSize=10;
    ax(i).YLim=[-2,2];
end
export_fig(fullfile(folderFigures,'4.18','Errorbar responsivity'),'-pdf','-transparent','-nocrop')
close

figure('units','normalized','outerposition',[0 0 0.15 1]);
ax(1)=subplot(3,1,1);
hold on
plot(ones(height(data(KORD&P9P13&resp1,:)),1),data.fano_Change(KORD&P9P13&resp1),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.fano_Change(KORD&P9P13&resp1),'omitnan'),sem(data.fano_Change(KORD&P9P13&resp1)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18&resp1,:)),1)*4,data.fano_Change(KORD&P14P18&resp1),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.fano_Change(KORD&P14P18&resp1),'omitnan'),sem(data.fano_Change(KORD&P14P18&resp1)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P9P13&resp1,:)),1)*2,data.fano_Change(GFP&P9P13&resp1),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.fano_Change(GFP&P9P13&resp1),'omitnan'),sem(data.fano_Change(GFP&P9P13&resp1)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18&resp1,:)),1)*5,data.fano_Change(GFP&P14P18&resp1),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.fano_Change(GFP&P14P18&resp1),'omitnan'),sem(data.fano_Change(GFP&P14P18&resp1)),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(1).YAxis.Label.String='Fano factor change (ratio log2)';

[h,p,stats] = my_ttest2(data.fano_Change(KORD&P14P18&resp1),data.fano_Change(GFP&P14P18&resp1))

ax(2)=subplot(3,1,2);
hold on
plot(ones(height(data(KORD&P9P13&SB_entrained,:)),1),data.sbSpikeProb_Change(KORD&P9P13&SB_entrained),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.sbSpikeProb_Change(KORD&P9P13&SB_entrained),'omitnan'),sem(data.sbSpikeProb_Change(KORD&P9P13&SB_entrained)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18&SB_entrained,:)),1)*4,data.sbSpikeProb_Change(KORD&P14P18&SB_entrained),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.sbSpikeProb_Change(KORD&P14P18&SB_entrained),'omitnan'),sem(data.sbSpikeProb_Change(KORD&P14P18&SB_entrained)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P9P13&SB_entrained,:)),1)*2,data.sbSpikeProb_Change(GFP&P9P13&SB_entrained),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.sbSpikeProb_Change(GFP&P9P13&SB_entrained),'omitnan'),sem(data.sbSpikeProb_Change(GFP&P9P13&SB_entrained)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18&SB_entrained,:)),1)*5,data.sbSpikeProb_Change(GFP&P14P18&SB_entrained),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.sbSpikeProb_Change(GFP&P14P18&SB_entrained),'omitnan'),sem(data.sbSpikeProb_Change(GFP&P14P18&SB_entrained)),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(2).YAxis.Label.String='SB probability change (ratio log2)';

[h,p,stats] = my_ttest2(data.sbSpikeProb_Change(KORD&P9P13&SB_entrained),data.sbSpikeProb_Change(GFP&P9P13&SB_entrained))

ax(3)=subplot(3,1,3);
hold on
plot(ones(height(data(KORD&P9P13&PPC_entrained,:)),1),data.PPC_Change(KORD&P9P13&PPC_entrained),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.PPC_Change(KORD&P9P13&PPC_entrained),'omitnan'),sem(data.PPC_Change(KORD&P9P13&PPC_entrained)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18&PPC_entrained,:)),1)*4,data.PPC_Change(KORD&P14P18&PPC_entrained),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.PPC_Change(KORD&P14P18&PPC_entrained),'omitnan'),sem(data.PPC_Change(KORD&P14P18&PPC_entrained)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P9P13&PPC_entrained,:)),1)*2,data.PPC_Change(GFP&P9P13&PPC_entrained),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.PPC_Change(GFP&P9P13&PPC_entrained),'omitnan'),sem(data.PPC_Change(GFP&P9P13&PPC_entrained)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18&PPC_entrained,:)),1)*5,data.PPC_Change(GFP&P14P18&PPC_entrained),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.PPC_Change(GFP&P14P18&PPC_entrained),'omitnan'),sem(data.PPC_Change(GFP&P14P18&PPC_entrained)),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(3).YAxis.Label.String='PPC beta change (ratio log2)';

[h,p,stats] = my_ttest2(data.PPC_Change(KORD&P9P13&PPC_entrained),data.PPC_Change(GFP&P9P13&PPC_entrained))

for i=1:numel(ax)
    ax(i).XLim=[0.5,5.5];
    ax(i).FontSize=10;
    ax(i).YLim=[-2,2];
end
export_fig(fullfile(folderFigures,'4.18','Errorbar spontaneous'),'-pdf','-transparent','-nocrop')
close

figure('units','normalized','outerposition',[0 0 0.15 1]);
ax(1)=subplot(3,1,1);
hold on
plot(ones(height(data(KORD&P14P18&resp1&RS,:)),1),data.peakVisualFast_Change(KORD&P14P18&resp1&RS),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.peakVisualFast_Change(KORD&P14P18&resp1&RS),'omitnan'),sem(data.peakVisualFast_Change(KORD&P14P18&resp1&RS)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18&resp1&FS,:)),1)*4,data.peakVisualFast_Change(KORD&P14P18&resp1&FS),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.peakVisualFast_Change(KORD&P14P18&resp1&FS),'omitnan'),sem(data.peakVisualFast_Change(KORD&P14P18&resp1&FS)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18&resp1&RS,:)),1)*2,data.peakVisualFast_Change(GFP&P14P18&resp1&RS),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.peakVisualFast_Change(GFP&P14P18&resp1&RS),'omitnan'),sem(data.peakVisualFast_Change(GFP&P14P18&resp1&RS)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18&resp1&FS,:)),1)*5,data.peakVisualFast_Change(GFP&P14P18&resp1&FS),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.peakVisualFast_Change(GFP&P14P18&resp1&FS),'omitnan'),sem(data.peakVisualFast_Change(GFP&P14P18&resp1&FS)),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(1).YAxis.Label.String='Max PSTH fast change (ratio log2)';

[h,p,stats] = my_ttest2(data.peakVisualFast_Change(KORD&P14P18&resp1&RS),data.peakVisualFast_Change(GFP&P14P18&resp1&RS))
[h,p,stats] = my_ttest2(data.peakVisualFast_Change(KORD&P14P18&resp1&FS),data.peakVisualFast_Change(GFP&P14P18&resp1&FS))

ax(2)=subplot(3,1,2);
hold on
plot(ones(height(data(KORD&P14P18&resp1&RS,:)),1),data.fano_Change(KORD&P14P18&resp1&RS),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.fano_Change(KORD&P14P18&resp1&RS),'omitnan'),sem(data.fano_Change(KORD&P14P18&resp1&RS)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18&resp1&FS,:)),1)*4,data.fano_Change(KORD&P14P18&resp1&FS),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.fano_Change(KORD&P14P18&resp1&FS),'omitnan'),sem(data.fano_Change(KORD&P14P18&resp1&FS)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18&resp1&RS,:)),1)*2,data.fano_Change(GFP&P14P18&resp1&RS),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.fano_Change(GFP&P14P18&resp1&RS),'omitnan'),sem(data.fano_Change(GFP&P14P18&resp1&RS)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18&resp1&FS,:)),1)*5,data.fano_Change(GFP&P14P18&resp1&FS),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.fano_Change(GFP&P14P18&resp1&FS),'omitnan'),sem(data.fano_Change(GFP&P14P18&resp1&FS)),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(2).YAxis.Label.String='Fano factor change (ratio log2)';

[h,p] = my_ttest2(data.fano_Change(KORD&P14P18&resp1&RS),data.fano_Change(GFP&P14P18&resp1&RS))
[h,p] = my_ttest2(data.fano_Change(KORD&P14P18&resp1&FS),data.fano_Change(GFP&P14P18&resp1&FS))

ax(3)=subplot(3,1,3);
hold on
plot(ones(height(data(KORD&P14P18&RS,:)),1),data.baseline_Change(KORD&P14P18&RS),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.baseline_Change(KORD&P14P18&RS),'omitnan'),sem(data.baseline_Change(KORD&P14P18&RS)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18&FS,:)),1)*4,data.baseline_Change(KORD&P14P18&FS),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.baseline_Change(KORD&P14P18&FS),'omitnan'),sem(data.baseline_Change(KORD&P14P18&FS)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18&RS,:)),1)*2,data.baseline_Change(GFP&P14P18&RS),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.baseline_Change(GFP&P14P18&RS),'omitnan'),sem(data.baseline_Change(GFP&P14P18&RS)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18&FS,:)),1)*5,data.baseline_Change(GFP&P14P18&FS),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.baseline_Change(GFP&P14P18&FS),'omitnan'),sem(data.baseline_Change(GFP&P14P18&FS)),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(3).YAxis.Label.String='Baseline change (ratio log2)';

[h,p] = my_ttest2(data.baseline_Change(KORD&P14P18&RS),data.baseline_Change(GFP&P14P18&RS))
[h,p] = my_ttest2(data.baseline_Change(KORD&P14P18&FS),data.baseline_Change(GFP&P14P18&FS))

for i=1:numel(ax)
    ax(i).XLim=[0.5,5.5];
    ax(i).FontSize=10;
    ax(i).YLim=[-2,2];
end
export_fig(fullfile(folderFigures,'4.18','Errorbar RSvsFS'),'-pdf','-transparent','-nocrop')
close
%% Bar plot spindle bursts
% SST
figure('units','normalized','outerposition',[0 0 0.2 1]);
clear ax
ax(1)=subplot(3,1,1);
yHistSB=[nnz(SB_entrained&SST&P9P13)/nnz(SST&P9P13),nnz(SB_entrained&Nkx&P9P13)/nnz(Nkx&P9P13);nnz(SB_entrained&SST&P14P18)/nnz(SST&P14P18),nnz(SB_entrained&Nkx&P14P18)/nnz(Nkx&P14P18)];
xHistDev=reordercats(categorical(cellstr({'P9-P13','P14-P18'})),{'P9-P13','P14-P18'});
b=bar(xHistDev,yHistSB);
for j=1:size(b,2)
    b(1,j).LineWidth=1;
end

ax(2)=subplot(3,1,2);
yHistSB=[nnz(PPC_entrained&SST&P9P13)/nnz(SST&P9P13),nnz(PPC_entrained&Nkx&P9P13)/nnz(Nkx&P9P13);nnz(PPC_entrained&SST&P14P18)/nnz(SST&P14P18),nnz(PPC_entrained&Nkx&P14P18)/nnz(Nkx&P14P18)];
xHistDev=reordercats(categorical(cellstr({'P9-P13','P14-P18'})),{'P9-P13','P14-P18'});
b=bar(xHistDev,yHistSB);
for j=1:size(b,2)
    b(1,j).LineWidth=1;
end

for i=1:numel(ax)
    ax(i).YLabel.String='Ratio responsive single units';
    ax(i).Box='off';
    ax(i).LineWidth = 1.5;
    ax(i).FontSize=15;
    legend('SST','Nkx','Location','northeast')
    legend('boxoff')
end

% ax.XLabel.String='Postnatal day';

% text(0.76,yHistSB(1,1)+0.05,strcat(int2str(nnz(SST&resp1&P9P13)),'/',int2str(nnz(SST&P9P13))),'FontSize',20)
% text(1.03,yHistSB(1,2)+0.05,strcat(int2str(nnz(SST&resp2&P9P13)),'/',int2str(nnz(SST&P9P13))),'FontSize',20)
% text(1.73,yHistSB(2,1)+0.05,strcat(int2str(nnz(SST&resp1&P14P18)),'/',int2str(nnz(SST&P14P18))),'FontSize',20)
% text(1.02,yHistSB(2,2)+0.05,strcat(int2str(nnz(SST&resp2&P14P18)),'/',int2str(nnz(SST&P14P18))),'FontSize',20)

export_fig(fullfile(folderFigures,'4.8','RatioSBentrainedUnits'),'-pdf','-transparent','-nocrop')
close

%%
data.cellIdentity(Nkx)=categorical(cellstr('Nkx2-1'));
clear v ax
figure('units','normalized','outerposition',[0 0 0.2 0.7])
ax(1)=subplot(3,2,1);
v(1,:)=violinplot(data.sb_spikeProb(SB_entrained&P9P13&(SST|Nkx)), removecats(data.cellIdentity(SB_entrained&P9P13&(SST|Nkx))));
ax(1).YLim=[0, 1];
% [h,p,ci,stats] = ttest2(data.sb_spikeProb(SST&SB_entrained&P9P13),data.sb_spikeProb(Nkx&SB_entrained&P9P13));

ax(2)=subplot(3,2,2);
v(2,:)=violinplot(data.sb_spikeProb(SB_entrained&P14P18&(SST|Nkx)), removecats(data.cellIdentity(SB_entrained&P14P18&(SST|Nkx))));
ax(2).YLim=[0, 1];

ax(3)=subplot(3,2,3);
v(3,:)=violinplot(data.PPC((SST|Nkx)&(P9P13|P5P8)&PPC_entrained,2), removecats(data.cellIdentity((SST|Nkx)&(P9P13|P5P8)&PPC_entrained)));
ax(3).YLim=[0, 1];

ax(3).YScale='log';
ax(3).YTick=[0,0.01,0.1,1];

% [h,p,ci,stats] = ttest2(data.PPC(SST&SB_entrained&(P9P13|P5P8),2),data.PPC(Nkx&SB_entrained&(P9P13|P5P8),2));

ax(4)=subplot(3,2,4);
v(4,:)=violinplot(data.PPC((SST|Nkx)&P14P18&PPC_entrained,2), removecats(data.cellIdentity((SST|Nkx)&P14P18&PPC_entrained)));
ax(4).YLim=[0, .8];


ax(5)=subplot(3,2,5);
v(5,:)=violinplot(wrapTo360(rad2deg(data.vectorAngle((SST|Nkx)&(P9P13|P5P8)&PPC_entrained,2))), removecats(data.cellIdentity((SST|Nkx)&(P9P13|P5P8)&PPC_entrained)));
ax(5).YLim=[0, 360];
% [h,p,ci,stats] = ttest2(data.PPC(SST&SB_entrained&P9P13),data.PPC(Nkx&SB_entrained&P9P13));
ax(5).YGrid = 'on';
ax(5).YTick=[0:90:360];
[h,p,ci,stats] = ttest2(wrapTo360(rad2deg(data.vectorAngle(SST&SB_entrained&(P9P13|P5P8),2))),wrapTo360(rad2deg(data.vectorAngle(Nkx&SB_entrained&(P9P13|P5P8),2))));

ax(6)=subplot(3,2,6);
v(6,:)=violinplot(wrapTo360(rad2deg(data.vectorAngle((SST|Nkx)&P14P18&PPC_entrained,2))), removecats(data.cellIdentity((SST|Nkx)&P14P18&PPC_entrained)));
ax(6).YLim=[0, 360];
ax(6).YGrid = 'on';
ax(6).YTick=[0:90:360];

for i=1:numel(ax)
    for j=1:size(v,2)
        v(i,j).ViolinColor=[100,100,100]/255;
        v(i,j).ScatterPlot.MarkerFaceColor=[37,37,37]/255;
        v(i,j).ScatterPlot.MarkerFaceAlpha=1;
        v(i,j).ScatterPlot.SizeData=2;
    end

    ax(i).FontSize=12;
    ax(i).LineWidth=1;
    ax(i).XLim=[0.5,3.5];
    
end



print(gcf,'-dpdf','C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\WT Doctoral Programme in Neuroscience\DPhil in Neuroscience\Manuscripts\V1 S1\Figures\Fig. 2\S1BF_PPC_SB_violinIN')
close


clear ax
figure('units','normalized','outerposition',[0 0 0.2 1])
subplot(3,2,1)
polarscatter(wrapTo360(rad2deg(data.vectorAngle(SST&P9P13&PPC_entrained,2))),data.vectorLength(SST&P9P13&PPC_entrained,2),'b.')
subplot(3,2,3)
polarscatter(wrapTo360(rad2deg(data.vectorAngle(Nkx&P9P13&PPC_entrained,2))),data.vectorLength(Nkx&P9P13&PPC_entrained,2),'r.')

subplot(3,2,2)
polarscatter(wrapTo360(rad2deg(data.vectorAngle(SST&P14P18&PPC_entrained,2))),data.vectorLength(SST&P14P18&PPC_entrained,2),'b.')
subplot(3,2,4)
polarscatter(wrapTo360(rad2deg(data.vectorAngle(Nkx&P14P18&PPC_entrained,2))),data.vectorLength(Nkx&P14P18&PPC_entrained,2),'r.')

print(gcf,'-dpdf',fullfile(folderFigures,'4.6','polarscatterplots_INs'))
close


%% Imagesc plots - All

%SST
figure('units','normalized','outerposition',[0 0 1 0.5]);
ax=subplot(1,2,1);
imagesc(PSTHbins,1:height(data(SST&P9P13,:)),sortrows(data(SST&P9P13,:),{'Layer','MaxPSTH_Idx'}).PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(SST&P9P13,:).Layer=='L2/3'),nnz(data(SST&P9P13,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(SST&P9P13,:).Layer=='L5/6')),nnz(not(data(SST&P9P13,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-100,1100];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P9-P13';

ax=subplot(1,2,2);
imagesc(PSTHbins,1:height(data(SST&P14P18,:)),sortrows(data(SST&P14P18,:),{'Layer','MaxPSTH_Idx'}).PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(SST&P14P18,:).Layer=='L2/3'),nnz(data(SST&P14P18,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(SST&P14P18,:).Layer=='L5/6')),nnz(not(data(SST&P14P18,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-100,1100];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P14-P18';

sgtitle('SST','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH_img_SST_all'),'-pdf','-tiff','-transparent','-nocrop')
close
    
%Nkx
figure('units','normalized','outerposition',[0 0 1 0.5]);
ax=subplot(1,2,1);
imagesc(PSTHbins,1:height(data(Nkx&P9P13,:)),sortrows(data(Nkx&P9P13,:),{'Layer','MaxPSTH_Idx'}).PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(Nkx&P9P13,:).Layer=='L2/3'),nnz(data(Nkx&P9P13,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(Nkx&P9P13,:).Layer=='L5/6')),nnz(not(data(Nkx&P9P13,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-100,1100];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P9-P13';

ax=subplot(1,2,2);
imagesc(PSTHbins,1:height(data(Nkx&P14P18,:)),sortrows(data(Nkx&P14P18,:),{'Layer','MaxPSTH_Idx'}).PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(Nkx&P14P18,:).Layer=='L2/3'),nnz(data(Nkx&P14P18,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(Nkx&P14P18,:).Layer=='L5/6')),nnz(not(data(Nkx&P14P18,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-100,1100];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P14-P18';

sgtitle('Nkx2-1','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH_img_Nkx2-1_all'),'-pdf','-tiff','-transparent','-nocrop')
close

%RS
figure('units','normalized','outerposition',[0 0 1 1]);
ax=subplot(1,2,1);
imagesc(PSTHbins,1:height(data(RS&P9P13,:)),sortrows(data(RS&P9P13,:),{'Layer','MaxPSTH_Idx'}).PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(RS&P9P13,:).Layer=='L2/3'),nnz(data(RS&P9P13,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(RS&P9P13,:).Layer=='L5/6')),nnz(not(data(RS&P9P13,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-100,1100];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P9-P13';

ax=subplot(1,2,2);
imagesc(PSTHbins,1:height(data(RS&P14P18,:)),sortrows(data(RS&P14P18,:),{'Layer','MaxPSTH_Idx'}).PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(RS&P14P18,:).Layer=='L2/3'),nnz(data(RS&P14P18,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(RS&P14P18,:).Layer=='L5/6')),nnz(not(data(RS&P14P18,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-100,1100];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P14-P18';

sgtitle('RS','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH_img_RS_all'),'-pdf','-tiff','-transparent','-nocrop')
close

%FS
figure('units','normalized','outerposition',[0 0 1 1]);
ax=subplot(1,2,2);
imagesc(PSTHbins,1:height(data(FS&P14P18,:)),sortrows(data(FS&P14P18,:),{'Layer','MaxPSTH_Idx'}).PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(FS&P14P18,:).Layer=='L2/3'),nnz(data(FS&P14P18,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(FS&P14P18,:).Layer=='L5/6')),nnz(not(data(FS&P14P18,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-100,1100];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P14-P18';

sgtitle('FS','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH_img_FS_all'),'-pdf','-tiff','-transparent','-nocrop')
close

%% PSTH clustering
% Set logical groups for responsive cells
RS_resp=RS(resp);
FS_resp=FS(resp);
SST_resp=SST(resp);
Nkx_resp=Nkx(resp);
% Nkx_FS_resp=Nkx_FS(resp);
% Nkx_RS_resp=Nkx_RS(resp);

P14P18_resp=P14P18(resp);
P9P13_resp=P9P13(resp);

L23=(data.Layer=='L2/3');
L4=(data.Layer=='L4');
L56=(data.Layer=='L5/6');
L23_resp=L23(resp);
L4_resp=L4(resp);
L56_resp=L56(resp);
%%
close all
smoothingval=1;
x1=smooth(data.PSTHvisual(data.MouseID=='SC10' & data.suid==145,:),smoothingval)';
x2=smooth(data.PSTHvisual(data.MouseID=='SC2' & data.suid==161,:),smoothingval)';
x3=smooth(data.PSTHvisual(data.MouseID=='SC23' & data.suid==45,:),smoothingval)';
figure; hold on; plot(x1); plot(x2); plot(x3);
distance='cosine';
distancevals=[pdist2(x1,x2,distance),pdist2(x1,x3,distance),pdist2(x2,x3,distance)]
%%
x=data.PSTHvisual(resp,:);
% x=zscore(x,[],'all');
% x=x-mean(x(:,PSTHbins>-1000 & PSTHbins<-500),2); 

tsneOptions=struct;
tsneOptions.MaxIter=10000;
tsneOptions.OutputFcn=[];
tsneOptions.TolFun=1e-10;
x_score=tsne(x,'Algorithm','exact','Distance','correlation','NumDimensions',2,'Perplexity',20,'Options',tsneOptions);


colours={[255,227,145]/255,[0,0,179]/255,[203,24,29]/255,[107,174,214]/255,[253,141,60]/255};
figure('units','normalized','outerposition',[0 0 1 1]);
ax=subplot(2,3,2);
scatter(x_score(:,1),x_score(:,2),'filled','k')
ax.XLabel.String='t-SNE 1';
ax.YLabel.String='t-SNE 2';
ax.FontSize=18;
ax.LineWidth=2;
ax.XLim=[-50, 50];
ax.YLim=[-50, 50];
ax.Title.String=strcat('Perlexity=',int2str(perplexityVals(i)),' - Distance=',distance{i});
grid


%True label: cell types
ax=subplot(2,3,4);
hold on
scatter(x_score(RS_resp,1),x_score(RS_resp,2),'filled','MarkerFaceColor',colours{1})
scatter(x_score(FS_resp,1),x_score(FS_resp,2),'filled','MarkerFaceColor',colours{3})
scatter(x_score(SST_resp,1),x_score(SST_resp,2),'filled','MarkerFaceColor',colours{2})
scatter(x_score(Nkx_resp,1),x_score(Nkx_resp,2),'filled','MarkerFaceColor',colours{4})
% scatter(x_score(Nkx_FS_resp,1),x_score(Nkx_FS_resp,2),'filled','k')
% scatter(x_score(Nkx_RS_resp,1),x_score(Nkx_RS_resp,2),'filled','k')
ax.Title.String='Cell identity';
ax.XLabel.String='t-SNE 1';
ax.YLabel.String='t-SNE 2';
ax.FontSize=18;
ax.LineWidth=2;
% ax.XLim=[-140, 140];
% ax.YLim=[-99, 99];
grid
lg=legend('RS','FS','SST','Nkx2-1','Location','northwest');
lg.Box='off';

%True label: age
ax=subplot(2,3,5);
hold on
scatter(x_score(P9P13_resp,1),x_score(P9P13_resp,2),'filled','MarkerFaceColor',colours{2})
scatter(x_score(P14P18_resp,1),x_score(P14P18_resp,2),'filled','MarkerFaceColor',colours{3})
ax.Title.String='Development';
ax.XLabel.String='t-SNE 1';
ax.YLabel.String='t-SNE 2';
ax.FontSize=18;
ax.LineWidth=2;
% ax.XLim=[-140, 140];
% ax.YLim=[-99, 99];
grid
lg=legend('P9-P13','P14-P18','Location','northwest');
lg.Box='off';

%True label: layers
ax=subplot(2,3,6);
hold on
scatter(x_score(L56_resp,1),x_score(L56_resp,2),'filled','MarkerFaceColor',colours{1})
scatter(x_score(L23_resp,1),x_score(L23_resp,2),'filled','MarkerFaceColor',colours{2})
scatter(x_score(L4_resp,1),x_score(L4_resp,2),'filled','MarkerFaceColor',colours{3})
ax.XLabel.String='t-SNE 1';
ax.YLabel.String='t-SNE 2';
ax.Title.String='Cell layer';
ax.FontSize=18;
ax.LineWidth=2;
% ax.XLim=[-140, 140];
% ax.YLim=[-99, 99];
grid
lg=legend('L2/3','L4','L5/6','Location','northwest');
lg.Box='off';

% export_fig(fullfile(folderFigures,'PSTH-tSNE'),'-tiff','-transparent','-nocrop')
% close
%%
figure('units','normalized','outerposition',[0 0 0.5 0.5]);
Klabel = IterateKMeansClustering(x, 2:15, 'NumIterations', 1000, 'Distance','correlation','Verbose',true);

figure('units','normalized','outerposition',[0 0 1 1]);
ax=subplot(2,3,2);
hold on
for i=1:numel(unique(Klabel))
    sc=scatter(x_score(Klabel==i,1),x_score(Klabel==i,2),'filled');
    sc.Marker='o';
    sc.MarkerFaceColor=colours{i+1};
end
ax.XLabel.String='t-SNE 1';
ax.YLabel.String='t-SNE 2';
ax.Title.String='K-Means clustering ';
ax.FontSize=18;
ax.LineWidth=2;
ax.XLim=[-140, 140];
ax.YLim=[-99, 99];
grid
% lg=legend('C1','C2','Location','northwest');
% lg.Box='off';

smoothingVal=5;

ax1=subplot(2,3,4:6);
patch([0 100 100 0],[-0.8 -0.8 -0.3 -0.3],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
for i=1:numel(unique(Klabel))
    plot(PSTHbins,smooth(mean(x(Klabel==i,:),1),smoothingVal),'Color',colours{i+1},'LineWidth',2)
end
ax1=applyFont(ax1,0);
ax1.Color=[237,237,237]/255;
% ax1.YLim=[-1 4];
ax1.XLim=[-50 4000];
ax1.Title.FontSize=23;
% legend('C1','C2')
% export_fig(fullfile(folderFigures,'PSTH-tSNE_kMeans'),'-tiff','-transparent','-nocrop')
% close

%% Coherence
data.cellIdentity = reordercats(data.cellIdentity,{'RS','SST','FS','Nkx2-1','Nkx2-1 - RS','Nkx2-1 - FS'});
coherenceBands={'Low','Beta','Gamma'};
developmentStage={'P5-P8','P9-P13','P14-P18'};
agesLogical=[P5P8,P9P13,P14P18];

%All
figure('units','normalized','outerposition',[0 0 1 1]);
c=0;
for k=1:3
for j=1:3
    c=c+1;
    ax1=subplot(3,3,c);
    v=violinplot(data.coherence(agesLogical(:,k),j), removecats(data.cellIdentity(agesLogical(:,k))));
    for i=1:numel(v)
        v(i).ViolinColor=colours{i};
        v(i).ScatterPlot.MarkerFaceColor=[189,189,189]/255;
        v(i).ScatterPlot.MarkerFaceAlpha=1;
        v(i).ScatterPlot.SizeData=20;
    end
    ax1.FontSize=18;
    ax1.LineWidth=2;
    ax1.YLim=[0,.8];
    ax1.YLabel.String='Coherence';
    ax1.Title.String=strcat(developmentStage{k},' - ',coherenceBands{j});
    ax1.Title.FontSize=23;
    ax1.XTickLabelRotation=20;
    grid
end
end
export_fig(fullfile(folderFigures,'Coherence_All'),'-tiff','-transparent','-nocrop')
close

%Only responsive
figure('units','normalized','outerposition',[0 0 1 1]);
c=0;
for k=2:3
for j=1:3
    c=c+1;
    ax1=subplot(2,3,c);
    v=violinplot(data.coherence(resp&agesLogical(:,k),j), removecats(data.cellIdentity(resp&agesLogical(:,k))));
    for i=1:numel(v)
        v(i).ViolinColor=colours{i};
        v(i).ScatterPlot.MarkerFaceColor=[189,189,189]/255;
        v(i).ScatterPlot.MarkerFaceAlpha=1;
        v(i).ScatterPlot.SizeData=20;
    end
    ax1.FontSize=18;
    ax1.LineWidth=2;
    ax1.YLim=[0,.8];
    ax1.YLabel.String='Coherence';
    ax1.Title.String=strcat(developmentStage{k},' - ',coherenceBands{j});
    ax1.Title.FontSize=23;
    ax1.XTickLabelRotation=20;
    grid
end
end
export_fig(fullfile(folderFigures,'Coherence_responsive'),'-tiff','-transparent','-nocrop')
close

%% Only SST
figure('units','normalized','outerposition',[0 0 1 1]);
ax=subplot(2,3,1);
hold on
plot(data.PPC(BetaC&RS&P9P13,:)','-o','Color',[.5,.5,.5])
errorbar([1,2,3]+0.1,nanmean(data.PPC(BetaC&RS&P9P13,:)),nanstd(data.PPC(BetaC&RS&P9P13,:)),'k','LineWidth',2)
ax.XLim=[.5 3.5];
ax.YLabel.String='PPC';
ax.XTick=[1,2,3];
ax.XTickLabel=coherenceBands;
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P9-P13 - RS';
ax.YLim=[0 .8];

ax=subplot(2,3,2);
hold on
plot(data.PPC(BetaC&SST&P9P13,:)','-o','Color',[.5,.5,.5])
errorbar([1,2,3]+0.1,nanmean(data.PPC(BetaC&SST&P9P13,:)),nanstd(data.PPC(BetaC&SST&P9P13,:)),'k','LineWidth',2)
ax.XLim=[.5 3.5];
ax.YLabel.String='PPC';
ax.XTick=[1,2,3];
ax.XTickLabel=coherenceBands;
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P9-P13 - SST';
ax.YLim=[0 .8];

ax=subplot(2,3,3);
hold on
plot(data.PPC(BetaC&Nkx_y&P9P13,:)','-o','Color',[.5,.5,.5])
errorbar([1,2,3]+0.1,nanmean(data.PPC(BetaC&Nkx_y&P9P13,:)),nanstd(data.PPC(BetaC&Nkx_y&P9P13,:)),'k','LineWidth',2)
ax.XLim=[.5 3.5];
ax.YLabel.String='PPC';
ax.XTick=[1,2,3];
ax.XTickLabel=coherenceBands;
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P9-P13 - Nkx2-1';
ax.YLim=[0 .8];

ax=subplot(2,5,6);
hold on
plot(data.PPC(BetaC&RS&P14P18,:)','-o','Color',[.5,.5,.5])
errorbar([1,2,3]+0.1,nanmean(data.PPC(BetaC&RS&P14P18,:)),nanstd(data.PPC(BetaC&RS&P14P18,:)),'k','LineWidth',2)
ax.XLim=[.5 3.5];
ax.YLabel.String='PPC';
ax.XTick=[1,2,3];
ax.XTickLabel=coherenceBands;
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P14-P18 - RS';
ax.YLim=[0 .8];

ax=subplot(2,5,7);
hold on
plot(data.PPC(BetaC&SST&P14P18,:)','-o','Color',[.5,.5,.5])
errorbar([1,2,3]+0.1,nanmean(data.PPC(BetaC&SST&P14P18,:)),nanstd(data.PPC(BetaC&SST&P14P18,:)),'k','LineWidth',2)
ax.XLim=[.5 3.5];
ax.YLabel.String='PPC';
ax.XTick=[1,2,3];
ax.XTickLabel=coherenceBands;
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P14-P18 - SST';
ax.YLim=[0 .8];

ax=subplot(2,5,8);
hold on
plot(data.PPC(BetaC&FS&P14P18,:)','-o','Color',[.5,.5,.5])
errorbar([1,2,3]+0.1,nanmean(data.PPC(BetaC&FS&P14P18,:)),nanstd(data.PPC(BetaC&FS&P14P18,:)),'k','LineWidth',2)
ax.XLim=[.5 3.5];
ax.YLabel.String='PPC';
ax.XTick=[1,2,3];
ax.XTickLabel=coherenceBands;
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P14-P18 - FS';
ax.YLim=[0 .8];

ax=subplot(2,5,9);
hold on
plot(data.PPC(BetaC&Nkx_RS&P14P18,:)','-o','Color',[.5,.5,.5])
errorbar([1,2,3]+0.1,nanmean(data.PPC(BetaC&Nkx_RS&P14P18,:)),nanstd(data.PPC(BetaC&Nkx_RS&P14P18,:)),'k','LineWidth',2)
ax.XLim=[.5 3.5];
ax.YLabel.String='PPC';
ax.XTick=[1,2,3];
ax.XTickLabel=coherenceBands;
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P14-P18 - Nkx-RS';
ax.YLim=[0 .8];

ax=subplot(2,5,10);
hold on
plot(data.PPC(BetaC&Nkx_FS&P14P18,:)','-o','Color',[.5,.5,.5])
errorbar([1,2,3]+0.1,nanmean(data.PPC(BetaC&Nkx_FS&P14P18,:)),nanstd(data.PPC(BetaC&Nkx_FS&P14P18,:)),'k','LineWidth',2)
ax.XLim=[.5 3.5];
ax.YLabel.String='PPC';
ax.XTick=[1,2,3];
ax.XTickLabel=coherenceBands;
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P14-P18 - Nkx-FS';
ax.YLim=[0 .8];
% export_fig(fullfile(folderFigures,'Coherence_non-responsive_comparison'),'-tiff','-transparent','-nocrop')
% close
%% PPC
%All
figure('units','normalized','outerposition',[0 0 1 1]);
c=0;
for k=1:3
for j=1:3
    c=c+1;
    ax1=subplot(3,3,c);
    v=violinplot(data.PPC(agesLogical(:,k)&data.pValuePPC(:,j)<=0.05,j), removecats(data.cellIdentity(agesLogical(:,k)&data.pValuePPC(:,j)<=0.05)));
    for i=1:numel(v)
        v(i).ViolinColor=colours{i};
        v(i).ScatterPlot.MarkerFaceColor=[189,189,189]/255;
        v(i).ScatterPlot.MarkerFaceAlpha=1;
        v(i).ScatterPlot.SizeData=20;
    end
    ax1.FontSize=18;
    ax1.LineWidth=2;
    ax1.YLim=[0,.5];
    ax1.YLabel.String='Coherence';
    ax1.Title.String=strcat(developmentStage{k},' - ',coherenceBands{j});
    ax1.Title.FontSize=23;
    ax1.XTickLabelRotation=20;
    grid
end
end
export_fig(fullfile(folderFigures,'PPC_All'),'-tiff','-transparent','-nocrop')
close

%Only responsive
figure('units','normalized','outerposition',[0 0 1 1]);
c=0;
for k=2:3
for j=1:3
    c=c+1;
    ax1=subplot(2,3,c);
    v=violinplot(data.PPC(resp&agesLogical(:,k)&data.pValuePPC(:,j)<=0.05,j), removecats(data.cellIdentity(resp&agesLogical(:,k)&data.pValuePPC(:,j)<=0.05)));
    for i=1:numel(v)
        v(i).ViolinColor=colours{i};
        v(i).ScatterPlot.MarkerFaceColor=[189,189,189]/255;
        v(i).ScatterPlot.MarkerFaceAlpha=1;
        v(i).ScatterPlot.SizeData=20;
    end
    ax1.FontSize=18;
    ax1.LineWidth=2;
    ax1.YLim=[0,.5];
    ax1.YLabel.String='Coherence';
    ax1.Title.String=strcat(developmentStage{k},' - ',coherenceBands{j});
    ax1.Title.FontSize=23;
    ax1.XTickLabelRotation=20;
    grid
end
end
export_fig(fullfile(folderFigures,'PPC_responsive'),'-tiff','-transparent','-nocrop')
close
%% SU PSTH with opto average
plotLaserOnly=0;
figure('units','normalized','outerposition',[0 0 1 0.5]);

ax=subplot(1,3,1);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P5P8 & ~SST,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P5P8 & ~SST,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(PSTHbins,mean(Z_PSTHlaser(P5P8 & ~SST,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisualOpto(P5P8 & ~SST,:),1))-1 max(mean(Z_PSTHvisual(P5P8 & ~SST,:),1))+1];

ax=subplot(1,3,2);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P9P13 & ~SST,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P9P13 & ~SST,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(PSTHbins,mean(Z_PSTHlaser(P9P13 & ~SST,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisualOpto(P9P13 & ~SST,:),1))-1 max(mean(Z_PSTHvisual(P9P13 & ~SST,:),1))+1];


ax=subplot(1,3,3);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P14P18 & ~SST,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P14P18 & ~SST,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(PSTHbins,mean(Z_PSTHlaser(P14P18 & ~SST,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisualOpto(P14P18 & ~SST,:),1))-1 max(mean(Z_PSTHvisual(P14P18 & ~SST,:),1))+1];

legend('Control','+ SST opto - laser only')
figname='SU_ControlVsOpto_Summary';
% sg=sgtitle(strcat('P',int2str(ops.age)));
% sg.FontSize=30;
export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ON vs OFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
plotLaserOnly=0;
figure('units','normalized','outerposition',[0 0 1 1]);

ax=subplot(2,3,1);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P5P8 & ~SST & ON,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P5P8 & ~SST & ON,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(PSTHbins,mean(Z_PSTHlaser(P5P8 & ~SST & ON,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
if ~any(isnan([min(mean(Z_PSTHvisualOpto(P5P8 & ~SST & ON,:),1))-1 max(mean(Z_PSTHvisual(P5P8 & ~SST & ON,:),1))+1]))
    ax.YLim=[min(mean(Z_PSTHvisualOpto(P5P8 & ~SST & ON,:),1))-1 max(mean(Z_PSTHvisual(P5P8 & ~SST & ON,:),1))+1];
end
ax.XLabel.String='';


ax=subplot(2,3,2);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P9P13 & ~SST & ON,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P9P13 & ~SST & ON,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(PSTHbins,mean(Z_PSTHlaser(P9P13 & ~SST & ON,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisualOpto(P9P13 & ~SST & ON,:),1))-1 max(mean(Z_PSTHvisual(P9P13 & ~SST & ON,:),1))+1];
ax.XLabel.String='';

ax=subplot(2,3,3);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P14P18 & ~SST & ON,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P14P18 & ~SST & ON,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(PSTHbins,mean(Z_PSTHlaser(P14P18 & ~SST & ON,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisualOpto(P14P18 & ~SST & ON,:),1))-1 max(mean(Z_PSTHvisual(P14P18 & ~SST & ON,:),1))+1];
ax.XLabel.String='';
legend('Control','+ SST opto - laser only')

ax=subplot(2,3,4);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P5P8 & ~SST & OFF,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P5P8 & ~SST & OFF,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(PSTHbins,mean(Z_PSTHlaser(P5P8 & ~SST & OFF,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
if ~any(isnan([min(mean(Z_PSTHvisualOpto(P5P8 & ~SST & OFF,:),1))-1 max(mean(Z_PSTHvisual(P5P8 & ~SST & OFF,:),1))+1]))
    ax.YLim=[min(mean(Z_PSTHvisualOpto(P5P8 & ~SST & OFF,:),1))-1 max(mean(Z_PSTHvisual(P5P8 & ~SST & OFF,:),1))+1];
end

ax=subplot(2,3,5);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P9P13 & ~SST & OFF,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P9P13 & ~SST & OFF,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(PSTHbins,mean(Z_PSTHlaser(P9P13 & ~SST & OFF,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisualOpto(P9P13 & ~SST & OFF,:),1))-1 max(mean(Z_PSTHvisual(P9P13 & ~SST & OFF,:),1))+1];

ax=subplot(2,3,6);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P14P18 & ~SST & OFF,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P14P18 & ~SST & OFF,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(PSTHbins,mean(Z_PSTHlaser(P14P18 & ~SST & OFF,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisualOpto(P14P18 & ~SST & OFF,:),1))-1 max(mean(Z_PSTHvisual(P14P18 & ~SST & OFF,:),1))+1];



figname='SU_ControlVsOpto_ONvsOFF';
% sg=sgtitle(strcat('P',int2str(ops.age)));
% sg.FontSize=30;
export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% late vs nonresponsive %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
plotLaserOnly=0;
figure('units','normalized','outerposition',[0 0 1 1]);

ax=subplot(2,3,1);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P5P8 & ~SST & late,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P5P8 & ~SST & late,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(PSTHbins,mean(Z_PSTHlaser(P5P8 & ~SST & late,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
if ~any(isnan([min(mean(Z_PSTHvisualOpto(P5P8 & ~SST & late,:),1))-1 max(mean(Z_PSTHvisual(P5P8 & ~SST & late,:),1))+1]))
    ax.YLim=[min(mean(Z_PSTHvisualOpto(P5P8 & ~SST & late,:),1))-1 max(mean(Z_PSTHvisual(P5P8 & ~SST & late,:),1))+1];
end
ax.XLabel.String='';


ax=subplot(2,3,2);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P9P13 & ~SST & late,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P9P13 & ~SST & late,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(PSTHbins,mean(Z_PSTHlaser(P9P13 & ~SST & late,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisualOpto(P9P13 & ~SST & late,:),1))-1 max(mean(Z_PSTHvisual(P9P13 & ~SST & late,:),1))+1];
ax.XLabel.String='';

ax=subplot(2,3,3);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P14P18 & ~SST & late,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P14P18 & ~SST & late,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(PSTHbins,mean(Z_PSTHlaser(P14P18 & ~SST & late,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisualOpto(P14P18 & ~SST & late,:),1))-1 max(mean(Z_PSTHvisual(P14P18 & ~SST & late,:),1))+1];
ax.XLabel.String='';
legend('Control','+ SST opto - laser only')

ax=subplot(2,3,4);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P5P8 & ~SST & nonResponsive,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P5P8 & ~SST & nonResponsive,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(PSTHbins,mean(Z_PSTHlaser(P5P8 & ~SST & nonResponsive,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
if ~any(isnan([min(mean(Z_PSTHvisualOpto(P5P8 & ~SST & nonResponsive,:),1))-1 max(mean(Z_PSTHvisual(P5P8 & ~SST & nonResponsive,:),1))+1]))
    ax.YLim=[min(mean(Z_PSTHvisualOpto(P5P8 & ~SST & nonResponsive,:),1))-1 max(mean(Z_PSTHvisual(P5P8 & ~SST & nonResponsive,:),1))+1];
end

ax=subplot(2,3,5);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P9P13 & ~SST & nonResponsive,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P9P13 & ~SST & nonResponsive,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(PSTHbins,mean(Z_PSTHlaser(P9P13 & ~SST & nonResponsive,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisualOpto(P9P13 & ~SST & nonResponsive,:),1))-1 max(mean(Z_PSTHvisual(P9P13 & ~SST & nonResponsive,:),1))+1];

ax=subplot(2,3,6);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P14P18 & ~SST & nonResponsive,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P14P18 & ~SST & nonResponsive,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(PSTHbins,mean(Z_PSTHlaser(P14P18 & ~SST & nonResponsive,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisualOpto(P14P18 & ~SST & nonResponsive,:),1))-1 max(mean(Z_PSTHvisual(P14P18 & ~SST & nonResponsive,:),1))+1];



figname='SU_ControlVsOpto_lateVsNonResponsive';
% sg=sgtitle(strcat('P',int2str(ops.age)));
% sg.FontSize=30;
export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')

% Z_PSTHvisualSub=Z_PSTHvisualOpto-Z_PSTHlaser;
% 
% %With subtraction of laser only
% figure('units','normalized','outerposition',[0 0 .5 1]);
% 
% ax=subplot(3,1,1);
% patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
% hold on;
% plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
% plot(PSTHbins,mean(Z_PSTHvisual(P5P8,:),1),'k','LineWidth',2)
% plot(PSTHbins,mean(Z_PSTHvisualSub(P5P8,:),1),'r','LineWidth',2)
% ax=applyFont(ax,1);
% ax.YLim=[min(mean(Z_PSTHvisualSub(P5P8,:),1))-1 max(mean(Z_PSTHvisual(P5P8,:),1))+1];
% 
% ax=subplot(3,1,2);
% patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
% hold on;
% plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
% plot(PSTHbins,mean(Z_PSTHvisual(P9P13,:),1),'k','LineWidth',2)
% plot(PSTHbins,mean(Z_PSTHvisualSub(P9P13,:),1),'r','LineWidth',2)
% ax=applyFont(ax,1);
% ax.YLim=[min(mean(Z_PSTHvisualSub(P9P13,:),1))-1 max(mean(Z_PSTHvisual(P9P13,:),1))+1];
% 
% 
% ax=subplot(3,1,3);
% patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
% hold on;
% plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
% plot(PSTHbins,mean(Z_PSTHvisual(P14P18,:),1),'k','LineWidth',2)
% plot(PSTHbins,mean(Z_PSTHvisualSub(P14P18,:),1),'r','LineWidth',2)
% ax=applyFont(ax,1);
% ax.YLim=[min(mean(Z_PSTHvisualSub(P14P18,:),1))-1 max(mean(Z_PSTHvisual(P14P18,:),1))+1];
% 
% legend('Control','+ SST opto - laser only')








%% Optotagging
figure('units','normalized','outerposition',[0 0 1 0.5]);

ax=subplot(1,3,1);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot(PSTHbins,mean(Z_PSTHvisual(P5P8 & SST,:),1),'k','LineWidth',2)
% plot(PSTHbins,Z_PSTHvisual(P5P8 & SST,:),'Color',[.8 .8 .8],'LineWidth',2)
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisual(P5P8 & SST,:),1))-1 max(mean(Z_PSTHvisual(P5P8 & SST,:),1))+1];

ax=subplot(1,3,2);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot(PSTHbins,mean(Z_PSTHvisual(P9P13 & SST,:),1),'k','LineWidth',2)
% plot(PSTHbins,Z_PSTHvisual(P9P13 & SST,:),'Color',[.8 .8 .8],'LineWidth',2)
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisual(P9P13 & SST,:),1))-1 max(mean(Z_PSTHvisual(P9P13 & SST,:),1))+1];

ax=subplot(1,3,3);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
% plot(PSTHbins,Z_PSTHvisual(P14P18 & SST,:),'Color',[.8 .8 .8],'LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(P14P18 & SST,:),1),'k','LineWidth',2)
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisual(P14P18 & SST,:),1))-1 max(mean(Z_PSTHvisual(P14P18 & SST,:),1))+1];

% legend('Control','+ SST opto - laser only')
figname='SU_Optotagging_Zoom';
% sg=sgtitle(strcat('P',int2str(ops.age)));
% sg.FontSize=30;
export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')

age={'P5-P8','P9-P13','P14-P18'};
condition=[P5P8,P9P13,P14P18];
%%%%%%%%%%%%%%%%%%%%%%%%% Plot individual cells
for i = 1:3
    figure('units','normalized','outerposition',[0 0 1 1]);
    numUnits=size(Z_PSTHvisual(condition(:,i) & SST,:),1);
    for j=1:numUnits
        UnitsToPlot=Z_PSTHvisual(condition(:,i) & SST,:);
        ax=subplot(4,4,j);
        patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
        hold on;
        plot(PSTHbins,UnitsToPlot(j,:),'k','LineWidth',2)
        ax=applyFont(ax,1);
        ax.YLim=[min(UnitsToPlot(j,:))-1 max(UnitsToPlot(j,:))+1];
    end
    
    figname=strcat('SU_Optotagging_',age{i});
    sg=sgtitle(age{i});
    sg.FontSize=30;
    export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')

end

%% Chemogenetics S1
close all
clear 
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 
load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SingleUnitData.mat')
folderFigures='C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\WT Doctoral Programme in Neuroscience\DPhil in Neuroscience\DPhil thesis\Figures\Chapter 4';

%% Set group logic arrays
data=data(~isnan(data.endSlope),:);
data=data(data.Tagging~='SST;NrgKO',:);
data=data(data.brainArea=='S1BF',:);
data=data(data.state=='Urethane',:);

%Layer
tmpLayers(data.Layer==1,1)=categorical(cellstr('L2/3'));
tmpLayers(data.Layer==2,1)=categorical(cellstr('L4'));
tmpLayers(data.Layer==3,1)=categorical(cellstr('L5/6'));
data.Layer=tmpLayers;
clear tmpLayers
L23=data.Layer=='L2/3';
L4=data.Layer=='L4';
L56=data.Layer=='L5/6';

%Age
P5P8=data.Age<9;
P9P13=data.Age>=9 & data.Age<14;
P14P18=data.Age>=14;
data.Dev(P5P8)=categorical(cellstr('P5-P8'));
data.Dev(P9P13)=categorical(cellstr('P9-P13'));
data.Dev(P14P18)=categorical(cellstr('P14-P18'));

%Cell Type
SST=(data.Tagging=='SST' & data.responseTag==1);
Nkx=(data.Tagging=='Nkx2-1' & data.responseTag==1);
Untagged=(~SST & ~Nkx);

SSTopto = data.Tagging=='SST';
Nkxopto = data.Tagging=='Nkx2-1';
isLaser = ~(isnan(data.PSTHlaser(:,1)));

%% Chemogenetics V1
isChemo = ~isnan(data.responseVisual_K);

ChemoRealV1 = {'K18','K24','K35','K36','K38','K39','K43','K44','K45','K47','K6'};
ChemoRealS1 = {'K51','K52','K54','K55','K56','K58','K60','K62','K64','K71','K72','K73','K74'};
ChemoReal = [ChemoRealV1(:)',ChemoRealS1(:)'];
ChemoCtrlGFP = {'K82','K83','K84','K85','K86','K87','K88'}; 
ChemoCtrlSaline = {'K75','K76'};
ChemoNoExpression = {'K48','K4','K9','K10','K17','K20','K21','K23','K37','K41','K42','K46','K50','K53','K59'};

KORD=ismember(data.MouseID,ChemoReal);
GFP=ismember(data.MouseID,ChemoCtrlGFP);
Saline=ismember(data.MouseID,ChemoCtrlSaline);

%% Preprocess data
% [Z_PSTHvisual, responsiveVisual]=zscoreBaseline(data.PSTHvisual);
% Z_PSTHvisualOpto=zscoreBaseline(data.PSTHvisualOpto);
% data.Z_PSTHlaser=zscore(data.PSTHlaser,[],'all');
% Z_PSTHoptotagging=zscoreBaseline(data.PSTHoptotagging);
Z_PSTHvisual=zscore(data.PSTHvisual,[],'all');
% responsive=any(Z_PSTHvisual(:,PSTHbins>0 & PSTHbins<200)>=3,2);
resp1=data.responseWhisker==1;


data.peakVisualFast = max(data.PSTHvisual(:,(PSTHbins>0 & PSTHbins<200)),[],2);
data.peakVisualFast_N = data.peakVisualFast-data.rw_baselineFiring;

data.peakWhiskerFast = max(data.PSTHwhisker(:,(PSTHbins>0 & PSTHbins<200)),[],2);
data.peakWhiskerFast_N = data.peakWhiskerFast-data.baseline_firing;
data.peakWhiskerSecond = max(data.PSTHwhisker(:,(PSTHbins>500 & PSTHbins<700)),[],2);
data.PPR = (data.peakWhiskerSecond+1)./(data.peakWhiskerFast+1);

data.peakVisualFast_K = max(data.PSTHvisual_K(:,(PSTHbins>0 & PSTHbins<200)),[],2);
data.peakVisualFast_K_N = data.peakVisualFast_K-data.rw_baselineFiring;

data.peakWhiskerFast_K = max(data.PSTHwhisker_K(:,(PSTHbins>0 & PSTHbins<200)),[],2);
data.peakWhiskerFast_K_N = data.peakWhiskerFast_K-data.baseline_firing_K;
data.peakWhiskerSecond_K = max(data.PSTHwhisker_K(:,(PSTHbins>500 & PSTHbins<700)),[],2);
data.PPR_K = (data.peakWhiskerSecond_K+1)./(data.peakWhiskerFast_K+1);

data.visualOptoBaseline = mean(data.PSTHvisualOpto(:,(PSTHbins>-200 & PSTHbins<0)),2);
data.peakVisualOptoFast = max(data.PSTHvisualOpto(:,(PSTHbins>0 & PSTHbins<200)),[],2);
data.peakVisualOptoFast_N = data.peakVisualOptoFast-data.visualOptoBaseline;

data.optoBaseline = mean(data.PSTHlaser(:,(PSTHbins>-150 & PSTHbins<-50)),2);
data.meanOpto = mean(data.PSTHlaser(:,(PSTHbins>-50 & PSTHbins<50)),2);

data.rw_firing_N=data.rw_firing-data.rw_baselineFiring;
data.rw_firing_K_N=data.rw_firing_K-data.rw_baselineFiring_K;
data.rw_firing_N(data.rw_firing_N>10)=NaN;
data.rw_meanPSTH = mean(data.PSTHvisual(:,PSTHbins>400 & PSTHbins < 4000),2);
data.rwOpto_meanPSTH = mean(data.PSTHvisualOpto(:,PSTHbins>400 & PSTHbins < 4000),2);
data.rw_meanPSTH_K = mean(data.PSTHvisual_K(:,PSTHbins>400 & PSTHbins < 4000),2);

SB_entrained=data.sb_entrainedP<0.05;
SB_entrained_K=data.sb_entrainedP_K<0.05;
PPC_entrained=data.pValuePPC(:,2)<0.05;
PPC_entrained_K=data.pValuePPC_K(:,2)<0.05;

LaserInhibited = data.responseLaser==2;
LaserExcited = data.responseLaser==1;

data.peakVisualFast_Change = log2((data.peakVisualFast_K_N+1)./(data.peakVisualFast_N+1));
data.peakVisualSlow_Change = log2((data.rw_firing_K_N+1)./(data.rw_firing_N+1));
data.baseline_Change = log2((data.baseline_firing_K+1)./(data.baseline_firing+1));
data.fano_Change = log2(data.fanoWhisker_K./data.fanoWhisker);
data.sbSpikeProb_Change = log2((data.sb_spikeProb_K+1)./(data.sb_spikeProb+1));
data.PPC_Change = log2((data.PPC_K(:,2)+1)./(data.PPC(:,2)+1));

data.peakWhiskerFast_Change = log2((data.peakWhiskerFast_K_N+1)./(data.peakWhiskerFast_N+1));
data.PPR_Change = log2((data.PPR_K+1)./(data.PPR+1));

%% Re-set cell identity
FS=(Untagged & P14P18 & data.troughPeakTime<0.75);
RS=(Untagged & ~FS);
Nkx_FS=(Nkx & P14P18 & data.troughPeakTime<0.75);
Nkx_RS=(Nkx & P14P18 & ~Nkx_FS);
Nkx_y=(Nkx & (P5P8|P9P13));
SST=(data.Tagging=='SST' & data.responseTag==1);

data.cellIdentity(RS)=categorical(cellstr('RS'));
data.cellIdentity(FS)=categorical(cellstr('FS'));
data.cellIdentity(Nkx_FS)=categorical(cellstr('Nkx2-1 - FS'));
data.cellIdentity(Nkx_RS)=categorical(cellstr('Nkx2-1 - RS'));
data.cellIdentity(Nkx_y)=categorical(cellstr('Nkx2-1'));
data.cellIdentity(SST)=categorical(cellstr('SST'));

if ~(nnz(FS)+nnz(RS)+nnz(Nkx_FS)+nnz(Nkx_RS)+nnz(SST)+nnz(Nkx_y)==height(data)); error('Problem with cell identity'); end

%% Chemogenetic whisker plots
chemoYBar=[nnz(data.responseWhisker(KORD&P5P8))/nnz(KORD&P5P8),nnz(data.responseWhisker_K(KORD&P5P8))/nnz(KORD&P5P8);...
    nnz(data.responseWhisker(KORD&P9P13))/nnz(KORD&P9P13),nnz(data.responseWhisker_K(KORD&P9P13))/nnz(KORD&P9P13);...
    nnz(data.responseWhisker(KORD&P14P18))/nnz(KORD&P14P18),nnz(data.responseWhisker_K(KORD&P14P18))/nnz(KORD&P14P18)];

    chemoYBar2=[nnz(SB_entrained(KORD&P5P8))/nnz(KORD&P5P8),nnz(SB_entrained_K(KORD&P5P8))/nnz(KORD&P5P8);...
        nnz(PPC_entrained(KORD&P5P8))/nnz(KORD&P5P8),nnz(PPC_entrained_K(KORD&P5P8))/nnz(KORD&P5P8);...
            nnz(SB_entrained(KORD&P9P13))/nnz(KORD&P9P13),nnz(SB_entrained_K(KORD&P9P13))/nnz(KORD&P9P13);...
            nnz(PPC_entrained(KORD&P9P13))/nnz(KORD&P9P13),nnz(PPC_entrained_K(KORD&P9P13))/nnz(KORD&P9P13);...
            nnz(SB_entrained(KORD&P14P18))/nnz(KORD&P14P18),nnz(SB_entrained_K(KORD&P14P18))/nnz(KORD&P14P18);...
            nnz(PPC_entrained(KORD&P14P18))/nnz(KORD&P14P18),nnz(PPC_entrained_K(KORD&P14P18))/nnz(KORD&P14P18);];

figure('units','normalized','outerposition',[0 0 0.15 1]);
ax=subplot(3,1,1);
b = bar(chemoYBar);
for i=1:size(b,2)
    b(1,i).LineWidth=1;
end
ax.YLabel.String='Ratio single units';
ax.Box='off';
ax.LineWidth = 1;
ax.FontSize=10;
ax.YLim=[0,1];
legend('Control','KORD-SalB')
legend('boxoff')

ax=subplot(3,1,2);
b = bar(chemoYBar2);
for i=1:size(b,2)
    b(1,i).LineWidth=1;
end
ax.YLabel.String='Ratio single units';
ax.Box='off';
ax.LineWidth = 1;
ax.FontSize=10;
ax.YLim=[0,1];
legend('Control','KORD-SalB')
legend('boxoff')
export_fig(fullfile(folderFigures,'4.20','Barplots'),'-pdf','-transparent','-nocrop')
close

figure('units','normalized','outerposition',[0 0 0.15 1]);
ax(1)=subplot(3,1,1);
hold on
plot(ones(height(data(KORD&P5P8&resp1,:)),1)*1,data.peakWhiskerFast_Change(KORD&P5P8&resp1),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.peakWhiskerFast_Change(KORD&P5P8&resp1),'omitnan'),sem(data.peakWhiskerFast_Change(KORD&P5P8&resp1)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P9P13&resp1,:)),1)*4,data.peakWhiskerFast_Change(KORD&P9P13&resp1),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.peakWhiskerFast_Change(KORD&P9P13&resp1),'omitnan'),sem(data.peakWhiskerFast_Change(KORD&P9P13&resp1)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18&resp1,:)),1)*7,data.peakWhiskerFast_Change(KORD&P14P18&resp1),'o','Color',[150,150,150]/255)
errorbar(7.2,mean(data.peakWhiskerFast_Change(KORD&P14P18&resp1),'omitnan'),sem(data.peakWhiskerFast_Change(KORD&P14P18&resp1)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P5P8&resp1,:)),1)*2,data.peakWhiskerFast_Change(Saline&P5P8&resp1),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.peakWhiskerFast_Change(Saline&P5P8&resp1),'omitnan'),sem(data.peakWhiskerFast_Change(Saline&P5P8&resp1)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P9P13&resp1,:)),1)*5,data.peakWhiskerFast_Change(Saline&P9P13&resp1),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.peakWhiskerFast_Change(Saline&P9P13&resp1),'omitnan'),sem(data.peakWhiskerFast_Change(Saline&P9P13&resp1)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P14P18&resp1,:)),1)*8,data.peakWhiskerFast_Change(Saline&P14P18&resp1),'o','Color',[150,150,150]/255)
errorbar(8.2,mean(data.peakWhiskerFast_Change(Saline&P14P18&resp1),'omitnan'),sem(data.peakWhiskerFast_Change(Saline&P14P18&resp1)),'ro','Linewidth',1,'CapSize',10)
plot([0,8],[0,0],'k--','LineWidth',0.75)
ax(1).YAxis.Label.String='Max PSTH fast change (ratio log2)';

ax(2)=subplot(3,1,2);
hold on
plot(ones(height(data(KORD&P5P8&resp1,:)),1)*1,data.PPR_Change(KORD&P5P8&resp1),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.PPR_Change(KORD&P5P8&resp1),'omitnan'),sem(data.PPR_Change(KORD&P5P8&resp1)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P9P13&resp1,:)),1)*4,data.PPR_Change(KORD&P9P13&resp1),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.PPR_Change(KORD&P9P13&resp1),'omitnan'),sem(data.PPR_Change(KORD&P9P13&resp1)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18&resp1,:)),1)*7,data.PPR_Change(KORD&P14P18&resp1),'o','Color',[150,150,150]/255)
errorbar(7.2,mean(data.PPR_Change(KORD&P14P18&resp1),'omitnan'),sem(data.PPR_Change(KORD&P14P18&resp1)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P5P8&resp1,:)),1)*2,data.PPR_Change(Saline&P5P8&resp1),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.PPR_Change(Saline&P5P8&resp1),'omitnan'),sem(data.PPR_Change(Saline&P5P8&resp1)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P9P13&resp1,:)),1)*5,data.PPR_Change(Saline&P9P13&resp1),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.PPR_Change(Saline&P9P13&resp1),'omitnan'),sem(data.PPR_Change(Saline&P9P13&resp1)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P14P18&resp1,:)),1)*8,data.PPR_Change(Saline&P14P18&resp1),'o','Color',[150,150,150]/255)
errorbar(8.2,mean(data.PPR_Change(Saline&P14P18&resp1),'omitnan'),sem(data.PPR_Change(Saline&P14P18&resp1)),'ro','Linewidth',1,'CapSize',10)
plot([0,8],[0,0],'k--','LineWidth',0.75)
ax(2).YAxis.Label.String='PPR change (ratio log2)';


ax(3)=subplot(3,1,3);
hold on
plot(ones(height(data(KORD&P5P8,:)),1)*1,data.baseline_Change(KORD&P5P8),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.baseline_Change(KORD&P5P8),'omitnan'),sem(data.baseline_Change(KORD&P5P8)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P9P13,:)),1)*4,data.baseline_Change(KORD&P9P13),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.baseline_Change(KORD&P9P13),'omitnan'),sem(data.baseline_Change(KORD&P9P13)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*7,data.baseline_Change(KORD&P14P18),'o','Color',[150,150,150]/255)
errorbar(7.2,mean(data.baseline_Change(KORD&P14P18),'omitnan'),sem(data.baseline_Change(KORD&P14P18)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P5P8,:)),1)*2,data.baseline_Change(Saline&P5P8),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.baseline_Change(Saline&P5P8),'omitnan'),sem(data.baseline_Change(Saline&P5P8)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P9P13,:)),1)*5,data.baseline_Change(Saline&P9P13),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.baseline_Change(Saline&P9P13),'omitnan'),sem(data.baseline_Change(Saline&P9P13)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P14P18,:)),1)*8,data.baseline_Change(Saline&P14P18),'o','Color',[150,150,150]/255)
errorbar(8.2,mean(data.baseline_Change(Saline&P14P18),'omitnan'),sem(data.baseline_Change(Saline&P14P18)),'ro','Linewidth',1,'CapSize',10)
plot([0,8],[0,0],'k--','LineWidth',0.75)
ax(3).YAxis.Label.String='Baseline change (ratio log2)';

for i=1:numel(ax)
    ax(i).XLim=[0.5,8.5];
    ax(i).FontSize=10;
    ax(i).YLim=[-2,2];
end
export_fig(fullfile(folderFigures,'4.20','Errorbar responsivity'),'-pdf','-transparent','-nocrop')
close

figure('units','normalized','outerposition',[0 0 0.15 1]);
ax(1)=subplot(3,1,1);
hold on
plot(ones(height(data(KORD&P5P8&resp1,:)),1)*1,data.fano_Change(KORD&P5P8&resp1),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.fano_Change(KORD&P5P8&resp1),'omitnan'),sem(data.fano_Change(KORD&P5P8&resp1)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P9P13&resp1,:)),1)*4,data.fano_Change(KORD&P9P13&resp1),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.fano_Change(KORD&P9P13&resp1),'omitnan'),sem(data.fano_Change(KORD&P9P13&resp1)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18&resp1,:)),1)*7,data.fano_Change(KORD&P14P18&resp1),'o','Color',[150,150,150]/255)
errorbar(7.2,mean(data.fano_Change(KORD&P14P18&resp1),'omitnan'),sem(data.fano_Change(KORD&P14P18&resp1)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P5P8&resp1,:)),1)*2,data.fano_Change(Saline&P5P8&resp1),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.fano_Change(Saline&P5P8&resp1),'omitnan'),sem(data.fano_Change(Saline&P5P8&resp1)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P9P13&resp1,:)),1)*5,data.fano_Change(Saline&P9P13&resp1),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.fano_Change(Saline&P9P13&resp1),'omitnan'),sem(data.fano_Change(Saline&P9P13&resp1)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P14P18&resp1,:)),1)*8,data.fano_Change(Saline&P14P18&resp1),'o','Color',[150,150,150]/255)
errorbar(8.2,mean(data.fano_Change(Saline&P14P18&resp1),'omitnan'),sem(data.fano_Change(Saline&P14P18&resp1)),'ro','Linewidth',1,'CapSize',10)
plot([0,8],[0,0],'k--','LineWidth',0.75)
ax(1).YAxis.Label.String='Fano Factor change (ratio log2)';

ax(2)=subplot(3,1,2);
hold on
plot(ones(height(data(KORD&P5P8&SB_entrained,:)),1)*1,data.sbSpikeProb_Change(KORD&P5P8&SB_entrained),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.sbSpikeProb_Change(KORD&P5P8&SB_entrained),'omitnan'),sem(data.sbSpikeProb_Change(KORD&P5P8&SB_entrained)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P9P13&SB_entrained,:)),1)*4,data.sbSpikeProb_Change(KORD&P9P13&SB_entrained),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.sbSpikeProb_Change(KORD&P9P13&SB_entrained),'omitnan'),sem(data.sbSpikeProb_Change(KORD&P9P13&SB_entrained)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18&SB_entrained,:)),1)*7,data.sbSpikeProb_Change(KORD&P14P18&SB_entrained),'o','Color',[150,150,150]/255)
errorbar(7.2,mean(data.sbSpikeProb_Change(KORD&P14P18&SB_entrained),'omitnan'),sem(data.sbSpikeProb_Change(KORD&P14P18&SB_entrained)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P5P8&SB_entrained,:)),1)*2,data.sbSpikeProb_Change(Saline&P5P8&SB_entrained),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.sbSpikeProb_Change(Saline&P5P8&SB_entrained),'omitnan'),sem(data.sbSpikeProb_Change(Saline&P5P8&SB_entrained)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P9P13&SB_entrained,:)),1)*5,data.sbSpikeProb_Change(Saline&P9P13&SB_entrained),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.sbSpikeProb_Change(Saline&P9P13&SB_entrained),'omitnan'),sem(data.sbSpikeProb_Change(Saline&P9P13&SB_entrained)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P14P18&SB_entrained,:)),1)*8,data.sbSpikeProb_Change(Saline&P14P18&SB_entrained),'o','Color',[150,150,150]/255)
errorbar(8.2,mean(data.sbSpikeProb_Change(Saline&P14P18&SB_entrained),'omitnan'),sem(data.sbSpikeProb_Change(Saline&P14P18&SB_entrained)),'ro','Linewidth',1,'CapSize',10)
plot([0,8],[0,0],'k--','LineWidth',0.75)
ax(2).YAxis.Label.String='SB spike probability change (ratio log2)';

ax(3)=subplot(3,1,3);
hold on
plot(ones(height(data(KORD&P5P8&PPC_entrained,:)),1)*1,data.PPC_Change(KORD&P5P8&PPC_entrained),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.PPC_Change(KORD&P5P8&PPC_entrained),'omitnan'),sem(data.PPC_Change(KORD&P5P8&PPC_entrained)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P9P13&PPC_entrained,:)),1)*4,data.PPC_Change(KORD&P9P13&PPC_entrained),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.PPC_Change(KORD&P9P13&PPC_entrained),'omitnan'),sem(data.PPC_Change(KORD&P9P13&PPC_entrained)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18&PPC_entrained,:)),1)*7,data.PPC_Change(KORD&P14P18&PPC_entrained),'o','Color',[150,150,150]/255)
errorbar(7.2,mean(data.PPC_Change(KORD&P14P18&PPC_entrained),'omitnan'),sem(data.PPC_Change(KORD&P14P18&PPC_entrained)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P5P8&PPC_entrained,:)),1)*2,data.PPC_Change(Saline&P5P8&PPC_entrained),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.PPC_Change(Saline&P5P8&PPC_entrained),'omitnan'),sem(data.PPC_Change(Saline&P5P8&PPC_entrained)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P9P13&PPC_entrained,:)),1)*5,data.PPC_Change(Saline&P9P13&PPC_entrained),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.PPC_Change(Saline&P9P13&PPC_entrained),'omitnan'),sem(data.PPC_Change(Saline&P9P13&PPC_entrained)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P14P18&PPC_entrained,:)),1)*8,data.PPC_Change(Saline&P14P18&PPC_entrained),'o','Color',[150,150,150]/255)
errorbar(8.2,mean(data.PPC_Change(Saline&P14P18&PPC_entrained),'omitnan'),sem(data.PPC_Change(Saline&P14P18&PPC_entrained)),'ro','Linewidth',1,'CapSize',10)
plot([0,8],[0,0],'k--','LineWidth',0.75)
ax(3).YAxis.Label.String='PPC beta change (ratio log2)';


for i=1:numel(ax)
    ax(i).XLim=[0.5,8.5];
    ax(i).FontSize=10;
    ax(i).YLim=[-2,2];
end
export_fig(fullfile(folderFigures,'4.20','Errorbar spontaneous'),'-pdf','-transparent','-nocrop')
close
%%
figure('units','normalized','outerposition',[0 0 0.15 1]);
ax(1)=subplot(3,1,1);
hold on
plot(ones(height(data(KORD&P14P18&resp1&RS,:)),1),data.peakWhiskerFast_Change(KORD&P14P18&resp1&RS),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.peakWhiskerFast_Change(KORD&P14P18&resp1&RS),'omitnan'),sem(data.peakWhiskerFast_Change(KORD&P14P18&resp1&RS)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18&resp1&FS,:)),1)*4,data.peakWhiskerFast_Change(KORD&P14P18&resp1&FS),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.peakWhiskerFast_Change(KORD&P14P18&resp1&FS),'omitnan'),sem(data.peakWhiskerFast_Change(KORD&P14P18&resp1&FS)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P14P18&resp1&RS,:)),1)*2,data.peakWhiskerFast_Change(Saline&P14P18&resp1&RS),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.peakWhiskerFast_Change(Saline&P14P18&resp1&RS),'omitnan'),sem(data.peakWhiskerFast_Change(Saline&P14P18&resp1&RS)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P14P18&resp1&FS,:)),1)*5,data.peakWhiskerFast_Change(Saline&P14P18&resp1&FS),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.peakWhiskerFast_Change(Saline&P14P18&resp1&FS),'omitnan'),sem(data.peakWhiskerFast_Change(Saline&P14P18&resp1&FS)),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(1).YAxis.Label.String='Max PSTH fast change (ratio log2)';

ax(2)=subplot(3,1,2);
hold on
plot(ones(height(data(KORD&P14P18&resp1&RS,:)),1),data.fano_Change(KORD&P14P18&resp1&RS),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.fano_Change(KORD&P14P18&resp1&RS),'omitnan'),sem(data.fano_Change(KORD&P14P18&resp1&RS)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18&resp1&FS,:)),1)*4,data.fano_Change(KORD&P14P18&resp1&FS),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.fano_Change(KORD&P14P18&resp1&FS),'omitnan'),sem(data.fano_Change(KORD&P14P18&resp1&FS)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P14P18&resp1&RS,:)),1)*2,data.fano_Change(Saline&P14P18&resp1&RS),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.fano_Change(Saline&P14P18&resp1&RS),'omitnan'),sem(data.fano_Change(Saline&P14P18&resp1&RS)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P14P18&resp1&FS,:)),1)*5,data.fano_Change(Saline&P14P18&resp1&FS),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.fano_Change(Saline&P14P18&resp1&FS),'omitnan'),sem(data.fano_Change(Saline&P14P18&resp1&FS)),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(2).YAxis.Label.String='Fano factor change (ratio log2)';

ax(3)=subplot(3,1,3);
hold on
plot(ones(height(data(KORD&P14P18&RS,:)),1),data.baseline_Change(KORD&P14P18&RS),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(data.baseline_Change(KORD&P14P18&RS),'omitnan'),sem(data.baseline_Change(KORD&P14P18&RS)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18&FS,:)),1)*4,data.baseline_Change(KORD&P14P18&FS),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(data.baseline_Change(KORD&P14P18&FS),'omitnan'),sem(data.baseline_Change(KORD&P14P18&FS)),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P14P18&RS,:)),1)*2,data.baseline_Change(Saline&P14P18&RS),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(data.baseline_Change(Saline&P14P18&RS),'omitnan'),sem(data.baseline_Change(Saline&P14P18&RS)),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P14P18&FS,:)),1)*5,data.baseline_Change(Saline&P14P18&FS),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(data.baseline_Change(Saline&P14P18&FS),'omitnan'),sem(data.baseline_Change(Saline&P14P18&FS)),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(3).YAxis.Label.String='Baseline change (ratio log2)';

for i=1:numel(ax)
    ax(i).XLim=[0.5,5.5];
    ax(i).FontSize=10;
    ax(i).YLim=[-2,2];
end
export_fig(fullfile(folderFigures,'4.20','Errorbar RSvsFS'),'-pdf','-transparent','-nocrop')
close


%% Plotting for paper
data.cellIdentity(Nkx)=categorical(cellstr('Nkx2-1'));

clear ax v b
figure('units','normalized','outerposition',[0 0 0.2 .7])

ax(1)=subplot(3,2,1);
v(1,:)=violinplot(data.vectorAngle(P9P13&PPC_entrained & ((RS&L4)|SST|Nkx),2),data.cellIdentity(P9P13&PPC_entrained & ((RS&L4)|SST|Nkx)));
ax(1).YLim=[0, 360];
ax(1).YTick = [0,90,180,270,360];
ax(1).YGrid = 'on';
ax(1).YLabel.String = 'Vector Angle';

ax(2)=subplot(3,2,2);
v(2,:)=violinplot(data.vectorAngle(RS&P14P18&PPC_entrained,2),data.Layer(RS&P14P18&PPC_entrained));
ax(2).YLim=[0, 360];
ax(2).YTick = [0,90,180,270,360];
ax(2).YGrid = 'on';
ax(2).YLabel.String = 'Vector Angle';

ax(3)=subplot(3,2,3);
b(1,:)=boxchart(removecats(data.cellIdentity(P9P13&PPC_entrained & ((RS&L4)|SST|Nkx))),data.PPC(P9P13&PPC_entrained & ((RS&L4)|SST|Nkx),2),'MarkerStyle','none');
ax(3).YLim=[0, 0.6];
ax(3).YLabel.String = 'PPC';
% ax(5).YScale = 'log';

ax(4)=subplot(3,2,4);
b(2,:)=boxchart(removecats(data.cellIdentity((SST|Nkx)&P14P18&PPC_entrained)),data.PPC((SST|Nkx)&P14P18&PPC_entrained,2),'MarkerStyle','none');
ax(4).YLim=[0, .4];
ax(4).YLabel.String = 'PPC';


ax(5)=subplot(3,2,5);
b(3,:)=boxchart(data.Layer(RS&P9P13&PPC_entrained),data.PPC(RS&P9P13&PPC_entrained,2),'MarkerStyle','none');
ax(5).YLim=[0, 0.6];
ax(5).YLabel.String = 'PPC';
% ax(5).YScale = 'log';

ax(6)=subplot(3,2,6);
b(4,:)=boxchart(data.Layer(RS&P14P18&PPC_entrained),data.PPC(RS&P14P18&PPC_entrained,2),'MarkerStyle','none');
ax(6).YLim=[0, .6];
ax(6).YLabel.String = 'PPC';

for i=1:numel(ax)
    if ismember(i,[1:2])
        for j=1:size(v,2)
            v(i,j).ViolinColor=[100,100,100]/255;
            v(i,j).ScatterPlot.MarkerFaceColor=[37,37,37]/255;
            v(i,j).ScatterPlot.MarkerFaceAlpha=1;
            v(i,j).ScatterPlot.SizeData=2;
        end
    ax(i).XLim=[0.5,3.5];
    end
    
  
    ax(i).FontSize=12;
    ax(i).LineWidth=1;
    
end
print(gcf,'-dpdf',fullfile('C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\WT Doctoral Programme in Neuroscience\DPhil in Neuroscience\Manuscripts\V1 S1\Figures\Fig. 2\S1_PPC'))
% close

swtest(data.vectorAngle(P9P13&PPC_entrained & ((RS&L4)|SST|Nkx),2))
[p,tbl,stats]=kruskalwallis(data.vectorAngle(P9P13&PPC_entrained & ((RS&L4)|SST|Nkx),2),removecats(data.cellIdentity(P9P13&PPC_entrained & ((RS&L4)|SST|Nkx))));
multcompare(stats,'Dimension',[1,2],'CType' ,'dunn-sidak')

swtest(data.PPC(P9P13&PPC_entrained & RS,2))
[p,tbl,stats]=kruskalwallis(data.PPC(P9P13&PPC_entrained & RS,2),removecats(data.Layer(P9P13&PPC_entrained & RS)));
multcompare(stats,'Dimension',[1,2],'CType' ,'dunn-sidak')

[h,p,stats]=my_ttest2 (data.PPC(P9P13&PPC_entrained & SST,2),data.PPC(P9P13&PPC_entrained & Nkx,2));



%% 
clear ax b
figure('units','normalized','outerposition',[0 0 0.2 1])

yfast=[nnz(responsive&P5P8&RS)/nnz(P5P8&RS),nnz(responsive&P5P8&SST)/nnz(P5P8&SST),nnz(responsive&P5P8&Nkx)/nnz(P5P8&Nkx);...
    nnz(responsive&P9P13&RS)/nnz(P9P13&RS),nnz(responsive&P9P13&SST)/nnz(P9P13&SST),nnz(responsive&P9P13&Nkx)/nnz(P9P13&Nkx);...
    nnz(responsive&P14P18&RS)/nnz(P14P18&RS),nnz(responsive&P14P18&SST)/nnz(P14P18&SST),nnz(responsive&P14P18&Nkx)/nnz(P14P18&Nkx)];
% yslow=[nnz(data.rw_responsive&P9P13&RS)/nnz(P9P13&RS),nnz(data.rw_responsive&P9P13&SST)/nnz(P9P13&SST),nnz(data.rw_responsive&P9P13&Nkx)/nnz(P9P13&Nkx);...
%     nnz(data.rw_responsive&P14P18&RS)/nnz(P14P18&RS),nnz(data.rw_responsive&P14P18&SST)/nnz(P14P18&SST),nnz(data.rw_responsive&P14P18&Nkx)/nnz(P14P18&Nkx)];

xDev=reordercats(categorical (cellstr({'P5-P8','P9-P13','P14-P18'})),{'P5-P8','P9-P13','P14-P18'});

ax(1)=subplot(3,1,1);
b(1,:)=bar(xDev,yfast);
% 
% ax(2)=subplot(3,1,2);
% b(2,:)=bar(xDev,yslow);
% titles={'Fast','Slow'};
for i=1:numel(ax)
    ax(i).YLabel.String='Ratio single units responsive';
    ax(i).Box='off';
    ax(i).LineWidth = 1.5;
    ax(i).FontSize=15;
    ax(i).YTick=[0:0.2:1];
    legend('RS','SST','Nkx2-1','Location','northeast')
    legend('boxoff')
    for j=1:size(b,2)
        b(i,j).LineWidth=1;
    end
    ax(i).YLim=[0,1];
%     ax(i).Title.String=titles{i};
end
export_fig('C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\WT Doctoral Programme in Neuroscience\DPhil in Neuroscience\Manuscripts\V1 S1\Figures\Fig. 3\S1BF_RS_Responsive_BarPlot.pdf','-pdf','-transparent','-nocrop')
close