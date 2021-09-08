close all
clear 
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 
load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SingleUnitData.mat')
data.PSTHvisualOpto=[];
data.PSTHlaser=[];
data.responseVisualOpto=[];
data.responseLaser=[];
% data1=data;
% load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SingleUnitData_Part2.mat')
% data=[data1;data];
folderFigures='C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\Conferences\11th Annual Oxford Neuroscience Symposium\MyInVivo';

%% Set group logic arrays
data=data(~isnan(data.endSlope),:);
data=data(data.Tagging~='SST;NrgKO',:);
data=data(data.brainArea=='V1',:);
data=data(data.Tagging~='SST;NrgKO',:);



%Layer
tmpLayers(data.Layer==1,1)=categorical(cellstr('L2/3'));
tmpLayers(data.Layer==2,1)=categorical(cellstr('L4'));
tmpLayers(data.Layer==3,1)=categorical(cellstr('L5/6'));
data.Layer=tmpLayers;
%Age
P5P8=data.Age<9;
P9P13=data.Age>=9 & data.Age<14;
P14P18=data.Age>=14;

%Visual response
% responsive=any(data.responseVisual(:,1:2)==1,2);
responsive=any(data.responseVisual==1,2);


SST=(data.Tagging=='SST' & data.responseTag==1);
Nkx=(data.Tagging=='Nkx2-1' & data.responseTag==1);
Untagged=(~SST & ~Nkx);

% SST_E=responseLaser(:,1)==1;
% SST_I=responseLaser(:,1)==2;
% SST_Nonresponsive=responseLaser(:,1)==0;

%% Z-scoring SU
[Z_PSTHvisual, responsiveVisual]=zscoreBaseline(data.PSTHvisual);
% Z_PSTHvisualOpto=zscoreBaseline(data.PSTHvisualOpto);
% Z_PSTHlaser=zscoreBaseline(data.PSTHlaser);
Z_PSTHoptotagging=zscoreBaseline(data.PSTHoptotagging);


%%
Z=zscore(data.PSTHvisual,[],'all');
responsive=any(Z(:,PSTHbins>0 & PSTHbins<200)>=3,2);

for i=6:18
    fastPerc(i)=nnz(responsive(data.Age==i,:))/numel(responsive(data.Age==i,:));
    slowPerc(i)=nnz(data.rw_responsive(data.Age==i,:))/numel(data.rw_responsive(data.Age==i,:));
    bothPerc(i)=nnz(responsive(data.Age==i,:) & data.rw_responsive(data.Age==i,:))/numel(data.rw_responsive(data.Age==i,:));
end

%% Single Units identity
%P14-P18
figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle('P14-P18', 'FontSize',30)
ax=subplot(1,4,1);
[value,edges]=histcounts(data.troughPeakTime(P14P18&Untagged),0:0.1:3);
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
scatter(data.halfWidth(P14P18&Untagged),data.troughPeakTime(P14P18&Untagged),100,'k', 'filled','HandleVisibility','off')
scatter(data.halfWidth(P14P18&Nkx),data.troughPeakTime(P14P18&Nkx),100,'filled','r');
scatter(data.halfWidth(P14P18&SST),data.troughPeakTime(P14P18&SST),100,'filled','g');
xlim([0 1])
ylim([0 2.6])
ax.XLabel.String='Spike Half-Width (ms)';
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
ax.YAxis.Visible='off';
grid
legend('Nkx2-1','SST')
% export_fig(fullfile(folderFigures,'WF_scatter_P14P18'),'-tiff','-transparent','-nocrop')
% close

% <P14
figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle('P8-P13', 'FontSize',30)
ax=subplot(1,4,1);
[value,edges]=histcounts(data.troughPeakTime((P5P8|P9P13)&Untagged),0:0.1:3);
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
scatter(data.halfWidth((P5P8|P9P13)&Untagged),data.troughPeakTime((P5P8|P9P13)&Untagged),100,'k', 'filled','HandleVisibility','off')
hold on
scatter(data.halfWidth((P5P8|P9P13)&Nkx),data.troughPeakTime((P5P8|P9P13)&Nkx),100,'filled','r');
scatter(data.halfWidth((P5P8|P9P13)&SST),data.troughPeakTime((P5P8|P9P13)&SST),100,'filled','g');
xlim([0 1])
ylim([0 2.6])
ax.XLabel.String='Spike Half-Width (ms)';
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
ax.YAxis.Visible='off';
legend('Nkx2-1','SST')
grid
% export_fig(fullfile(folderFigures,'WF_scatter_P8P13'),'-tiff','-transparent','-nocrop')
% close


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
SST=(data.Tagging=='SST' & data.responseTag==1);
Nkx=(data.Tagging=='Nkx2-1' & data.responseTag==1);
Untagged=(~SST & ~Nkx);


FS=(Untagged & P14P18 & data.troughPeakTime<0.7);
RS=(Untagged & ~FS);
Nkx_FS=(Nkx & P14P18 & data.troughPeakTime<0.7);
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

%% All units summary imagesc

SummaryTableVisual=[nnz(P5P8 & ON) nnz(P9P13 & ON) nnz(P14P18 & ON);...
    nnz(P5P8 & OFF) nnz(P9P13 & OFF) nnz(P14P18 & OFF);...
    nnz(P5P8 & late) nnz(P9P13 & late) nnz(P14P18 & late);...
    nnz(P5P8 & nonResponsive) nnz(P9P13 & nonResponsive) nnz(P14P18 & nonResponsive);...
    nnz(P5P8 & ONOFF) nnz(P9P13 & ONOFF) nnz(P14P18 & ONOFF);...
    nnz(P5P8 & ONlate) nnz(P9P13 & ONlate) nnz(P14P18 & ONlate);...
    nnz(P5P8 & OFFlate) nnz(P9P13 & OFFlate) nnz(P14P18 & OFFlate);...
    nnz(P5P8 & ONOFFlate) nnz(P9P13 & ONOFFlate) nnz(P14P18 & ONOFFlate)];


SummaryTableSST=[nnz(P5P8 & SST_E) nnz(P9P13 & SST_E) nnz(P14P18 & SST_E);...
    nnz(P5P8 & SST_I) nnz(P9P13 & SST_I) nnz(P14P18 & SST_I);...
    nnz(P5P8 & SST_Nonresponsive) nnz(P9P13 & SST_Nonresponsive) nnz(P14P18 & SST_Nonresponsive);...
    nnz(P5P8 & SST) nnz(P9P13 & SST) nnz(P14P18 & SST)];

figure
[b,i]=sort(data.Age);
PSTHtoplot=Z_PSTHvisual(i,:);
imagesc('XData',PSTHbins,'YData',1:numel(data.suid),'CData',PSTHtoplot,[0 40])
hold on
plot([min(PSTHbins) max(PSTHbins)],[nnz(P5P8) nnz(P5P8)],'w--','LineWidth',2)
plot([min(PSTHbins) max(PSTHbins)],[nnz(P5P8)+nnz(P9P13) nnz(P5P8)+nnz(P9P13)],'w--','LineWidth',2)
ax=gca;
ax.YDir='reverse';
colormap('hot')
ax.XLim=[-500 5000];
ax.YLim=[1 numel(data.suid)];
ax.XLabel.String='Time (ms)';
ax.YLabel.String=('Unit #');
ax.FontSize=20;

% figname='SU_Imagesc';
% sg=sgtitle(strcat('P',int2str(ops.age)));
% sg.FontSize=30;
% export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')

%% Responsive units stats
% resp=data.rw_responsive;
resp=responsive;
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

% export_fig(fullfile(folderFigures,'PercentageResponsiveCells'),'-tiff','-transparent','-nocrop')
% close

%% Visual response stats
[MaxPSTH_Val,MaxPSTH_Idx]=max(data.PSTHvisual(:,PSTHbins>0&PSTHbins<3000),[],2);
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

export_fig(fullfile(folderFigures,'ViolinPlotPSTHAmplitudeLatency'),'-tiff','-transparent','-nocrop')
ax3.YLim=[0,400];
ax4.YLim=[0,400];
export_fig(fullfile(folderFigures,'ViolinPlotPSTHAmplitudeLatency_Zoom'),'-tiff','-transparent','-nocrop')
close

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

%SST young
subset=data(resp & SST & P9P13,:);
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
sgtitle('SST - P9-P13 - Responsive','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-SST-P9P13-Responsive'),'-tiff','-transparent','-nocrop')
close

subset=data(~resp & SST & P9P13,:);
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
sgtitle('Nkx2-1 - P9-P13 - Responsive','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-Nkx21-P9P13-Responsive'),'-tiff','-transparent','-nocrop')
close

subset=data(~resp & Nkx_y & P9P13,:);
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
sgtitle('SST - P14-P18 - Responsive','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-SST-P14P18-Responsive'),'-tiff','-transparent','-nocrop')
close

subset=data(~resp & SST & P14P18,:);
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
sgtitle('Nkx-FS - P14-P18 - Responsive','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-Nkx_FS-P14P18-Responsive'),'-tiff','-transparent','-nocrop')
close

subset=data(~resp & Nkx_FS & P14P18,:);
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
sgtitle('Nkx-RS - P14-P18 - Responsive','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-Nkx_RS-P14P18-Responsive'),'-tiff','-transparent','-nocrop')
close

subset=data(~resp & Nkx_RS & P14P18,:);
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
        sgtitle('FS - P14-P18 - Responsive','FontSize',25)
        export_fig(fullfile(folderFigures,strcat('PSTH-FS-P14P18-Responsive-',int2str(f))),'-tiff','-transparent','-nocrop')
        close
        c=0; 
    end
end
if c~=0
    export_fig(fullfile(folderFigures,strcat('PSTH-FS-P14P18-Responsive-',int2str(f+1))),'-tiff','-transparent','-nocrop')
    close
end

subset=data(~resp & FS & P14P18,:);
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
sgtitle('FS - P14-P18 - Non-responsive','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH-FS-P14P18-NonResponsive'),'-tiff','-transparent','-nocrop')
close

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
        sgtitle('RS - P9-P13 - Responsive','FontSize',25)
        export_fig(fullfile(folderFigures,strcat('PSTH-RS-P9P13-Responsive-',int2str(f))),'-tiff','-transparent','-nocrop')
        close
        c=0; 
    end
end
if c~=0
    sgtitle('RS - P9-P13 - Responsive','FontSize',25)
    export_fig(fullfile(folderFigures,strcat('PSTH-RS-P9P13-Responsive-',int2str(f+1))),'-tiff','-transparent','-nocrop')
    close
end

subset=data(~resp & RS & P9P13,:);
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
    sgtitle('RS - P14-P18 - Responsive','FontSize',25)
    export_fig(fullfile(folderFigures,strcat('PSTH-RS-P14P18-Responsive-',int2str(f+1))),'-tiff','-transparent','-nocrop')
    close
end

subset=data(~resp & RS & P14P18,:);
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

%SST
figure('units','normalized','outerposition',[0 0 1 0.5]);
ax=subplot(1,2,1);
imagesc(PSTHbins,1:height(data(resp&SST&P9P13,:)),sortrows(data(resp&SST&P9P13,:),{'Layer','MaxPSTH_Idx'}).PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp&SST&P9P13,:).Layer=='L2/3'),nnz(data(resp&SST&P9P13,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp&SST&P9P13,:).Layer=='L5/6')),nnz(not(data(resp&SST&P9P13,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-100,400];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P9-P13';
ax.CLim=[0,20];

ax=subplot(1,2,2);
imagesc(PSTHbins,1:height(data(resp&SST&P14P18,:)),sortrows(data(resp&SST&P14P18,:),{'Layer','MaxPSTH_Idx'}).PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp&SST&P14P18,:).Layer=='L2/3'),nnz(data(resp&SST&P14P18,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp&SST&P14P18,:).Layer=='L5/6')),nnz(not(data(resp&SST&P14P18,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-100,400];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P14-P18';

sgtitle('SST','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH_img_SST'),'-tiff','-pdf','-transparent','-nocrop')
close
    
%Nkx
figure('units','normalized','outerposition',[0 0 1 0.5]);
ax=subplot(1,2,1);
imagesc(PSTHbins,1:height(data(resp&Nkx&P9P13,:)),sortrows(data(resp&Nkx&P9P13,:),{'Layer','MaxPSTH_Idx'}).PSTHvisual)
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
imagesc(PSTHbins,1:height(data(resp&Nkx&P14P18,:)),sortrows(data(resp&Nkx&P14P18,:),{'Layer','MaxPSTH_Idx'}).PSTHvisual)
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
export_fig(fullfile(folderFigures,'PSTH_img_Nkx2-1'),'-tiff','-pdf','-transparent','-nocrop')
close

%RS
figure('units','normalized','outerposition',[0 0 1 1]);
ax=subplot(1,2,1);
imagesc(PSTHbins,1:height(data(resp&RS&P9P13,:)),sortrows(data(resp&RS&P9P13,:),{'Layer','MaxPSTH_Idx'}).PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp&RS&P9P13,:).Layer=='L2/3'),nnz(data(resp&RS&P9P13,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp&RS&P9P13,:).Layer=='L5/6')),nnz(not(data(resp&RS&P9P13,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-100,400];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P9-P13';

ax=subplot(1,2,2);
imagesc(PSTHbins,1:height(data(resp&RS&P14P18,:)),sortrows(data(resp&RS&P14P18,:),{'Layer','MaxPSTH_Idx'}).PSTHvisual)
colormap('hot')
colorbar
hold on
plot([-1000,5000],[nnz(data(resp&RS&P14P18,:).Layer=='L2/3'),nnz(data(resp&RS&P14P18,:).Layer=='L2/3')]+0.5,'w--','LineWidth',2)
plot([-1000,5000],[nnz(not(data(resp&RS&P14P18,:).Layer=='L5/6')),nnz(not(data(resp&RS&P14P18,:).Layer=='L5/6'))]+0.5,'w--','LineWidth',2)
ax.XLim=[-100,400];
ax.YAxis.Visible='off';
ax.XLabel.String='Time (ms)';
ax.FontSize=18;
ax.Title.String='P14-P18';

sgtitle('RS','FontSize',25)
export_fig(fullfile(folderFigures,'PSTH_img_RS'),'-tiff','-pdf','-transparent','-nocrop')
close

%FS
figure('units','normalized','outerposition',[0 0 1 1]);
ax=subplot(1,2,2);
imagesc(PSTHbins,1:height(data(resp&FS&P14P18,:)),sortrows(data(resp&FS&P14P18,:),{'Layer','MaxPSTH_Idx'}).PSTHvisual)
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
export_fig(fullfile(folderFigures,'PSTH_img_FS'),'-tiff','-pdf','-transparent','-nocrop')
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
plot(data.coherence(resp&RS&P9P13,:)','-o','Color',[.5,.5,.5])
errorbar([1,2,3]+0.1,nanmean(data.coherence(resp&RS&P9P13,:)),nanstd(data.coherence(resp&RS&P9P13,:)),'k','LineWidth',2)
ax.XLim=[.5 3.5];
ax.YLabel.String='Coherence';
ax.XTick=[1,2,3];
ax.XTickLabel=coherenceBands;
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P9-P13 - RS';
ax.YLim=[0 .8];

ax=subplot(2,3,2);
hold on
plot(data.coherence(resp&SST&P9P13,:)','-o','Color',[.5,.5,.5])
errorbar([1,2,3]+0.1,nanmean(data.coherence(resp&SST&P9P13,:)),nanstd(data.coherence(resp&SST&P9P13,:)),'k','LineWidth',2)
ax.XLim=[.5 3.5];
ax.YLabel.String='Coherence';
ax.XTick=[1,2,3];
ax.XTickLabel=coherenceBands;
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P9-P13 - SST';
ax.YLim=[0 .8];

ax=subplot(2,3,3);
hold on
plot(data.coherence(resp&Nkx_y&P9P13,:)','-o','Color',[.5,.5,.5])
errorbar([1,2,3]+0.1,nanmean(data.coherence(resp&Nkx_y&P9P13,:)),nanstd(data.coherence(resp&Nkx_y&P9P13,:)),'k','LineWidth',2)
ax.XLim=[.5 3.5];
ax.YLabel.String='Coherence';
ax.XTick=[1,2,3];
ax.XTickLabel=coherenceBands;
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P9-P13 - Nkx2-1';
ax.YLim=[0 .8];

ax=subplot(2,5,6);
hold on
plot(data.coherence(resp&RS&P14P18,:)','-o','Color',[.5,.5,.5])
errorbar([1,2,3]+0.1,nanmean(data.coherence(resp&RS&P14P18,:)),nanstd(data.coherence(resp&RS&P14P18,:)),'k','LineWidth',2)
ax.XLim=[.5 3.5];
ax.YLabel.String='Coherence';
ax.XTick=[1,2,3];
ax.XTickLabel=coherenceBands;
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P14-P18 - RS';
ax.YLim=[0 .8];

ax=subplot(2,5,7);
hold on
plot(data.coherence(resp&SST&P14P18,:)','-o','Color',[.5,.5,.5])
errorbar([1,2,3]+0.1,nanmean(data.coherence(resp&SST&P14P18,:)),nanstd(data.coherence(resp&SST&P14P18,:)),'k','LineWidth',2)
ax.XLim=[.5 3.5];
ax.YLabel.String='Coherence';
ax.XTick=[1,2,3];
ax.XTickLabel=coherenceBands;
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P14-P18 - SST';
ax.YLim=[0 .8];

ax=subplot(2,5,8);
hold on
plot(data.coherence(resp&FS&P14P18,:)','-o','Color',[.5,.5,.5])
errorbar([1,2,3]+0.1,nanmean(data.coherence(resp&FS&P14P18,:)),nanstd(data.coherence(resp&FS&P14P18,:)),'k','LineWidth',2)
ax.XLim=[.5 3.5];
ax.YLabel.String='Coherence';
ax.XTick=[1,2,3];
ax.XTickLabel=coherenceBands;
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P14-P18 - FS';
ax.YLim=[0 .8];

ax=subplot(2,5,9);
hold on
plot(data.coherence(resp&Nkx_RS&P14P18,:)','-o','Color',[.5,.5,.5])
errorbar([1,2,3]+0.1,nanmean(data.coherence(resp&Nkx_RS&P14P18,:)),nanstd(data.coherence(resp&Nkx_RS&P14P18,:)),'k','LineWidth',2)
ax.XLim=[.5 3.5];
ax.YLabel.String='Coherence';
ax.XTick=[1,2,3];
ax.XTickLabel=coherenceBands;
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P14-P18 - Nkx-RS';
ax.YLim=[0 .8];

ax=subplot(2,5,10);
hold on
plot(data.coherence(resp&Nkx_FS&P14P18,:)','-o','Color',[.5,.5,.5])
errorbar([1,2,3]+0.1,nanmean(data.coherence(resp&Nkx_FS&P14P18,:)),nanstd(data.coherence(resp&Nkx_FS&P14P18,:)),'k','LineWidth',2)
ax.XLim=[.5 3.5];
ax.YLabel.String='Coherence';
ax.XTick=[1,2,3];
ax.XTickLabel=coherenceBands;
ax.FontSize=18;
ax.LineWidth=2;
ax.Title.String='P14-P18 - Nkx-FS';
ax.YLim=[0 .8];
export_fig(fullfile(folderFigures,'Coherence_non-responsive_comparison'),'-tiff','-transparent','-nocrop')
close
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


    

%% Searching for PV
%Trough Peak Time vs Spike Half-width
figure('units','normalized','outerposition',[0 0 1 1]);
ax=subplot(1,4,1);
[value,edges]=histcounts(data.troughPeakTime,0:0.1:3);
h=barh(edges(1:end-1),value,'w','LineWidth',2);
set(get(h,'Parent'),'xdir','r')
ylim([0 2.2])
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.XLabel.String='Cell #';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
grid

ax=subplot(1,4,(2:4));
scatter(data.halfWidth(P14P18),data.troughPeakTime(P14P18),100,'k', 'filled','HandleVisibility','off')
hold on
scatter(data.halfWidth(P9P13),data.troughPeakTime(P9P13),100,'filled','r');
scatter(data.halfWidth(P5P8),data.troughPeakTime(P5P8),100,'filled','g');
xlim([0 1.1])
ylim([0 2.2])
ax.XLabel.String='Spike Half-Width (ms)';
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
ax.YAxis.Visible='off';
grid

%End slope peak troguh ratio
figure('units','normalized','outerposition',[0 0 1 1]);
ax=subplot(1,4,1);
[value,edges]=histcounts(data.endSlope);
h=barh(edges(1:end-1),value,'w','LineWidth',2);
set(get(h,'Parent'),'xdir','r')
ylim([-17 1])
ax.YLabel.String='End Slope';
ax.XLabel.String='Cell #';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
grid

ax=subplot(1,4,(2:4));
scatter(data.peakTroughRatio(P14P18),data.endSlope(P14P18),100,'k', 'filled','HandleVisibility','off')
hold on
scatter(data.peakTroughRatio(P9P13),data.endSlope(P9P13),100,'filled','r');
scatter(data.peakTroughRatio(P5P8),data.endSlope(P5P8),100,'filled','g');
xlim([-0.4 0.6])
ylim([-17 1])
ax.XLabel.String='Peak Trough Height Ratio';
ax.YLabel.String='End Slope';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
ax.YAxis.Visible='off';
grid

%Trough Peak Time vs peak troguh ratio
figure('units','normalized','outerposition',[0 0 1 1]);
ax=subplot(1,4,1);
[value,edges]=histcounts(data.troughPeakTime,0:0.1:3);
h=barh(edges(1:end-1),value,'w','LineWidth',2);
set(get(h,'Parent'),'xdir','r')
ylim([0 2.2])
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.XLabel.String='Cell #';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
grid

ax=subplot(1,4,(2:4));
scatter(data.peakTroughRatio(P14P18),data.troughPeakTime(P14P18),100,'k', 'filled','HandleVisibility','off')
hold on
scatter(data.peakTroughRatio(P9P13),data.troughPeakTime(P9P13),100,'filled','r');
scatter(data.peakTroughRatio(P5P8),data.troughPeakTime(P5P8),100,'filled','g');
xlim([-0.4 0.6])
ylim([0 2.2])
ax.XLabel.String='Peak Trough Height Ratio';
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
ax.YAxis.Visible='off';
grid

figname=strcat('WaveformComponents');
export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')





kMeansTab=[halfWidth(P14P18),troughPeakTime(P14P18)];
kMeansIdx = kmeans(kMeansTab,2,'Replicates',100);


ax=subplot(1,4,(2:4));
scatter(kMeansTab(kMeansIdx==1,1),kMeansTab(kMeansIdx==1,2),'k','filled')
hold on
scatter(kMeansTab(kMeansIdx==2,1),kMeansTab(kMeansIdx==2,2),'r','filled')

figname=strcat('WaveformComponents_KmeansClustering');
export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')
% hold on
% scatter(halfWidth(P9P13),troughPeakTime(P9P13),'filled','b');
% scatter(halfWidth(P5P8),troughPeakTime(P5P8),'filled','g');
% ax=get(gca);
% ax.XLim=[0 1.1];
% ax.YLim=[0 2.2];
% ax.XLabel.String='Spike Half-Width (ms)';
% ax.YLabel.String='Trough to Peak Latency (ms)';
% ax.FontSize=20;
ax=subplot(1,4,(2:4));
hold on
scatter(halfWidth(P9P13),troughPeakTime(P9P13),100,'filled','c');
scatter(halfWidth(P5P8),troughPeakTime(P5P8),100,'filled','g');
legend('P14-P18 - RS','P14-P18 - FS','P9-P13','P5-P8')

figname=strcat('WaveformComponents_KmeansClustering+Young');
export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')

ax=subplot(1,4,(2:4));
hold on
scatter(halfWidth(~SST),troughPeakTime(~SST),100,'filled','k');
scatter(halfWidth(SST),troughPeakTime(SST),100,'filled','r');
legend('Non-SST','SST')

figname=strcat('WaveformComponents_SST');
export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')


%DBSCAN clustering


dbscanIdx = dbscan(kMeansTab,0.04,10);

figure
scatter(kMeansTab(dbscanIdx==1,1),kMeansTab(dbscanIdx==1,2),'k','filled')
hold on
scatter(kMeansTab(dbscanIdx==-1,1),kMeansTab(dbscanIdx==-1,2),'m','filled')
scatter(kMeansTab(dbscanIdx==2,1),kMeansTab(dbscanIdx==2,2),'r','filled')
scatter(kMeansTab(dbscanIdx==3,1),kMeansTab(dbscanIdx==3,2),'g','filled')
scatter(kMeansTab(dbscanIdx==4,1),kMeansTab(dbscanIdx==4,2),'c','filled')
scatter(kMeansTab(dbscanIdx==5,1),kMeansTab(dbscanIdx==5,2),'b','filled')
scatter(kMeansTab(dbscanIdx==6,1),kMeansTab(dbscanIdx==6,2),'y','filled')

kD = pdist2(kMeansTab,kMeansTab,'euc','Smallest',10); % The minpts smallest pairwise distances
figure
plot(sort(kD(end,:)));
title('k-distance graph')
xlabel('Points sorted with 50th nearest distances')
ylabel('50th nearest distances')
grid

%% Fast spiking vs. regular spiking
FS=troughPeakTime<=0.7;
nnz(FS&P5P8)

figure
scatter3(halfWidth(~FS),troughPeakTime(~FS),peakTroughRatio(~FS),'k')
hold on
scatter3(halfWidth(FS),troughPeakTime(FS),peakTroughRatio(FS),'r')

%Trough Peak Time vs Spike Half-width
figure('units','normalized','outerposition',[0 0 1 1]);
ax=subplot(1,4,1);
[value,edges]=histcounts(troughPeakTime,0:0.1:3);
h=barh(edges(1:end-1),value,'w','LineWidth',2);
set(get(h,'Parent'),'xdir','r')
ylim([0 2.2])
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.XLabel.String='Cell #';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
grid

ax=subplot(1,4,(2:4));
scatter(halfWidth(~FSPutative),troughPeakTime(~FSPutative),100,'k', 'filled','HandleVisibility','off')
hold on
scatter(halfWidth(FSPutative),troughPeakTime(FSPutative),100,'filled','r');
xlim([0 1.1])
ylim([0 2.2])
ax.XLabel.String='Spike Half-Width (ms)';
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
ax.YAxis.Visible='off';
grid

%Trough Peak Time vs peak troguh ratio
figure('units','normalized','outerposition',[0 0 1 1]);
ax=subplot(1,4,1);
[value,edges]=histcounts(troughPeakTime,0:0.1:3);
h=barh(edges(1:end-1),value,'w','LineWidth',2);
set(get(h,'Parent'),'xdir','r')
ylim([0 2.2])
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.XLabel.String='Cell #';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
grid

ax=subplot(1,4,(2:4));
scatter(peakTroughRatio(~FS),troughPeakTime(~FS),100,'k', 'filled','HandleVisibility','off')
hold on
scatter(peakTroughRatio(FS),troughPeakTime(FS),100,'filled','r');
xlim([-0.4 0.6])
ylim([0 2.2])
ax.XLabel.String='Peak Trough Height Ratio';
ax.YLabel.String='Trough to Peak Latency (ms)';
ax.FontSize=25;
ax.Box='off';
ax.LineWidth = 1.5;
ax.YAxis.Visible='off';
grid

clos%%
[reduction,umap]=run_umap([halfWidth,troughPeakTime,peakTroughRatio]);

score=tsne([halfWidth,troughPeakTime,peakTroughRatio],'Distance','euclidean');
figure
hold on
scatter(score(~FS,1),score(~FS,2),'filled','k')
scatter(score(FS,1),score(FS,2),'filled','k')


kClu=kmeans([halfWidth,troughPeakTime,peakTroughRatio],2,'Distance','sqeuclidean','Replicates',1000);
figure
hold on
scatter(score(kClu==1,1),score(kClu==1,2),'filled','k')
scatter(score(kClu==2,1),score(kClu==2,2),'filled','r')
%% 

figure('units','normalized','outerposition',[0 0 1 0.5]);

ax=subplot(1,3,1);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot(PSTHbins,mean(Z_PSTHvisual(P5P8 & ~FS,:),1),'k','LineWidth',2)
% plot(PSTHbins,Z_PSTHvisual(P5P8 & SST,:),'Color',[.8 .8 .8],'LineWidth',2)
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisual(P5P8 & ~FS,:),1))-1 max(mean(Z_PSTHvisual(P5P8 & ~FS,:),1))+1];

ax=subplot(1,3,2);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot(PSTHbins,mean(Z_PSTHvisual(P9P13 & ~FS,:),1),'k','LineWidth',2)
% plot(PSTHbins,Z_PSTHvisual(P9P13 & SST,:),'Color',[.8 .8 .8],'LineWidth',2)
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisual(P9P13 & ~FS,:),1))-1 max(mean(Z_PSTHvisual(P9P13 & ~FS,:),1))+1];

ax=subplot(1,3,3);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
% plot(PSTHbins,Z_PSTHvisual(P14P18 & SST,:),'Color',[.8 .8 .8],'LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(P14P18 & ~FS,:),1),'k','LineWidth',2)
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisual(P14P18 & ~FS,:),1))-1 max(mean(Z_PSTHvisual(P14P18 & ~FS,:),1))+1];

% legend('Control','+ SST opto - laser only')
figname='SU_RS';
% sg=sgtitle(strcat('P',int2str(ops.age)));
% sg.FontSize=30;
export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')


figure('units','normalized','outerposition',[0 0 1 0.5]);

ax=subplot(1,3,1);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P5P8 & ~FS & ~SST,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P5P8 & ~FS & ~SST,:),1),'r','LineWidth',2)
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisualOpto(P5P8 & ~FS & ~SST,:),1))-1 max(mean(Z_PSTHvisualOpto(P5P8 & ~FS & ~SST,:),1))+1];

ax=subplot(1,3,2);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P9P13 & ~FS & ~SST,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P9P13 & ~FS & ~SST,:),1),'r','LineWidth',2)
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisual(P9P13 & ~FS & ~SST,:),1))-1 max(mean(Z_PSTHvisual(P9P13 & ~FS & ~SST,:),1))+1];

ax=subplot(1,3,3);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(Z_PSTHvisual(P14P18 & ~FS & ~SST,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisualOpto(P14P18 & ~FS & ~SST,:),1),'r','LineWidth',2)
ax=applyFont(ax,1);
ax.YLim=[min(mean(Z_PSTHvisual(P14P18 & ~FS & ~SST,:),1))-1 max(mean(Z_PSTHvisual(P14P18 & ~FS & ~SST,:),1))+1];
legend('Control','SST optogenetics')

% legend('Control','+ SST opto - laser only')
figname='SU_RS_WithOpto';
% sg=sgtitle(strcat('P',int2str(ops.age)));
% sg.FontSize=30;
export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')


figure('units','normalized','outerposition',[0 0 1 0.5]);

ax=subplot(1,3,1);
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(PSTHlaser(P5P8 & FS & ~SST,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[min(mean(PSTHlaser(P5P8 & FS & ~SST,:),1))-1 max(mean(PSTHlaser(P5P8 & FS & ~SST,:),1))+1];

ax=subplot(1,3,2);
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(PSTHlaser(P9P13 & FS & ~SST,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[min(mean(PSTHlaser(P9P13 & FS & ~SST,:),1))-1 max(mean(PSTHlaser(P9P13 & FS & ~SST,:),1))+1];

ax=subplot(1,3,3);
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(PSTHbins,mean(PSTHlaser(P14P18 & FS & ~SST,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[min(mean(PSTHlaser(P14P18 & FS & ~SST,:),1))-1 max(mean(PSTHlaser(P14P18 & FS & ~SST,:),1))+1];

% legend('Control','+ SST opto - laser only')
figname='SU_FS_Laser only';
% sg=sgtitle(strcat('P',int2str(ops.age)));
% sg.FontSize=30;
export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')

%% Bullshit analysis for Simon
FSnew=(FS&~SST);
RS=((~FSnew)&(~SST));
FS=FSnew;
% nnz(RS)+nnz(FSnew)+nnz(SST)


% Percent responsive per class per development window
nnz(responsive&RS&P5P8)/nnz(RS&P5P8)*100
nnz(responsive&RS&P9P13)/nnz(RS&P9P13)*100
nnz(responsive&RS&P14P18)/nnz(RS&P14P18)*100

nnz(responsive&FS&P5P8)/nnz(FS&P5P8)*100
nnz(responsive&FS&P9P13)/nnz(FS&P9P13)*100
nnz(responsive&FS&P14P18)/nnz(FS&P14P18)*100

nnz(responsive&SST&P5P8)/nnz(SST&P5P8)*100
nnz(responsive&SST&P9P13)/nnz(SST&P9P13)*100
nnz(responsive&SST&P14P18)/nnz(SST&P14P18)*100


%Calculate peak firing and peak latency
[PeakFiring,LatencyFiring]=max(PSTHvisual(:,PSTHbins>0&PSTHbins<300),[],2);
[~,minIdx]=min(abs(PSTHbins));
LatencyFiring=PSTHbins(LatencyFiring+minIdx)';
bullshit=table;

%Firing
%RS
bullshit.RSfiring_mean(1)=mean(PeakFiring(responsive&RS&P9P13));
bullshit.RSfiring_sem(1)=std(PeakFiring(responsive&RS&P9P13))/sqrt(numel(PeakFiring(responsive&RS&P9P13)));
bullshit.RSfiring_N(1)=numel(PeakFiring(responsive&RS&P9P13));

bullshit.RSfiring_mean(2)=mean(PeakFiring(responsive&RS&P14P18));
bullshit.RSfiring_sem(2)=std(PeakFiring(responsive&RS&P14P18))/sqrt(numel(PeakFiring(responsive&RS&P14P18)));
bullshit.RSfiring_N(2)=numel(PeakFiring(responsive&RS&P14P18));

%FS
bullshit.FSfiring_mean(1)=mean(PeakFiring(responsive&FS&P9P13));
bullshit.FSfiring_sem(1)=std(PeakFiring(responsive&FS&P9P13))/sqrt(numel(PeakFiring(responsive&FS&P9P13)));
bullshit.FSfiring_N(1)=numel(PeakFiring(responsive&FS&P9P13));

bullshit.FSfiring_mean(2)=mean(PeakFiring(responsive&FS&P14P18));
bullshit.FSfiring_sem(2)=std(PeakFiring(responsive&FS&P14P18))/sqrt(numel(PeakFiring(responsive&RS&P14P18)));
bullshit.FSfiring_N(2)=numel(PeakFiring(responsive&FS&P14P18));

%SST
bullshit.SSTfiring_mean(1)=mean(PeakFiring(responsive&SST&P9P13));
bullshit.SSTfiring_sem(1)=std(PeakFiring(responsive&SST&P9P13))/sqrt(numel(PeakFiring(responsive&SST&P9P13)));
bullshit.SSTfiring_N(1)=numel(PeakFiring(responsive&SST&P9P13));

bullshit.SSTfiring_mean(2)=mean(PeakFiring(responsive&SST&P14P18));
bullshit.SSTfiring_sem(2)=std(PeakFiring(responsive&SST&P14P18))/sqrt(numel(PeakFiring(responsive&SST&P14P18)));
bullshit.SSTfiring_N(2)=numel(PeakFiring(responsive&SST&P14P18));


%Latency
%RS
bullshit.RSlatency_mean(1)=mean(LatencyFiring(responsive&RS&P9P13));
bullshit.RSlatency_sem(1)=std(LatencyFiring(responsive&RS&P9P13))/sqrt(numel(LatencyFiring(responsive&RS&P9P13)));
bullshit.RSlatency_N(1)=numel(LatencyFiring(responsive&RS&P9P13));

bullshit.RSlatency_mean(2)=mean(LatencyFiring(responsive&RS&P14P18));
bullshit.RSlatency_sem(2)=std(LatencyFiring(responsive&RS&P14P18))/sqrt(numel(LatencyFiring(responsive&RS&P14P18)));
bullshit.RSlatency_N(2)=numel(LatencyFiring(responsive&RS&P14P18));

%FS
bullshit.FSlatency_mean(1)=mean(LatencyFiring(responsive&FS&P9P13));
bullshit.FSlatency_sem(1)=std(LatencyFiring(responsive&FS&P9P13))/sqrt(numel(LatencyFiring(responsive&FS&P9P13)));
bullshit.FSlatency_N(1)=numel(LatencyFiring(responsive&FS&P9P13));

bullshit.FSlatency_mean(2)=mean(LatencyFiring(responsive&FS&P14P18));
bullshit.FSlatency_sem(2)=std(LatencyFiring(responsive&FS&P14P18))/sqrt(numel(LatencyFiring(responsive&RS&P14P18)));
bullshit.FSlatency_N(2)=numel(LatencyFiring(responsive&FS&P14P18));

%SST
bullshit.SSTlatency_mean(1)=mean(LatencyFiring(responsive&SST&P9P13));
bullshit.SSTlatency_sem(1)=std(LatencyFiring(responsive&SST&P9P13))/sqrt(numel(LatencyFiring(responsive&SST&P9P13)));
bullshit.SSTlatency_N(1)=numel(LatencyFiring(responsive&SST&P9P13));

bullshit.SSTlatency_mean(2)=mean(LatencyFiring(responsive&SST&P14P18));
bullshit.SSTlatency_sem(2)=std(LatencyFiring(responsive&SST&P14P18))/sqrt(numel(LatencyFiring(responsive&SST&P14P18)));
bullshit.SSTlatency_N(2)=numel(LatencyFiring(responsive&SST&P14P18));


%Compare inhibited by SST vs non inhibited
% Percent inhibited per class per development window
nnz(SST_I&RS&P5P8)/nnz(RS&P5P8)*100
nnz(SST_I&RS&P9P13)/nnz(RS&P9P13)*100
nnz(SST_I&RS&P14P18)/nnz(RS&P14P18)*100

nnz(SST_I&FS&P5P8)/nnz(FS&P5P8)*100
nnz(SST_I&FS&P9P13)/nnz(FS&P9P13)*100
nnz(SST_I&FS&P14P18)/nnz(FS&P14P18)*100

bullshitInhibited=table;
bullshitNonInhibited=table;

%Firing
%RS inhibited
bullshitInhibited.RSfiring_mean(1)=mean(PeakFiring(responsive&RS&P9P13&SST_I));
bullshitInhibited.RSfiring_sem(1)=std(PeakFiring(responsive&RS&P9P13&SST_I))/sqrt(numel(PeakFiring(responsive&RS&P9P13&SST_I)));
bullshitInhibited.RSfiring_N(1)=numel(PeakFiring(responsive&RS&P9P13&SST_I));

bullshitInhibited.RSfiring_mean(2)=mean(PeakFiring(responsive&RS&P14P18&SST_I));
bullshitInhibited.RSfiring_sem(2)=std(PeakFiring(responsive&RS&P14P18&SST_I))/sqrt(numel(PeakFiring(responsive&RS&P14P18&SST_I)));
bullshitInhibited.RSfiring_N(2)=numel(PeakFiring(responsive&RS&P14P18&SST_I));

%FSinhibited
bullshitInhibited.FSfiring_mean(1)=mean(PeakFiring(responsive&FS&P9P13&SST_I));
bullshitInhibited.FSfiring_sem(1)=std(PeakFiring(responsive&FS&P9P13&SST_I))/sqrt(numel(PeakFiring(responsive&FS&P9P13&SST_I)));
bullshitInhibited.FSfiring_N(1)=numel(PeakFiring(responsive&FS&P9P13&SST_I));

bullshitInhibited.FSfiring_mean(2)=mean(PeakFiring(responsive&FS&P14P18&SST_I));
bullshitInhibited.FSfiring_sem(2)=std(PeakFiring(responsive&FS&P14P18&SST_I))/sqrt(numel(PeakFiring(responsive&RS&P14P18&SST_I)));
bullshitInhibited.FSfiring_N(2)=numel(PeakFiring(responsive&FS&P14P18&SST_I));

%RS non inhibited
bullshitNonInhibited.RSfiring_mean(1)=mean(PeakFiring(responsive&RS&P9P13&SST_Nonresponsive));
bullshitNonInhibited.RSfiring_sem(1)=std(PeakFiring(responsive&RS&P9P13&SST_Nonresponsive))/sqrt(numel(PeakFiring(responsive&RS&P9P13&SST_Nonresponsive)));
bullshitNonInhibited.RSfiring_N(1)=numel(PeakFiring(responsive&RS&P9P13&SST_Nonresponsive));

bullshitNonInhibited.RSfiring_mean(2)=mean(PeakFiring(responsive&RS&P14P18&SST_Nonresponsive));
bullshitNonInhibited.RSfiring_sem(2)=std(PeakFiring(responsive&RS&P14P18&SST_Nonresponsive))/sqrt(numel(PeakFiring(responsive&RS&P14P18&SST_Nonresponsive)));
bullshitNonInhibited.RSfiring_N(2)=numel(PeakFiring(responsive&RS&P14P18&SST_Nonresponsive));

%FS non inhibited
bullshitNonInhibited.FSfiring_mean(1)=mean(PeakFiring(responsive&FS&P9P13&SST_Nonresponsive));
bullshitNonInhibited.FSfiring_sem(1)=std(PeakFiring(responsive&FS&P9P13&SST_Nonresponsive))/sqrt(numel(PeakFiring(responsive&FS&P9P13&SST_Nonresponsive)));
bullshitNonInhibited.FSfiring_N(1)=numel(PeakFiring(responsive&FS&P9P13&SST_Nonresponsive));

bullshitNonInhibited.FSfiring_mean(2)=mean(PeakFiring(responsive&FS&P14P18&SST_Nonresponsive));
bullshitNonInhibited.FSfiring_sem(2)=std(PeakFiring(responsive&FS&P14P18&SST_Nonresponsive))/sqrt(numel(PeakFiring(responsive&RS&P14P18&SST_Nonresponsive)));
bullshitNonInhibited.FSfiring_N(2)=numel(PeakFiring(responsive&FS&P14P18&SST_Nonresponsive));


%Latency
%RS
bullshitInhibited.RSlatency_mean(1)=mean(LatencyFiring(responsive&RS&P9P13&SST_I));
bullshitInhibited.RSlatency_sem(1)=std(LatencyFiring(responsive&RS&P9P13&SST_I))/sqrt(numel(LatencyFiring(responsive&RS&P9P13&SST_I)));
bullshitInhibited.RSlatency_N(1)=numel(LatencyFiring(responsive&RS&P9P13&SST_I));

bullshitInhibited.RSlatency_mean(2)=mean(LatencyFiring(responsive&RS&P14P18&SST_I));
bullshitInhibited.RSlatency_sem(2)=std(LatencyFiring(responsive&RS&P14P18&SST_I))/sqrt(numel(LatencyFiring(responsive&RS&P14P18&SST_I)));
bullshitInhibited.RSlatency_N(2)=numel(LatencyFiring(responsive&RS&P14P18&SST_I));

%FS
bullshitInhibited.FSlatency_mean(1)=mean(LatencyFiring(responsive&FS&P9P13&SST_I));
bullshitInhibited.FSlatency_sem(1)=std(LatencyFiring(responsive&FS&P9P13&SST_I))/sqrt(numel(LatencyFiring(responsive&FS&P9P13&SST_I)));
bullshitInhibited.FSlatency_N(1)=numel(LatencyFiring(responsive&FS&P9P13&SST_I));

bullshitInhibited.FSlatency_mean(2)=mean(LatencyFiring(responsive&FS&P14P18&SST_I));
bullshitInhibited.FSlatency_sem(2)=std(LatencyFiring(responsive&FS&P14P18&SST_I))/sqrt(numel(LatencyFiring(responsive&RS&P14P18&SST_I)));
bullshitInhibited.FSlatency_N(2)=numel(LatencyFiring(responsive&FS&P14P18&SST_I));

%RS non inhibited
bullshitNonInhibited.RSlatency_mean(1)=mean(LatencyFiring(responsive&RS&P9P13&SST_Nonresponsive));
bullshitNonInhibited.RSlatency_sem(1)=std(LatencyFiring(responsive&RS&P9P13&SST_Nonresponsive))/sqrt(numel(LatencyFiring(responsive&RS&P9P13&SST_Nonresponsive)));
bullshitNonInhibited.RSlatency_N(1)=numel(LatencyFiring(responsive&RS&P9P13&SST_Nonresponsive));

bullshitNonInhibited.RSlatency_mean(2)=mean(LatencyFiring(responsive&RS&P14P18&SST_Nonresponsive));
bullshitNonInhibited.RSlatency_sem(2)=std(LatencyFiring(responsive&RS&P14P18&SST_Nonresponsive))/sqrt(numel(LatencyFiring(responsive&RS&P14P18&SST_Nonresponsive)));
bullshitNonInhibited.RSlatency_N(2)=numel(LatencyFiring(responsive&RS&P14P18&SST_Nonresponsive));

%FS non inhibited
bullshitNonInhibited.FSlatency_mean(1)=mean(LatencyFiring(responsive&FS&P9P13&SST_Nonresponsive));
bullshitNonInhibited.FSlatency_sem(1)=std(LatencyFiring(responsive&FS&P9P13&SST_Nonresponsive))/sqrt(numel(LatencyFiring(responsive&FS&P9P13&SST_Nonresponsive)));
bullshitNonInhibited.FSlatency_N(1)=numel(LatencyFiring(responsive&FS&P9P13&SST_Nonresponsive));

bullshitNonInhibited.FSlatency_mean(2)=mean(LatencyFiring(responsive&FS&P14P18&SST_Nonresponsive));
bullshitNonInhibited.FSlatency_sem(2)=std(LatencyFiring(responsive&FS&P14P18&SST_Nonresponsive))/sqrt(numel(LatencyFiring(responsive&RS&P14P18&SST_Nonresponsive)));
bullshitNonInhibited.FSlatency_N(2)=numel(LatencyFiring(responsive&FS&P14P18&SST_Nonresponsive));


%PSTH plotting
figure('units','normalized','outerposition',[0 0 1 0.5]);

ax=subplot(1,3,1);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P9P13,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P9P13,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P9P13,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(PSTHvisual(responsive&FS&P9P13,:),1))+1];
ax.XLim=[-200,1000];

ax=subplot(1,3,2);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P14P18,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P14P18,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P14P18,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(PSTHvisual(responsive&FS&P14P18,:),1))+1];
ax.XLim=[-200,1000];

legend('RS','FS','SST')
figname='SU_RSvsFSvsSST';
export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')

%% FENS poster
figure
hold on;
% patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,smooth(mean(Z_PSTHvisual(SST&P5P8,:),1),5),'r','LineWidth',2)
plot(PSTHbins,smooth(mean(Z_PSTHvisual(SST&P9P13,:),1),5),'y','LineWidth',2)
plot(PSTHbins,smooth(mean(Z_PSTHvisual(SST&P14P18,:),1),5),'b','LineWidth',2)
ax=gca;
ax=applyFont(ax,1);
ax.YLim=[-1 9];
ax.XLim=[-50,3000];

L4=(sulayer==2);
L23=(sulayer==1);
L56=(sulayer==3);

figure
hold on
plot(PSTHbins,smooth(mean(PSTHvisual(SST&P9P13&L23,:),1),5),'r','LineWidth',2)
plot(PSTHbins,smooth(mean(PSTHvisual(SST&P9P13&L4,:),1),5),'y','LineWidth',2)
plot(PSTHbins,smooth(mean(PSTHvisual(SST&P9P13&L56,:),1),5),'b','LineWidth',2)

figure
hold on
plot(PSTHbins,smooth(mean(PSTHvisual(SST&P14P18&L23,:),1),5),'r','LineWidth',2)
plot(PSTHbins,smooth(mean(PSTHvisual(SST&P14P18&L4,:),1),5),'y','LineWidth',2)
plot(PSTHbins,smooth(mean(PSTHvisual(SST&P14P18&L56,:),1),5),'b','LineWidth',2)

