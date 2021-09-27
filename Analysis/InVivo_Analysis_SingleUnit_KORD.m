close all
clear 
clc

addpath(genpath('C:\Users\Filippo\Documents\GitHub\InVivoEphys_Analysis')) 
load('SingleUnitData.mat')
data.PSTHvisualOpto=[];
data.PSTHlaser=[];
data.responseVisualOpto=[];
data.responseLaser=[];
% data1=data;
% load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SingleUnitData_Part2.mat')
% data=[data1;data];
folderFigures='C:\Users\Filippo\Desktop\Figures_InVivo';

%% Set group logic arrays
goodRecs=categorical(cellstr({'K6', 'K24', 'K35', 'K36', 'K38', 'K39', 'K43', 'K44', 'K45', 'K47', 'K48', 'K51', 'K52', 'K54', 'K55', 'K56', 'K58'}));
data=data(~isnan(data.endSlope),:);
data=data(data.Tagging=='SST;KORD',:);
data=data(ismember(data.MouseID,goodRecs),:);


% data=data(data.brainArea=='V1',:);


%Layer
tmpLayers(data.Layer==1,1)=categorical(cellstr('L2/3'));
tmpLayers(data.Layer==2,1)=categorical(cellstr('L4'));
tmpLayers(data.Layer==3,1)=categorical(cellstr('L5/6'));
data.Layer=tmpLayers;

%Age
P5P8=data.Age<9;
P9P13=data.Age>=9 & data.Age<14;
P14P18=data.Age>=14;

%Brain area
V1=data.brainArea=='V1';
S1=data.brainArea=='S1BF';

%Visual response
% responsive=any(data.responseVisual(:,1:2)==1,2);
% responsive=any(data.responseVisual==1,2);
sizeS1young=height(data(S1&P9P13,:));
sizeS1old=height(data(S1&P14P18,:));
sizeV1young=height(data(V1&P9P13,:));
sizeV1old=height(data(V1&P14P18,:));

%% Baseline firing
cats={'V1 P9-P13','V1 P14-P18','S1 P9-P13','S1 P14-P18'};

figure
ax1=subplot(2,2,1);
hold on
plot([data(V1&P9P13,:).baseline_firing,data(V1&P9P13,:).baseline_firing_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(V1&P9P13,:).baseline_firing),mean(data(V1&P9P13,:).baseline_firing_K)]',[std(data(V1&P9P13,:).baseline_firing),std(data(V1&P9P13,:).baseline_firing_K)]','k-o')
ax1.XLim=[.5,2.5];
ax1.Title.String='V1 P9-P13';
ax1.YLabel.String='Baseline firing (Hz)';

ax2=subplot(2,2,2);
hold on
plot([data(V1&P14P18,:).baseline_firing,data(V1&P14P18,:).baseline_firing_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(V1&P14P18,:).baseline_firing),mean(data(V1&P14P18,:).baseline_firing_K)]',[std(data(V1&P14P18,:).baseline_firing),std(data(V1&P14P18,:).baseline_firing_K)]','k-o')
ax2.XLim=[.5,2.5];
ax2.Title.String='V1 P14-P18';
ax2.YLabel.String='Baseline firing (Hz)';

ax1=subplot(2,2,3);
hold on
plot([data(S1&P9P13,:).baseline_firing,data(S1&P9P13,:).baseline_firing_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(S1&P9P13,:).baseline_firing),mean(data(S1&P9P13,:).baseline_firing_K)]',[std(data(S1&P9P13,:).baseline_firing),std(data(S1&P9P13,:).baseline_firing_K)]','k-o')
ax1.XLim=[.5,2.5];
ax1.Title.String='V1 P9-P13';
ax1.YLabel.String='Baseline firing (Hz)';

ax2=subplot(2,2,4);
hold on
plot([data(S1&P14P18,:).baseline_firing,data(S1&P14P18,:).baseline_firing_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(S1&P14P18,:).baseline_firing),mean(data(S1&P14P18,:).baseline_firing_K)]',[std(data(S1&P14P18,:).baseline_firing),std(data(S1&P14P18,:).baseline_firing_K)]','k-o')
ax2.XLim=[.5,2.5];
ax2.Title.String='V1 P14-P18';
ax2.YLabel.String='Baseline firing (Hz)';

figure
data.baseline_firing_norm=log2(data.baseline_firing_K./(data.baseline_firing+0.001));
hold on
plot(ones(sizeV1young,1),data(V1&P9P13,:).baseline_firing_norm,'o','Color',[189,189,189]/255)
plot(ones(sizeV1old,1)*2,data(V1&P14P18,:).baseline_firing_norm,'o','Color',[189,189,189]/255)
plot(ones(sizeS1young,1)*3,data(S1&P9P13,:).baseline_firing_norm,'o','Color',[189,189,189]/255)
plot(ones(sizeS1old,1)*4,data(S1&P14P18,:).baseline_firing_norm,'o','Color',[189,189,189]/255)
errorbar([1,2,3,4]+0.1,[mean(data(V1&P9P13,:).baseline_firing_norm),mean(data(V1&P14P18,:).baseline_firing_norm),mean(data(S1&P9P13,:).baseline_firing_norm),mean(data(S1&P14P18,:).baseline_firing_norm)]',[std(data(V1&P9P13,:).baseline_firing_norm),std(data(V1&P14P18,:).baseline_firing_norm),std(data(S1&P9P13,:).baseline_firing_norm),std(data(S1&P14P18,:).baseline_firing_norm)]','ko')
ax=gca;
ax.XLim=[.5,4.5];
ax.XTick=[1,2,3,4];
ax.XTickLabel=cats;
ax.YLabel.String='Baseline firing Change (log_2 ratio)';

%% RW firing
cats={'V1 P9-P13','V1 P14-P18','S1 P9-P13','S1 P14-P18'};

figure
ax1=subplot(2,2,1);
hold on
plot([data(V1&P9P13,:).rw_firing,data(V1&P9P13,:).rw_firing_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(V1&P9P13,:).rw_firing),mean(data(V1&P9P13,:).rw_firing_K)]',[std(data(V1&P9P13,:).rw_firing),std(data(V1&P9P13,:).rw_firing_K)]','k-o')
ax1.XLim=[.5,2.5];
ax1.Title.String='V1 P9-P13';
ax1.YLabel.String='RW firing (Hz)';

ax2=subplot(2,2,2);
hold on
plot([data(V1&P14P18,:).rw_firing,data(V1&P14P18,:).rw_firing_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(V1&P14P18,:).rw_firing),mean(data(V1&P14P18,:).rw_firing_K)]',[std(data(V1&P14P18,:).rw_firing),std(data(V1&P14P18,:).rw_firing_K)]','k-o')
ax2.XLim=[.5,2.5];
ax2.Title.String='V1 P14-P18';
ax2.YLabel.String='RW firing (Hz)';


figure
data.rw_firing_norm=log2(data.rw_firing_K./(data.rw_firing+0.001));
hold on
plot(ones(sizeV1young,1),data(V1&P9P13,:).rw_firing_norm,'o','Color',[189,189,189]/255)
plot(ones(sizeV1old,1)*2,data(V1&P14P18,:).rw_firing_norm,'o','Color',[189,189,189]/255)
errorbar([1,2]+0.1,[mean(data(V1&P9P13,:).rw_firing_norm),mean(data(V1&P14P18,:).rw_firing_norm)]',[std(data(V1&P9P13,:).rw_firing_norm),std(data(V1&P14P18,:).rw_firing_norm)]','ko')
ax=gca;
ax.XLim=[.5,2.5];
ax.XTick=[1,2];
ax.XTickLabel=cats(1:2);
ax.YLabel.String='RW firing Change (log_2 ratio)';


data.delta_rw_firing=data.rw_firing-data.rw_baselineFiring;
data.delta_rw_firing_K=data.rw_firing_K-data.rw_baselineFiring_K;
figure
ax1=subplot(2,2,1);
hold on
plot([data(V1&P9P13,:).delta_rw_firing,data(V1&P9P13,:).delta_rw_firing_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(V1&P9P13,:).delta_rw_firing),mean(data(V1&P9P13,:).delta_rw_firing_K)]',[std(data(V1&P9P13,:).delta_rw_firing),std(data(V1&P9P13,:).delta_rw_firing_K)]','k-o')
ax1.XLim=[.5,2.5];
ax1.Title.String='V1 P9-P13';
ax1.YLabel.String='Norm RW firing (Hz)';

ax2=subplot(2,2,2);
hold on
plot([data(V1&P14P18,:).delta_rw_firing,data(V1&P14P18,:).delta_rw_firing_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(V1&P14P18,:).delta_rw_firing),mean(data(V1&P14P18,:).delta_rw_firing_K)]',[std(data(V1&P14P18,:).delta_rw_firing),std(data(V1&P14P18,:).delta_rw_firing_K)]','k-o')
ax2.XLim=[.5,2.5];
ax2.Title.String='V1 P14-P18';
ax2.YLabel.String='Norm RW firing (Hz)';

figure
data.delta_rw_firing_norm=log2(data.delta_rw_firing_K./(data.delta_rw_firing+0.001));
hold on
plot(ones(sizeV1young,1),data(V1&P9P13,:).delta_rw_firing_norm,'o','Color',[189,189,189]/255)
plot(ones(sizeV1old,1)*2,data(V1&P14P18,:).delta_rw_firing_norm,'o','Color',[189,189,189]/255)
errorbar([1,2]+0.1,[mean(data(V1&P9P13,:).delta_rw_firing_norm),mean(data(V1&P14P18,:).delta_rw_firing_norm)]',[std(data(V1&P9P13,:).delta_rw_firing_norm),std(data(V1&P14P18,:).delta_rw_firing_norm)]','ko')
ax=gca;
ax.XLim=[.5,2.5];
ax.XTick=[1,2];
ax.XTickLabel=cats(1:2);
ax.YLabel.String='Norm RW firing Change (log_2 ratio)';

%% Fast response to stimulus
responsive_V1=any(data.responseVisual(:,1:2)==1,2);
responsive_S1=data.responseWhisker==1;
responsive_V1_K=any(data.responseVisual_K(:,1:2)==1,2);
responsive_S1_K=data.responseWhisker_K==1;

maxResponseV1_y=max(data.PSTHvisual(responsive_V1&P9P13,PSTHbins>0&PSTHbins<=200),[],2);
maxResponseV1_o=max(data.PSTHvisual(responsive_V1&P14P18,PSTHbins>0&PSTHbins<=200),[],2);
maxResponseV1_y_K=max(data.PSTHvisual_K(responsive_V1&P9P13,PSTHbins>0&PSTHbins<=200),[],2);
maxResponseV1_o_K=max(data.PSTHvisual_K(responsive_V1&P14P18,PSTHbins>0&PSTHbins<=200),[],2);

maxResponseS1_y=max(data.PSTHwhisker(responsive_S1&P9P13,PSTHbins>0&PSTHbins<=200),[],2);
maxResponseS1_o=max(data.PSTHwhisker(responsive_S1&P14P18,PSTHbins>0&PSTHbins<=200),[],2);
maxResponseS1_y_K=max(data.PSTHwhisker_K(responsive_S1&P9P13,PSTHbins>0&PSTHbins<=200),[],2);
maxResponseS1_o_K=max(data.PSTHwhisker_K(responsive_S1&P14P18,PSTHbins>0&PSTHbins<=200),[],2);

figure
ax1=subplot(2,2,1);
hold on
plot([maxResponseV1_y,maxResponseV1_y_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(maxResponseV1_y),mean(maxResponseV1_y_K)]',[std(maxResponseV1_y),std(maxResponseV1_y_K)]','k-o')
ax1.XLim=[.5,2.5];
ax1.Title.String='V1 P9-P13';
ax1.YLabel.String='PSTH peak firing (Hz)';

ax2=subplot(2,2,2);
hold on
plot([maxResponseV1_o,maxResponseV1_o_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(maxResponseV1_o),mean(maxResponseV1_o_K)]',[std(maxResponseV1_o),std(maxResponseV1_o_K)]','k-o')
ax2.XLim=[.5,2.5];
ax2.Title.String='V1 P14-P18';
ax1.YLabel.String='PSTH peak firing (Hz)';

ax1=subplot(2,2,3);
hold on
plot([maxResponseS1_y,maxResponseS1_y_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(maxResponseS1_y),mean(maxResponseS1_y_K)]',[std(maxResponseS1_y),std(maxResponseS1_y_K)]','k-o')
ax1.XLim=[.5,2.5];
ax1.Title.String='S1 P9-P13';
ax1.YLabel.String='PSTH peak firing (Hz)';

ax2=subplot(2,2,4);
hold on
plot([maxResponseS1_o,maxResponseS1_o_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(maxResponseS1_o),mean(maxResponseS1_o)]',[std(maxResponseS1_o),std(maxResponseS1_o)]','k-o')
ax2.XLim=[.5,2.5];
ax2.Title.String='S1 P14-P18';
ax1.YLabel.String='PSTH peak firing (Hz)';

figure
data.baseline_firing_norm=log2(data.baseline_firing_K./(data.baseline_firing+0.001));
hold on
plot(ones(numel(maxResponseV1_y),1),log2(maxResponseV1_y_K./maxResponseV1_y),'o','Color',[189,189,189]/255)
plot(ones(numel(maxResponseV1_o),1)*2,log2(maxResponseV1_o_K./maxResponseV1_o),'o','Color',[189,189,189]/255)
plot(ones(numel(maxResponseS1_y),1)*3,log2(maxResponseS1_y_K./maxResponseS1_y),'o','Color',[189,189,189]/255)
plot(ones(numel(maxResponseS1_o),1)*4,log2(maxResponseS1_o_K./maxResponseS1_o),'o','Color',[189,189,189]/255)
errorbar([1,2,3,4]+0.1,[mean(log2(maxResponseV1_y_K./maxResponseV1_y)),mean(log2(maxResponseV1_o_K./maxResponseV1_o)),mean(log2(maxResponseS1_y_K./maxResponseS1_y)),mean(log2(maxResponseS1_o_K./maxResponseS1_o))]',[std(log2(maxResponseV1_y_K./maxResponseV1_y)),std(log2(maxResponseV1_o_K./maxResponseV1_o)),std(log2(maxResponseS1_y_K./maxResponseS1_y)),std(log2(maxResponseS1_o_K./maxResponseS1_o))]','ko')
ax=gca;
ax.XLim=[.5,4.5];
ax.XTick=[1,2,3,4];
ax.XTickLabel=cats;
ax.YLabel.String='PSTH peak firing Change (log_2 ratio)';

%% PPR
data.PPR=max(data.PSTHwhisker(:,PSTHbins>=500 & PSTHbins<=550),[],2)./max(data.PSTHwhisker(:,PSTHbins>=0 & PSTHbins<=50),[],2);
data.PPR_K=max(data.PSTHwhisker_K(:,PSTHbins>=500 & PSTHbins<=550),[],2)./max(data.PSTHwhisker_K(:,PSTHbins>=0 & PSTHbins<=50),[],2);
data.PPR(isnan(data.PPR))=0;
data.PPR_K(isnan(data.PPR_K))=0;
data.PPR(isinf(data.PPR))=0;
data.PPR_K(isinf(data.PPR_K))=0;

figure
ax1=subplot(2,2,1);
hold on
plot([data(S1&P9P13&responsive_S1,:).PPR,data(S1&P9P13&responsive_S1,:).PPR_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(S1&P9P13&responsive_S1,:).PPR),mean(data(S1&P9P13&responsive_S1,:).PPR_K)]',[std(data(S1&P9P13&responsive_S1,:).PPR),std(data(S1&P9P13&responsive_S1,:).PPR_K)]','k-o')
ax1.XLim=[.5,2.5];
ax1.YLim=[-0,2];
ax1.Title.String='S1 P9-P13';
ax1.YLabel.String='PPR';

ax1=subplot(2,2,2);
hold on
plot([data(S1&P14P18&responsive_S1,:).PPR,data(S1&P14P18&responsive_S1,:).PPR_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(S1&P14P18&responsive_S1,:).PPR),mean(data(S1&P14P18&responsive_S1,:).PPR_K)]',[std(data(S1&P14P18&responsive_S1,:).PPR),std(data(S1&P14P18&responsive_S1,:).PPR_K)]','k-o')
ax1.XLim=[.5,2.5];
ax1.YLim=[-0,2];
ax1.Title.String='S1 P14-P18';
ax1.YLabel.String='PPR';


