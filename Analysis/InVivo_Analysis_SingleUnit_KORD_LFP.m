close all
clear 
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 
load('LFPData.mat')
folderFigures='C:\Users\Butt Lab\Desktop\Figures_InVivo';

%% Set group logic arrays
goodRecs=categorical(cellstr({'K6', 'K24', 'K35', 'K36', 'K38', 'K39', 'K43', 'K44', 'K45', 'K47', 'K48', 'K51', 'K52', 'K54', 'K55', 'K56', 'K58'}));
controlRecs=categorical(cellstr({'K20','K21','K23','K37','K41','K42','K46','K50','K53'}));

data=data(data.Tagging=='SST;KORD',:);
% data=data(ismember(data.MouseID,goodRecs),:);

%Age
P5P8=data.mouseAge<9;
P9P13=data.mouseAge>=9 & data.mouseAge<14;
P14P18=data.mouseAge>=14;

%Brain area
V1=data.brainArea=='V1';
S1=data.brainArea=='S1BF';
sizeS1young=height(data(S1&P9P13,:));
sizeS1old=height(data(S1&P14P18,:));
sizeV1young=height(data(V1&P9P13,:));
sizeV1old=height(data(V1&P14P18,:));

%
KORD=ismember(data.MouseID,goodRecs);
control=ismember(data.MouseID,controlRecs);

%% V1 responses
figure
ax1=subplot(2,2,1);
hold on
plot([data(V1&P9P13&KORD,:).MUA_L4_peakFast,data(V1&P9P13&KORD,:).MUA_L4_peakFast_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(V1&P9P13&KORD,:).MUA_L4_peakFast),mean(data(V1&P9P13&KORD,:).MUA_L4_peakFast_K)]',[std(data(V1&P9P13&KORD,:).MUA_L4_peakFast),std(data(V1&P9P13&KORD,:).MUA_L4_peakFast_K)]','k-o','LineWidth',2)
plot([data(V1&P9P13&control,:).MUA_L4_peakFast,data(V1&P9P13&control,:).MUA_L4_peakFast_K]','-o','Color',[254,227,145]/255)
errorbar([1,2]+.1,[mean(data(V1&P9P13&control,:).MUA_L4_peakFast),mean(data(V1&P9P13&control,:).MUA_L4_peakFast_K)]',[std(data(V1&P9P13&control,:).MUA_L4_peakFast),std(data(V1&P9P13&control,:).MUA_L4_peakFast_K)]','-o','LineWidth',2,'Color',[227,26,28]/255)
ax1.XLim=[.5,2.5];
ax1.Title.String='V1 P9-P13 - Primary response';
ax1.YLabel.String='MUA L4 peak firing (Hz)';

ax2=subplot(2,2,2);
hold on
plot([data(V1&P14P18&KORD,:).MUA_L4_peakFast,data(V1&P14P18&KORD,:).MUA_L4_peakFast_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(V1&P14P18&KORD,:).MUA_L4_peakFast),mean(data(V1&P14P18&KORD,:).MUA_L4_peakFast_K)]',[std(data(V1&P14P18&KORD,:).MUA_L4_peakFast),std(data(V1&P14P18&KORD,:).MUA_L4_peakFast_K)]','k-o','LineWidth',2)
plot([data(V1&P14P18&control,:).MUA_L4_peakFast,data(V1&P14P18&control,:).MUA_L4_peakFast_K]','-o','Color',[254,227,145]/255)
errorbar([1,2]+.1,[mean(data(V1&P14P18&control,:).MUA_L4_peakFast),mean(data(V1&P14P18&control,:).MUA_L4_peakFast_K)]',[std(data(V1&P14P18&control,:).MUA_L4_peakFast),std(data(V1&P14P18&control,:).MUA_L4_peakFast_K)]','-o','LineWidth',2,'Color',[227,26,28]/255)
ax2.XLim=[.5,2.5];
ax2.Title.String='V1 P14-P18 - Primary response';
ax1.YLabel.String='MUA L4 peak firing (Hz)';

ax1=subplot(2,2,3);
hold on
plot([data(V1&P9P13&KORD,:).MUA_L4_peakSlow,data(V1&P9P13&KORD,:).MUA_L4_peakSlow_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(V1&P9P13&KORD,:).MUA_L4_peakSlow),mean(data(V1&P9P13&KORD,:).MUA_L4_peakSlow_K)]',[std(data(V1&P9P13&KORD,:).MUA_L4_peakSlow),std(data(V1&P9P13&KORD,:).MUA_L4_peakSlow_K)]','k-o','LineWidth',2)
plot([data(V1&P9P13&control,:).MUA_L4_peakSlow,data(V1&P9P13&control,:).MUA_L4_peakFast_K]','-o','Color',[254,227,145]/255)
errorbar([1,2]+.1,[mean(data(V1&P9P13&control,:).MUA_L4_peakSlow),mean(data(V1&P9P13&control,:).MUA_L4_peakSlow_K)]',[std(data(V1&P9P13&control,:).MUA_L4_peakSlow),std(data(V1&P9P13&control,:).MUA_L4_peakSlow_K)]','-o','LineWidth',2,'Color',[227,26,28]/255)
ax1.XLim=[.5,2.5];
ax1.Title.String='V1 P9-P13 - Secondary response';
ax1.YLabel.String='MUA L4 peak firing (Hz)';

ax2=subplot(2,2,4);
hold on
plot([data(V1&P14P18&KORD,:).MUA_L4_peakSlow,data(V1&P14P18&KORD,:).MUA_L4_peakSlow_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(V1&P14P18&KORD,:).MUA_L4_peakSlow),mean(data(V1&P14P18&KORD,:).MUA_L4_peakSlow_K)]',[std(data(V1&P14P18&KORD,:).MUA_L4_peakSlow),std(data(V1&P14P18&KORD,:).MUA_L4_peakSlow_K)]','k-o','LineWidth',2)
plot([data(V1&P14P18&control,:).MUA_L4_peakSlow,data(V1&P14P18&control,:).MUA_L4_peakSlow_K]','-o','Color',[254,227,145]/255)
errorbar([1,2]+.1,[mean(data(V1&P14P18&control,:).MUA_L4_peakSlow),mean(data(V1&P14P18&control,:).MUA_L4_peakSlow_K)]',[std(data(V1&P14P18&control,:).MUA_L4_peakSlow),std(data(V1&P14P18&control,:).MUA_L4_peakSlow_K)]','-o','LineWidth',2,'Color',[227,26,28]/255)
ax2.XLim=[.5,2.5];
ax2.Title.String='V1 P14-P18 - Secondary response';
ax1.YLabel.String='MUA L4 peak firing (Hz)';

figure
cats={'V1 P9-P13 Fast','V1 P14-P18 Fast','V1 P9-P13 Slow','V1 P14-P18 Slow'};
hold on
plot(ones(height(data(V1&P9P13&KORD,:)),1),log2(data(V1&P9P13&KORD,:).MUA_L4_peakFast_K./(data(V1&P9P13&KORD,:).MUA_L4_peakFast+0.001)),'o','Color',[189,189,189]/255)
plot(ones(height(data(V1&P14P18&KORD,:)),1)*2,log2(data(V1&P14P18&KORD,:).MUA_L4_peakFast_K./(data(V1&P14P18&KORD,:).MUA_L4_peakFast+0.001)),'o','Color',[189,189,189]/255)
plot(ones(height(data(V1&P9P13&KORD,:)),1)*3,log2(data(V1&P9P13&KORD,:).MUA_L4_peakSlow_K./(data(V1&P9P13&KORD,:).MUA_L4_peakSlow+0.001)),'o','Color',[189,189,189]/255)
plot(ones(height(data(V1&P14P18&KORD,:)),1)*4,log2(data(V1&P14P18&KORD,:).MUA_L4_peakSlow_K./(data(V1&P14P18&KORD,:).MUA_L4_peakSlow+0.001)),'o','Color',[189,189,189]/255)
% errorbar([1,2,3,4]+0.1,[mean(data(V1&P9P13,:).baseline_firing_norm),mean(data(V1&P14P18,:).baseline_firing_norm),mean(data(S1&P9P13,:).baseline_firing_norm),mean(data(S1&P14P18,:).baseline_firing_norm)]',[std(data(V1&P9P13,:).baseline_firing_norm),std(data(V1&P14P18,:).baseline_firing_norm),std(data(S1&P9P13,:).baseline_firing_norm),std(data(S1&P14P18,:).baseline_firing_norm)]','ko')
plot(ones(height(data(V1&P9P13&control,:)),1),log2(data(V1&P9P13&control,:).MUA_L4_peakFast_K./(data(V1&P9P13&control,:).MUA_L4_peakFast+0.001)),'o','Color',[227,26,28]/255)
plot(ones(height(data(V1&P14P18&control,:)),1)*2,log2(data(V1&P14P18&control,:).MUA_L4_peakFast_K./(data(V1&P14P18&control,:).MUA_L4_peakFast+0.001)),'o','Color',[227,26,28]/255)
plot(ones(height(data(V1&P9P13&control,:)),1)*3,log2(data(V1&P9P13&control,:).MUA_L4_peakSlow_K./(data(V1&P9P13&control,:).MUA_L4_peakSlow+0.001)),'o','Color',[227,26,28]/255)
plot(ones(height(data(V1&P14P18&control,:)),1)*4,log2(data(V1&P14P18&control,:).MUA_L4_peakSlow_K./(data(V1&P14P18&control,:).MUA_L4_peakSlow+0.001)),'o','Color',[227,26,28]/255)

ax=gca;
ax.XLim=[.5,4.5];
ax.XTick=[1,2,3,4];
ax.XTickLabel=cats;
ax.YLabel.String='Mua L4 peak firing Change (log_2 ratio)';

%% S1 
figure
ax1=subplot(2,2,1);
hold on
plot([data(S1&P9P13&KORD,:).MUA_L4_peakFast,data(S1&P9P13&KORD,:).MUA_L4_peakFast_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(S1&P9P13&KORD,:).MUA_L4_peakFast),mean(data(S1&P9P13&KORD,:).MUA_L4_peakFast_K)]',[std(data(S1&P9P13&KORD,:).MUA_L4_peakFast),std(data(S1&P9P13&KORD,:).MUA_L4_peakFast_K)]','k-o','LineWidth',2)
plot([data(S1&P9P13&control,:).MUA_L4_peakFast,data(S1&P9P13&control,:).MUA_L4_peakFast_K]','-o','Color',[254,227,145]/255)
errorbar([1,2]+.1,[mean(data(S1&P9P13&control,:).MUA_L4_peakFast),mean(data(S1&P9P13&control,:).MUA_L4_peakFast_K)]',[std(data(S1&P9P13&control,:).MUA_L4_peakFast),std(data(S1&P9P13&control,:).MUA_L4_peakFast_K)]','-o','LineWidth',2,'Color',[227,26,28]/255)
ax1.XLim=[.5,2.5];
ax1.Title.String='V1 P9-P13 - Primary response';
ax1.YLabel.String='MUA L4 peak firing (Hz)';

ax2=subplot(2,2,2);
hold on
plot([data(S1&P14P18&KORD,:).MUA_L4_peakFast,data(S1&P14P18&KORD,:).MUA_L4_peakFast_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(S1&P14P18&KORD,:).MUA_L4_peakFast),mean(data(S1&P14P18&KORD,:).MUA_L4_peakFast_K)]',[std(data(S1&P14P18&KORD,:).MUA_L4_peakFast),std(data(S1&P14P18&KORD,:).MUA_L4_peakFast_K)]','k-o','LineWidth',2)
plot([data(S1&P14P18&control,:).MUA_L4_peakFast,data(S1&P14P18&control,:).MUA_L4_peakFast_K]','-o','Color',[254,227,145]/255)
errorbar([1,2]+.1,[mean(data(S1&P14P18&control,:).MUA_L4_peakFast),mean(data(S1&P14P18&control,:).MUA_L4_peakFast_K)]',[std(data(S1&P14P18&control,:).MUA_L4_peakFast),std(data(S1&P14P18&control,:).MUA_L4_peakFast_K)]','-o','LineWidth',2,'Color',[227,26,28]/255)
ax2.XLim=[.5,2.5];
ax2.Title.String='V1 P14-P18 - Primary response';
ax1.YLabel.String='MUA L4 peak firing (Hz)';

ax1=subplot(2,2,3);
hold on
plot([data(S1&P9P13&KORD,:).MUA_L4_peakSlow,data(S1&P9P13&KORD,:).MUA_L4_peakSlow_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(S1&P9P13&KORD,:).MUA_L4_peakSlow),mean(data(S1&P9P13&KORD,:).MUA_L4_peakSlow_K)]',[std(data(S1&P9P13&KORD,:).MUA_L4_peakSlow),std(data(S1&P9P13&KORD,:).MUA_L4_peakSlow_K)]','k-o','LineWidth',2)
plot([data(S1&P9P13&control,:).MUA_L4_peakSlow,data(S1&P9P13&control,:).MUA_L4_peakFast_K]','-o','Color',[254,227,145]/255)
errorbar([1,2]+.1,[mean(data(S1&P9P13&control,:).MUA_L4_peakSlow),mean(data(S1&P9P13&control,:).MUA_L4_peakSlow_K)]',[std(data(S1&P9P13&control,:).MUA_L4_peakSlow),std(data(S1&P9P13&control,:).MUA_L4_peakSlow_K)]','-o','LineWidth',2,'Color',[227,26,28]/255)
ax1.XLim=[.5,2.5];
ax1.Title.String='V1 P9-P13 - Secondary response';
ax1.YLabel.String='MUA L4 peak firing (Hz)';

ax2=subplot(2,2,4);
hold on
plot([data(S1&P14P18&KORD,:).MUA_L4_peakSlow,data(S1&P14P18&KORD,:).MUA_L4_peakSlow_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(S1&P14P18&KORD,:).MUA_L4_peakSlow),mean(data(S1&P14P18&KORD,:).MUA_L4_peakSlow_K)]',[std(data(S1&P14P18&KORD,:).MUA_L4_peakSlow),std(data(S1&P14P18&KORD,:).MUA_L4_peakSlow_K)]','k-o','LineWidth',2)
plot([data(S1&P14P18&control,:).MUA_L4_peakSlow,data(S1&P14P18&control,:).MUA_L4_peakSlow_K]','-o','Color',[254,227,145]/255)
errorbar([1,2]+.1,[mean(data(S1&P14P18&control,:).MUA_L4_peakSlow),mean(data(S1&P14P18&control,:).MUA_L4_peakSlow_K)]',[std(data(S1&P14P18&control,:).MUA_L4_peakSlow),std(data(S1&P14P18&control,:).MUA_L4_peakSlow_K)]','-o','LineWidth',2,'Color',[227,26,28]/255)
ax2.XLim=[.5,2.5];
ax2.Title.String='V1 P14-P18 - Secondary response';
ax1.YLabel.String='MUA L4 peak firing (Hz)';

%%PPR
data.MUA_L4_PPR=data.MUA_L4_peakSlow./data.MUA_L4_peakFast;
data.MUA_L4_PPR_K=data.MUA_L4_peakSlow_K./data.MUA_L4_peakFast_K;

figure
ax1=subplot(2,2,1);
hold on
plot([data(S1&P9P13&KORD,:).MUA_L4_PPR,data(S1&P9P13&KORD,:).MUA_L4_PPR_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(S1&P9P13&KORD,:).MUA_L4_PPR),mean(data(S1&P9P13&KORD,:).MUA_L4_PPR_K)]',[std(data(S1&P9P13&KORD,:).MUA_L4_PPR),std(data(S1&P9P13&KORD,:).MUA_L4_PPR_K)]','k-o')
plot([data(S1&P9P13&control,:).MUA_L4_PPR,data(S1&P9P13&control,:).MUA_L4_PPR_K]','-o','Color',[254,227,145]/255)
errorbar([1,2]+.1,[mean(data(S1&P9P13&control,:).MUA_L4_PPR),mean(data(S1&P9P13&control,:).MUA_L4_PPR_K)]',[std(data(S1&P9P13&control,:).MUA_L4_PPR),std(data(S1&P9P13&control,:).MUA_L4_PPR_K)]','-o','Color',[227,26,28]/255)
ax1.XLim=[.5,2.5];
ax1.Title.String='S1 P9-P13';
ax1.YLabel.String='PPR';

ax2=subplot(2,2,2);
hold on
plot([data(S1&P14P18&KORD,:).MUA_L4_PPR,data(S1&P14P18&KORD,:).MUA_L4_PPR_K]','-o','Color',[189,189,189]/255)
errorbar([1,2]+.1,[mean(data(S1&P14P18&KORD,:).MUA_L4_PPR),mean(data(S1&P14P18&KORD,:).MUA_L4_PPR_K)]',[std(data(S1&P14P18&KORD,:).MUA_L4_PPR),std(data(S1&P14P18&KORD,:).MUA_L4_PPR_K)]','k-o')
plot([data(S1&P14P18&control,:).MUA_L4_PPR,data(S1&P14P18&control,:).MUA_L4_PPR_K]','-o','Color',[254,227,145]/255)
errorbar([1,2]+.1,[mean(data(S1&P14P18&control,:).MUA_L4_PPR),mean(data(S1&P14P18&control,:).MUA_L4_PPR_K)]',[std(data(S1&P14P18&control,:).MUA_L4_PPR),std(data(S1&P14P18&control,:).MUA_L4_PPR_K)]','-o','Color',[227,26,28]/255)
ax2.XLim=[.5,2.5];
ax2.Title.String='S1 P14-P18';
ax2.YLabel.String='PPR';
ax2.YLabel.String='PPR';
