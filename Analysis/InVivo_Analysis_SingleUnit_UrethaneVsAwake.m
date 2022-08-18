close all
clear 
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 
load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SingleUnitData.mat')
folderFigures='C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\WT Doctoral Programme in Neuroscience\DPhil in Neuroscience\DPhil thesis\Figures\Chapter 4';

%% Set group logic arrays
data=data(~isnan(data.endSlope),:);
data=data(data.Tagging~='SST;NrgKO',:);
data=data(data.brainArea=='V1',:);

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
P9P13=data.Age<14;
P14P18=data.Age>=14;
data.Dev(P5P8)=categorical(cellstr('P5-P8'));
data.Dev(P9P13)=categorical(cellstr('P9-P13'));
data.Dev(P14P18)=categorical(cellstr('P14-P18'));

%Cell Type
SST=(data.Tagging=='SST' & data.responseTag==1);
Nkx=(data.Tagging=='Nkx2-1' & data.responseTag==1);
Untagged=(~SST & ~Nkx);

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

aw = data.state=='Awake';
ur = data.state=='Urethane';


%% Preprocess data
% [Z_PSTHvisual, responsiveVisual]=zscoreBaseline(data.PSTHvisual);
% Z_PSTHvisualOpto=zscoreBaseline(data.PSTHvisualOpto);
% data.Z_PSTHlaser=zscore(data.PSTHlaser,[],'all');
% Z_PSTHoptotagging=zscoreBaseline(data.PSTHoptotagging);
Z_PSTHvisual=zscore(data.PSTHvisual,[],'all');
% responsive=any(Z_PSTHvisual(:,PSTHbins>0 & PSTHbins<200)>=3,2);
responsive=data.responseVisual==1;
resp1=responsive;
resp2=data.rw_responsive;

data.peakVisualFast = max(data.PSTHvisual(:,(PSTHbins>0 & PSTHbins<200)),[],2);
data.peakVisualFast_N = data.peakVisualFast-data.rw_baselineFiring;

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

%% Responsive bar plots
ybar=[nnz(RS&resp1&P9P13&ur)/nnz(RS&P9P13&ur),nnz(RS&resp1&P9P13&aw)/nnz(RS&aw&P9P13);...
    nnz(RS&resp2&P9P13&ur)/nnz(RS&P9P13&ur),nnz(RS&resp2&P9P13&aw)/nnz(RS&aw&P9P13)];
ybar2=[nnz(SB_entrained(RS&P9P13&ur))/nnz(RS&P9P13&ur),nnz(SB_entrained(RS&aw&P9P13))/nnz(RS&aw&P9P13);...
            nnz(PPC_entrained(RS&P9P13&ur))/nnz(RS&P9P13&ur),nnz(PPC_entrained(RS&aw&P9P13))/nnz(RS&aw&P9P13);];

figure('units','normalized','outerposition',[0 0 0.15 1]);
ax=subplot(3,1,1);
b = bar([1,2],ybar);
for i=1:size(b,2)
    b(1,i).LineWidth=1;
end
ax.YLabel.String='Ratio single units';
ax.Box='off';
ax.LineWidth = 1;
ax.FontSize=10;
ax.YLim=[0,1];
legend('Urethane','Awake')
legend('boxoff')

ax=subplot(3,1,2);
b = bar([1,2],ybar2);
for i=1:size(b,2)
    b(1,i).LineWidth=1;
end
ax.YLabel.String='Ratio single units';
ax.Box='off';
ax.LineWidth = 1;
ax.FontSize=10;
ax.YLim=[0,1];
legend('Urethane','Awake')
legend('boxoff')
export_fig(fullfile(folderFigures,'4.21','Barplots'),'-pdf','-transparent','-nocrop')
% close

%% RS units firing
clear ax v
figure('units','normalized','outerposition',[0 0 0.2 1])
ax(1)=subplot(3,2,1);
v(1,:)=violinplot(data.peakVisualFast_N(RS&P9P13&resp1),data.state(RS&P9P13&resp1));
ax(1).YLim=[0, 80];
ax(1).YAxis.Label.String='Max fast PSTH (spike/s)';
[h,p,stats] = my_ttest2 (data.peakVisualFast_N(RS&P9P13&resp1&ur),data.peakVisualFast_N(RS&P9P13&resp1&aw))

ax(2)=subplot(3,2,2);
v(2,:)=violinplot(data.rw_firing_N(RS&P9P13&resp2),data.state(RS&P9P13&resp2));
ax(2).YLim=[0, 10];
ax(2).YAxis.Label.String='Average slow PSTH (spike/s)';
[h,p,stats] = my_ttest2 (data.rw_firing_N(RS&P9P13&resp2&ur),data.rw_firing_N(RS&P9P13&resp2&aw))

ax(3)=subplot(3,2,3);
v(3,:)=violinplot(data.fanoVisual(RS&P9P13&resp1),data.state(RS&P9P13&resp1));
ax(3).YLim=[0, 3];
ax(3).YAxis.Label.String='Fano factor';
[h,p,stats] = my_ttest2 (data.fanoVisual(RS&P9P13&resp1&ur),data.fanoVisual(RS&P9P13&resp1&aw))

ax(4)=subplot(3,2,4);
v(4,:)=violinplot(data.baseline_firing(RS&P9P13),data.state(RS&P9P13));
ax(4).YLim=[0, 5];
ax(4).YAxis.Label.String='Spontaneous firing (spike/s)';
[h,p,stats] = my_ttest2 (data.baseline_firing(RS&P9P13&ur),data.baseline_firing(RS&P9P13&aw))

ax(5)=subplot(3,2,5);
v(5,:)=violinplot(data.sb_spikeProb(RS&P9P13&SB_entrained),data.state(RS&P9P13&SB_entrained));
ax(5).YLim=[0, 1];
ax(5).YAxis.Label.String='Probability of spike within SB';
[h,p,stats] = my_ttest2 (data.sb_spikeProb(RS&P9P13&ur&SB_entrained),data.sb_spikeProb(RS&P9P13&aw&SB_entrained))

ax(6)=subplot(3,2,6);
v(6,:)=violinplot(data.PPC(RS&P9P13&PPC_entrained,2),data.state(RS&P9P13&PPC_entrained));
ax(6).YLim=[0, 0.4];
ax(6).YAxis.Label.String='PPC';
[h,p,stats] = my_ttest2 (data.PPC(RS&P9P13&ur&PPC_entrained,2),data.PPC(RS&P9P13&aw&PPC_entrained,2))

color = [200,200,200;20,20,20];
for i=1:numel(ax)
    for j=1:size(v,2)
        v(i,j).ViolinColor=color(j,:)/255;
        v(i,j).ScatterPlot.MarkerFaceColor=[37,37,37]/255;
        v(i,j).ScatterPlot.MarkerFaceAlpha=1;
        v(i,j).ScatterPlot.SizeData=2;
    end

    ax(i).FontSize=12;
    ax(i).LineWidth=1;
    ax(i).XLim=[0.5,2.5];
    
end
print(gcf,'-dpdf',fullfile(folderFigures,'4.21','ViolinPlots'))
% close 