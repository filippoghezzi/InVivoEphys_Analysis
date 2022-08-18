close all
clear 
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 
addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\npy-matlab')) 

load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SingleUnitData.mat')
folderFigures='C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\WT Doctoral Programme in Neuroscience\DPhil in Neuroscience\DPhil thesis\Figures\Chapter 5';

%% Set group logic arrays
data=data(data.Tagging=='NTSR1',:);
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
P9P13=data.Age>=9 & data.Age<14;
P14P18=data.Age>=14;
data.Dev(P5P8)=categorical(cellstr('P5-P8'));
data.Dev(P9P13)=categorical(cellstr('P9-P13'));
data.Dev(P14P18)=categorical(cellstr('P14-P18'));

isLaser = ~(isnan(data.PSTHlaser(:,1)));

CT = data.responseTag==1;
CC = (data.Depth<max(data.Depth(CT))) & ~CT;
L5 = L56 & ~CT & ~CC;

FS = (data.troughPeakTime<0.75) & P14P18;
RS = ~FS;

data.CellIdentity(CT&RS)=categorical(cellstr('CT'));
data.CellIdentity(CC&RS)=categorical(cellstr('CC'));
data.CellIdentity(L5&RS)=categorical(cellstr('L5'));
data.CellIdentity(L23&RS)=categorical(cellstr('L23'));
data.CellIdentity(L4&RS)=categorical(cellstr('L4'));
data.CellIdentity(FS)=categorical(cellstr('FS'));

resp1=data.responseVisual==1;
resp2=data.rw_responsive;


%% Preprocess data
% [Z_PSTHvisual, responsiveVisual]=zscoreBaseline(data.PSTHvisual);
% Z_PSTHvisualOpto=zscoreBaseline(data.PSTHvisualOpto);
% data.Z_PSTHlaser=zscore(data.PSTHlaser,[],'all');
% Z_PSTHoptotagging=zscoreBaseline(data.PSTHoptotagging);
Z_PSTHvisual=zscore(data.PSTHvisual,[],'all');
% responsive=any(Z_PSTHvisual(:,PSTHbins>0 & PSTHbins<200)>=3,2);
responsive=data.responseVisual==1;


data.peakVisualFast = max(data.PSTHvisual(:,(PSTHbins>0 & PSTHbins<150)),[],2);
data.peakVisualFast_N = data.peakVisualFast-data.rw_baselineFiring;

data.peakVisualFast_K = max(data.PSTHvisual_K(:,(PSTHbins>0 & PSTHbins<150)),[],2);
data.peakVisualFast_K_N = data.peakVisualFast_K-data.rw_baselineFiring;

data.visualOptoBaseline = mean(data.PSTHvisualOpto(:,(PSTHbins>-150 & PSTHbins<0)),2);
data.peakVisualOptoFast = max(data.PSTHvisualOpto(:,(PSTHbins>0 & PSTHbins<150)),[],2);
data.peakVisualOptoFast_N = data.peakVisualOptoFast-data.visualOptoBaseline;

data.peakVisualOptoFast_Change = log2((data.peakVisualOptoFast_N+1)./(data.peakVisualFast_N+1));

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

data.peakVisualOptoSlow_Change = log2((data.rwOpto_meanPSTH)./(data.rw_meanPSTH));

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

%% Plot scatter plot depth waveform
figure('units','normalized','outerposition',[0 0 0.2 1]);
hold on
plot(data.troughPeakTime(L23,:), data.Depth(L23,:),'Color',[150,150,150]/255,'LineStyle','none','Marker','o','MarkerFaceColor',[150,150,150]/255)
plot(data.troughPeakTime(L4,:), data.Depth(L4,:),'Color',[150,150,150]/255,'LineStyle','none','Marker','o','MarkerFaceColor',[150,150,150]/255)
plot(data.troughPeakTime(L5,:), data.Depth(L5,:),'Color',[150,150,150]/255,'LineStyle','none','Marker','o','MarkerFaceColor',[150,150,150]/255)
plot(data.troughPeakTime(CT,:), data.Depth(CT,:),'Color',[215,25,28]/255,'LineStyle','none','Marker','o','MarkerFaceColor',[215,25,28]/255)
plot(data.troughPeakTime(CC,:), data.Depth(CC,:),'Color',[44,123,182]/255,'LineStyle','none','Marker','o','MarkerFaceColor',[44,123,182]/255)
plot([0.75,0.75],[-700,300],'k--')
ax=gca;
% ax.Color = [240,240,240]/255;
% print(gcf,'-dpdf','C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\WT Doctoral Programme in Neuroscience\DPhil in Neuroscience\DPhil thesis\Figures\Chapter 5\5.10\ScatterPlot.pdf')

%% Plot waveform Ntsr1
x = readNPY('G:\InVivo_SpikeSorting_V1\NTV1\waveforms\82.npy');
time = linspace(0,200/30,size(x,2));
x=squeeze(x(:,:,1));
[~,Idx]=max(max(x,[],2));
x(Idx,:) = [];
[~,Idx]=min(min(x,[],2));
x(Idx,:) = [];

figure
hold on
plot(time,x(1:1000,:),'Color',[200,200,200]/255);
plot(time,mean(x,1),'Color',[215,25,28]/255);

plot([4,4],[600,1200],'k')
plot([4,4.5],[600,600],'k')
xlim([2 5])
ylim([-1000,1000])
export_fig('C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\WT Doctoral Programme in Neuroscience\DPhil in Neuroscience\DPhil thesis\Figures\Chapter 5\5.10\U82_waveform.pdf','-pdf','-transparent','-nocrop')

FS = readNPY('G:\InVivo_SpikeSorting_V1\NTV1\waveforms\97.npy');
NTSRcell = readNPY('G:\InVivo_SpikeSorting_V1\NTV1\waveforms\119.npy');
time = linspace(0,200/30,size(FS,2));

FS=squeeze(FS(:,:,1));
NTSRcell=squeeze(NTSRcell(:,:,1));

figure
hold on
plot(time,mean(NTSRcell,1),'Color',[215,25,28]/255);
plot(time,mean(FS,1),'Color',[37 37 37]/255);

xlim([2 5])
ylim([-800,400])

plot([4,4],[-600,-400],'k')
plot([4,4.5],[-600,-600],'k')

export_fig('C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\WT Doctoral Programme in Neuroscience\DPhil in Neuroscience\DPhil thesis\Figures\Chapter 5\5.11\U97-119_waveform.pdf','-pdf','-transparent','-nocrop')

%% PSTH plot all single units
smoothingVal=5;

subset=data(CC&P9P13,:);
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


%% SB bar plots
clear ax
yHist=[nnz(L23&SB_entrained&P9P13)/nnz(L23&P9P13),nnz(L4&SB_entrained&P9P13)/nnz(L4&P9P13),nnz(L5&SB_entrained&P9P13)/nnz(L5&P9P13),nnz(CC&SB_entrained&P9P13)/nnz(CC&P9P13),nnz(CT&SB_entrained&P9P13)/nnz(CT&P9P13);...
    nnz(L23&PPC_entrained&P9P13)/nnz(L23&P9P13),nnz(L4&PPC_entrained&P9P13)/nnz(L4&P9P13),nnz(L5&PPC_entrained&P9P13)/nnz(L5&P9P13),nnz(CC&PPC_entrained&P9P13)/nnz(CC&P9P13),nnz(CT&PPC_entrained&P9P13)/nnz(CT&P9P13)];

xHist=reordercats(categorical(cellstr({'L2/3','L4','L5','CC','CT'})),{'L2/3','L4','L5','CC','CT'});
figure('units','normalized','outerposition',[0 0 0.2 1]);
ax=subplot(3,1,1);
b=bar(xHist,yHist);
b(1,1).LineWidth=1;
b(1,2).LineWidth=1;
ax=gca;
% ax.XLabel.String='Postnatal day';
ax.YLabel.String='Ratio entrained single units';
ax.Box='off';
ax.LineWidth = 1.5;
ax.FontSize=10;
ax.YLim=[0,1];
legend('SB','PPC','Location','northwest')
legend('boxoff')

%% Responsive bar plots
clear ax
yHist=[nnz(L23&resp1&P9P13)/nnz(L23&P9P13),nnz(L4&resp1&P9P13)/nnz(L4&P9P13),nnz(L5&resp1&P9P13)/nnz(L5&P9P13),nnz(CC&resp1&P9P13)/nnz(CC&P9P13),nnz(CT&resp1&P9P13)/nnz(CT&P9P13);...
    nnz(L23&resp2&P9P13)/nnz(L23&P9P13),nnz(L4&resp2&P9P13)/nnz(L4&P9P13),nnz(L5&resp2&P9P13)/nnz(L5&P9P13),nnz(CC&resp2&P9P13)/nnz(CC&P9P13),nnz(CT&resp2&P9P13)/nnz(CT&P9P13)];

xHist=reordercats(categorical(cellstr({'L2/3','L4','L5','CC','CT'})),{'L2/3','L4','L5','CC','CT'});
figure('units','normalized','outerposition',[0 0 0.2 1]);
ax=subplot(3,1,1);
b=bar(xHist,yHist);
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

clear ax
yHist=[nnz(L23&resp1&P14P18)/nnz(L23&P14P18),nnz(L4&resp1&P14P18)/nnz(L4&P14P18),nnz(L5&resp1&P14P18)/nnz(L5&P14P18),nnz(CC&resp1&P14P18)/nnz(CC&P14P18),nnz(CT&resp1&P14P18)/nnz(CT&P14P18),nnz(FS&resp1&P14P18)/nnz(FS&P14P18);...
    nnz(L23&resp2&P14P18)/nnz(L23&P14P18),nnz(L4&resp2&P14P18)/nnz(L4&P14P18),nnz(L5&resp2&P14P18)/nnz(L5&P14P18),nnz(CC&resp2&P14P18)/nnz(CC&P14P18),nnz(CT&resp2&P14P18)/nnz(CT&P14P18),nnz(FS&resp2&P14P18)/nnz(FS&P14P18)];

xHist=reordercats(categorical(cellstr({'L2/3','L4','L5','CC','CT','FS'})),{'L2/3','L4','L5','CC','CT','FS'});

% figure('units','normalized','outerposition',[0 0 0.2 1]);
ax=subplot(3,1,2);
b=bar(xHist,yHist);
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

export_fig(fullfile(folderFigures,'5.12','BarPlots'),'-pdf','-transparent','-nocrop')

%%
figure('units','normalized','outerposition',[0 0 1 1]);
smoothingVal=5;
examplePSTH=[data.PSTHvisual(data.MouseID=='NTV19' & data.suid==71,:);data.PSTHvisual(data.MouseID=='NTV19' & data.suid==28,:);data.PSTHvisual(data.MouseID=='NTV9' & data.suid==27,:);data.PSTHvisual(data.MouseID=='NTV9' & data.suid==33,:)];
titlesExamplePSTH={'L2/3','L4','L6 CT','L6 CC'};
for i=1:4
    ax=subplot(2,2,i);
    hold on
%     patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(examplePSTH(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,49];
    ax.XLim=[-500,3000];
    ax.Title.String=titlesExamplePSTH{i};
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('P9-P13','FontSize',25)
export_fig(fullfile(folderFigures,'5.12','PSTH-RS-P9P13'),'-pdf','-transparent','-nocrop')
close

figure('units','normalized','outerposition',[0 0 1 1]);
smoothingVal=5;
examplePSTH=[data.PSTHvisual(data.MouseID=='NTV11' & data.suid==9,:);data.PSTHvisual(data.MouseID=='NTV1' & data.suid==111,:);data.PSTHvisual(data.MouseID=='NTV1' & data.suid==82,:);data.PSTHvisual(data.MouseID=='NTV1' & data.suid==122,:)];
titlesExamplePSTH={'L2/3','L4','L6 CT','L6 CC'};
for i=1:4
    ax=subplot(2,2,i);
    hold on
%     patch([0 100 100 0],[-2 -2 100 100],'y','EdgeColor','none','HandleVisibility','off') 
    plot(PSTHbins,smooth(examplePSTH(i,:),smoothingVal),'k','LineWidth',2)
	ax=applyFont(ax,0);
    ax.XLabel.String=[];
    ax.YLabel.String=[];
    ax.YLim=[-3,49];
    ax.XLim=[-500,3000];
    ax.Title.String=titlesExamplePSTH{i};
    ax.Title.FontSize=16;
    ax.Color=[247,247,247]/255;
end
sgtitle('P14-P18','FontSize',25)
export_fig(fullfile(folderFigures,'5.12','PSTH-RS-P14P18'),'-pdf','-transparent','-nocrop')
close


%% Violin resp
clear ax v
figure('units','normalized','outerposition',[0 0 0.2 1])
ax(1)=subplot(3,2,1);
v(1,:)=violinplot(data.peakVisualFast_N(P9P13&resp1),data.CellIdentity(P9P13&resp1));
ax(1).YLim=[0, 150];
ax(1).YAxis.Label.String='Max fast PSTH (spike/s)';

ax(2)=subplot(3,2,2);
v(2,:)=violinplot(data.peakVisualFast_N(P14P18&resp1),data.CellIdentity(P14P18&resp1));
ax(2).YLim=[0, 150];

[p,tbl,stats]=anovan(data.peakVisualFast_N(RS&responsive&(P9P13|P14P18)),{data.Layer(RS&responsive&(P9P13|P14P18)),removecats(data.Dev(RS&responsive&(P9P13|P14P18)))},'model','interaction','varnames',{'Cell Layer','Development'});
multcompare(stats,'Dimension',[1])

ax(3)=subplot(3,2,3);
v(3,:)=violinplot(data.fanoVisual(P9P13&resp1),data.CellIdentity(P9P13&resp1));
% ax(3).YLim=[-0.15, 1.2];
ax(3).YAxis.Label.String='Fano factor';

ax(4)=subplot(3,2,4);
v(4,:)=violinplot(data.fanoVisual(P14P18&resp1),data.CellIdentity(P14P18&resp1));
% ax(4).YLim=[-0.15, 1.2];

[p,tbl,stats]=anovan(data.fanoVisual(RS&responsive&(P9P13|P14P18)),{data.Layer(RS&responsive&(P9P13|P14P18)),removecats(data.Dev(RS&responsive&(P9P13|P14P18)))},'model','interaction','varnames',{'Cell Layer','Development'});
figure
multcompare(stats,'Dimension',[1,2])

ax(5)=subplot(3,2,5);
v(5,:)=violinplot(data.rw_firing_N(P9P13&resp2),data.CellIdentity(P9P13&resp2));
ax(5).YLim=[0, 10];
ax(5).YAxis.Label.String='Average slow PSTH (spike/s)';

ax(6)=subplot(3,2,6);
v(6,:)=violinplot(data.rw_firing_N(P14P18&resp2),data.CellIdentity(P14P18&resp2));
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

%% Opto baseline
clear ax
figure('units','normalized','outerposition',[0 0 0.2 1]);
ax(1)=subplot(3,1,1);
hold on
plot(-3,data.optoChange(RS&P9P13&L23&isLaser),'o','Color',[150,150,150]/255)
errorbar(-3+0.2,mean(data.optoChange(RS&P9P13&L23&isLaser),'omitnan'),sem(data.optoChange(RS&P9P13&L23&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P9P13&L23&isLaser)) 

plot(-2,data.optoChange(RS&P9P13&L4&isLaser),'o','Color',[150,150,150]/255)
errorbar(-2+0.2,mean(data.optoChange(RS&P9P13&L4&isLaser),'omitnan'),sem(data.optoChange(RS&P9P13&L4&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P9P13&L4&isLaser))

plot(-1,data.optoChange(RS&P9P13&CC&isLaser),'o','Color',[150,150,150]/255)
errorbar(-1+0.2,mean(data.optoChange(RS&P9P13&CC&isLaser),'omitnan'),sem(data.optoChange(RS&P9P13&CC&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P9P13&CC&isLaser))


plot(1,data.optoChange(RS&P14P18&L23&isLaser),'o','Color',[150,150,150]/255)
errorbar(1+0.2,mean(data.optoChange(RS&P14P18&L23&isLaser),'omitnan'),sem(data.optoChange(RS&P14P18&L23&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P14P18&L23&isLaser))

plot(2,data.optoChange(RS&P14P18&L4&isLaser),'o','Color',[150,150,150]/255)
errorbar(2+0.2,mean(data.optoChange(RS&P14P18&L4&isLaser),'omitnan'),sem(data.optoChange(RS&P14P18&L4&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P14P18&L4&isLaser))

plot(3,data.optoChange(RS&P14P18&L5&isLaser),'o','Color',[150,150,150]/255)
errorbar(3+0.2,mean(data.optoChange(RS&P14P18&L5&isLaser),'omitnan'),sem(data.optoChange(RS&P14P18&L5&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P14P18&L5&isLaser))

plot(4,data.optoChange(RS&P14P18&CC&isLaser),'o','Color',[150,150,150]/255)
errorbar(4+0.2,mean(data.optoChange(RS&P14P18&CC&isLaser),'omitnan'),sem(data.optoChange(RS&P14P18&CC&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(RS&P14P18&CC&isLaser))

plot(6,data.optoChange(L23&FS&isLaser),'o','Color',[150,150,150]/255)
errorbar(6+0.2,mean(data.optoChange(L23&FS&isLaser),'omitnan'),sem(data.optoChange(L23&FS&isLaser)),'ko','Linewidth',1,'CapSize',10)
% [h,p,stats]=my_ttest(data.optoChange(L23&FS&isLaser))

plot(7,data.optoChange(L4&FS&isLaser),'o','Color',[150,150,150]/255)
errorbar(7+0.2,mean(data.optoChange(L4&FS&isLaser),'omitnan'),sem(data.optoChange(L4&FS&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(L4&FS&isLaser))

plot(8,data.optoChange(L56&FS&isLaser),'o','Color',[150,150,150]/255)
errorbar(8+0.2,mean(data.optoChange(L56&FS&isLaser),'omitnan'),sem(data.optoChange(L56&FS&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.optoChange(L56&FS&isLaser))


ax(2)=subplot(3,1,2);
hold on
plot(-3,data.peakVisualOptoFast_Change(RS&P9P13&L23&isLaser),'o','Color',[150,150,150]/255)
errorbar(-3+0.2,mean(data.peakVisualOptoFast_Change(RS&P9P13&L23&isLaser),'omitnan'),sem(data.peakVisualOptoFast_Change(RS&P9P13&L23&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoFast_Change(RS&P9P13&L23&isLaser)) 

plot(-2,data.peakVisualOptoFast_Change(RS&P9P13&L4&isLaser),'o','Color',[150,150,150]/255)
errorbar(-2+0.2,mean(data.peakVisualOptoFast_Change(RS&P9P13&L4&isLaser),'omitnan'),sem(data.peakVisualOptoFast_Change(RS&P9P13&L4&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoFast_Change(RS&P9P13&L4&isLaser))

plot(-1,data.peakVisualOptoFast_Change(RS&P9P13&CC&isLaser),'o','Color',[150,150,150]/255)
errorbar(-1+0.2,mean(data.peakVisualOptoFast_Change(RS&P9P13&CC&isLaser),'omitnan'),sem(data.peakVisualOptoFast_Change(RS&P9P13&CC&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoFast_Change(RS&P9P13&CC&isLaser))


plot(1,data.peakVisualOptoFast_Change(RS&P14P18&L23&isLaser),'o','Color',[150,150,150]/255)
errorbar(1+0.2,mean(data.peakVisualOptoFast_Change(RS&P14P18&L23&isLaser),'omitnan'),sem(data.peakVisualOptoFast_Change(RS&P14P18&L23&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoFast_Change(RS&P14P18&L23&isLaser))

plot(2,data.peakVisualOptoFast_Change(RS&P14P18&L4&isLaser),'o','Color',[150,150,150]/255)
errorbar(2+0.2,mean(data.peakVisualOptoFast_Change(RS&P14P18&L4&isLaser),'omitnan'),sem(data.peakVisualOptoFast_Change(RS&P14P18&L4&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoFast_Change(RS&P14P18&L4&isLaser))

plot(3,data.peakVisualOptoFast_Change(RS&P14P18&L5&isLaser),'o','Color',[150,150,150]/255)
errorbar(3+0.2,mean(data.peakVisualOptoFast_Change(RS&P14P18&L5&isLaser),'omitnan'),sem(data.peakVisualOptoFast_Change(RS&P14P18&L5&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoFast_Change(RS&P14P18&L5&isLaser&resp1))

plot(4,data.peakVisualOptoFast_Change(RS&P14P18&CC&isLaser),'o','Color',[150,150,150]/255)
errorbar(4+0.2,mean(data.peakVisualOptoFast_Change(RS&P14P18&CC&isLaser),'omitnan'),sem(data.peakVisualOptoFast_Change(RS&P14P18&CC&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoFast_Change(RS&P14P18&CC&isLaser))

plot(6,data.peakVisualOptoFast_Change(L23&FS&isLaser),'o','Color',[150,150,150]/255)
errorbar(6+0.2,mean(data.peakVisualOptoFast_Change(L23&FS&isLaser),'omitnan'),sem(data.peakVisualOptoFast_Change(L23&FS&isLaser)),'ko','Linewidth',1,'CapSize',10)
% [h,p,stats]=my_ttest(data.peakVisualOptoFast_Change(L23&FS&isLaser))

plot(7,data.peakVisualOptoFast_Change(L4&FS&isLaser),'o','Color',[150,150,150]/255)
errorbar(7+0.2,mean(data.peakVisualOptoFast_Change(L4&FS&isLaser),'omitnan'),sem(data.peakVisualOptoFast_Change(L4&FS&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoFast_Change(L4&FS&isLaser))

plot(8,data.peakVisualOptoFast_Change(L56&FS&isLaser),'o','Color',[150,150,150]/255)
errorbar(8+0.2,mean(data.peakVisualOptoFast_Change(L56&FS&isLaser),'omitnan'),sem(data.peakVisualOptoFast_Change(L56&FS&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoFast_Change(L56&FS&isLaser))

ax(3)=subplot(3,1,3);
hold on
plot(-3,data.peakVisualOptoSlow_Change(RS&P9P13&L23&isLaser),'o','Color',[150,150,150]/255)
errorbar(-3+0.2,mean(data.peakVisualOptoSlow_Change(RS&P9P13&L23&isLaser),'omitnan'),sem(data.peakVisualOptoSlow_Change(RS&P9P13&L23&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoSlow_Change(RS&P9P13&L23&isLaser)) 

plot(-2,data.peakVisualOptoSlow_Change(RS&P9P13&L4&isLaser),'o','Color',[150,150,150]/255)
errorbar(-2+0.2,mean(data.peakVisualOptoSlow_Change(RS&P9P13&L4&isLaser),'omitnan'),sem(data.peakVisualOptoSlow_Change(RS&P9P13&L4&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoSlow_Change(RS&P9P13&L4&isLaser))

plot(-1,data.peakVisualOptoSlow_Change(RS&P9P13&CC&isLaser),'o','Color',[150,150,150]/255)
errorbar(-1+0.2,mean(data.peakVisualOptoSlow_Change(RS&P9P13&CC&isLaser),'omitnan'),sem(data.peakVisualOptoSlow_Change(RS&P9P13&CC&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoSlow_Change(RS&P9P13&CC&isLaser))


plot(1,data.peakVisualOptoSlow_Change(RS&P14P18&L23&isLaser),'o','Color',[150,150,150]/255)
errorbar(1+0.2,mean(data.peakVisualOptoSlow_Change(RS&P14P18&L23&isLaser),'omitnan'),sem(data.peakVisualOptoSlow_Change(RS&P14P18&L23&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoSlow_Change(RS&P14P18&L23&isLaser))

plot(2,data.peakVisualOptoSlow_Change(RS&P14P18&L4&isLaser),'o','Color',[150,150,150]/255)
errorbar(2+0.2,mean(data.peakVisualOptoSlow_Change(RS&P14P18&L4&isLaser),'omitnan'),sem(data.peakVisualOptoSlow_Change(RS&P14P18&L4&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoSlow_Change(RS&P14P18&L4&isLaser))

plot(3,data.peakVisualOptoSlow_Change(RS&P14P18&L5&isLaser),'o','Color',[150,150,150]/255)
errorbar(3+0.2,mean(data.peakVisualOptoSlow_Change(RS&P14P18&L5&isLaser),'omitnan'),sem(data.peakVisualOptoSlow_Change(RS&P14P18&L5&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoSlow_Change(RS&P14P18&L5&isLaser))

plot(4,data.peakVisualOptoSlow_Change(RS&P14P18&CC&isLaser),'o','Color',[150,150,150]/255)
errorbar(4+0.2,mean(data.peakVisualOptoSlow_Change(RS&P14P18&CC&isLaser),'omitnan'),sem(data.peakVisualOptoSlow_Change(RS&P14P18&CC&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoSlow_Change(RS&P14P18&CC&isLaser))

plot(6,data.peakVisualOptoSlow_Change(L23&FS&isLaser),'o','Color',[150,150,150]/255)
errorbar(6+0.2,mean(data.peakVisualOptoSlow_Change(L23&FS&isLaser),'omitnan'),sem(data.peakVisualOptoSlow_Change(L23&FS&isLaser)),'ko','Linewidth',1,'CapSize',10)
% [h,p,stats]=my_ttest(data.peakVisualOptoSlow_Change(L23&FS&isLaser))

plot(7,data.peakVisualOptoSlow_Change(L4&FS&isLaser),'o','Color',[150,150,150]/255)
errorbar(7+0.2,mean(data.peakVisualOptoSlow_Change(L4&FS&isLaser),'omitnan'),sem(data.peakVisualOptoSlow_Change(L4&FS&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoSlow_Change(L4&FS&isLaser))

plot(8,data.peakVisualOptoSlow_Change(L56&FS&isLaser),'o','Color',[150,150,150]/255)
errorbar(8+0.2,mean(data.peakVisualOptoSlow_Change(L56&FS&isLaser),'omitnan'),sem(data.peakVisualOptoSlow_Change(L56&FS&isLaser)),'ko','Linewidth',1,'CapSize',10)
[h,p,stats]=my_ttest(data.peakVisualOptoSlow_Change(L56&FS&isLaser))

for i=1:numel(ax)
    ax(i).XLim=[-3.5,8.5];
    ax(i).YAxis.Label.String='Mean firing change (ratio log2)';
    ax(i).FontSize=10;
    ax(i).YLim=[-4,4];
    ax(i).XTick=[-3,-2,-1,1,2,3,4,6,7,8];
    ax(i).XTickLabel={'L2/3','L4','L5/6','L2/3','L4','L5/6','CC','L2/3','L4','L5/6'};
end
ax(1).YLim=[-3,3];

export_fig(fullfile(folderFigures,'5.13','ChangePlots'),'-pdf','-transparent','-nocrop')
