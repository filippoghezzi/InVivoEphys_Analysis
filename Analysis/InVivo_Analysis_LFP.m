close all
clear 
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 
load('LFPData.mat')
folderFigures='C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\WT Doctoral Programme in Neuroscience\DPhil in Neuroscience\DPhil thesis\Figures\Chapter 4';


%% Set group logic arrays
data=data(data.state=='Urethane',:);
data=data(data.brainArea=='S1BF',:);
data=data(data.Tagging~='SST;NrgKO',:);
% data(data.MouseID=='NTV10',:)=[];
% data(data.MouseID=='NTV24',:)=[];
% data(data.MouseID=='NK14',:)=[];
% data(data.MouseID=='NK15',:)=[];
% data(data.MouseID=='SC35',:)=[];

%Age
P5P8=data.mouseAge<9;
P9P13=data.mouseAge>=9 & data.mouseAge<14;
P14P18=data.mouseAge>=14;
Dev=categorical;
Dev(P5P8,1)='P5-P8';
Dev(P9P13,1)='P9-P13';
Dev(P14P18,1)='P14-P18';

data.PPR=data.MUA_L4_peakSlow./data.MUA_L4_peakFast;
data.PPR_K=data.MUA_L4_peakSlow_K./data.MUA_L4_peakFast_K;

    %% PSD plots
figure('units','normalized','outerposition',[0 0 0.2 1])

ax(1) = subplot(3,1,1);
plot(PSD_f,data.PSD(data.MouseID=='SC101b',:),'k','LineWidth',1)
hold on
plot(PSD_f,data.PSD(data.MouseID=='SC99b',:),'r','LineWidth',1)

ax(2) = subplot(3,1,2);
plot(PSD_f,data.PSD(data.MouseID=='NK32',:),'k','LineWidth',1)
hold on
plot(PSD_f,data.PSD(data.MouseID=='SC99b',:),'r','LineWidth',1)

ax(3) = subplot(3,1,3);
plot(PSD_f,data.PSD(data.MouseID=='SC46',:),'k','LineWidth',1)

for i=1:numel(ax)
    ax(i).XLabel.String='Frequency (Hz)';
    ax(i).YLabel.String='Normalised PSD';
    ax(i).FontSize=12;
    ax(i).FontName='Arial';
    ax(i).Box='off';
    ax(i).LineWidth=1;
    ax(i).XLim=[1,50];
    ax(i).YLim=[0,5];
end

export_fig(fullfile(folderFigures,'4.1','PSD'),'-pdf','-transparent','-nocrop')
close

maxBetaPSD=max(data.PSD(:,(PSD_f>10) & (PSD_f<30)),[],2);


figure('units','normalized','outerposition',[0 0 0.2 1])

ax(1)=subplot(3,1,1);
v(1,:)=violinplot(maxBetaPSD, Dev);
ax(1).YLabel.String='Max \beta normalised power';
% [p,tbl,stats]=anovan(maxBetaPSD,{Dev},'model','linear','varnames',{'Development'});
% multcompare(stats)

ax(2)=subplot(3,1,2);
v(2,:)=violinplot(data.SB_L4_duration, Dev);
ax(2).YLabel.String='Spindle burst duration (s)';
% [p,tbl,stats]=anovan(data.SB_L4_duration,{Dev},'model','linear','varnames',{'Development'});
% multcompare(stats)

ax(3)=subplot(3,1,3);
v(3,:)=violinplot(data.SB_L4_frequency, Dev);
ax(3).YLabel.String='Spindle burst occurrence (Hz)';
% [p,tbl,stats]=anovan(data.SB_L4_frequency,{Dev},'model','linear','varnames',{'Development'});
% multcompare(stats)

for i=1:numel(ax)
    for j=1:size(v,2)
        v(i,j).ViolinColor=[100,100,100]/255;
        v(i,j).ScatterPlot.MarkerFaceColor=[37,37,37]/255;
        v(i,j).ScatterPlot.MarkerFaceAlpha=1;
        v(i,j).ScatterPlot.SizeData=10;
    end

    ax(i).FontSize=12;
    ax(i).LineWidth=1;
    ax(i).XLim=[0.5,3.5];
end

print(gcf, '-dpdf', fullfile(folderFigures,'4.1','PSD_SB_violin'))
close

%% LFP evoked
figure('units','normalized','outerposition',[0 0 0.2 1])
MeanLFPlatency=[];
SEMLFPlatency=[];

for i=9:18
    MeanLFPlatency(i-8) = mean(data.LFP_L4_latency(data.mouseAge==i));
    SEMLFPlatency(i-8) = sem(data.LFP_L4_latency(data.mouseAge==i));
    MeanLFPamplitude(i-8) = mean(data.LFP_L4_peak(data.mouseAge==i));
    SEMLFPamplitude(i-8) = sem(data.LFP_L4_peak(data.mouseAge==i));
    MeanMUA_L4_peakFast(i-8) = mean(data.MUA_L4_peakFast(data.mouseAge==i));
    SEMMUA_L4_peakFast(i-8) = sem(data.MUA_L4_peakFast(data.mouseAge==i));
    MeanMUA_L4_peakSlow(i-8) = mean(data.MUA_L4_peakSlow(data.mouseAge==i));
    SEMMUA_L4_peakSlow(i-8) = sem(data.MUA_L4_peakSlow(data.mouseAge==i));
end
  
ax(1)=subplot(3,1,1);
scatter(data.mouseAge(P9P13|P14P18),data.LFP_L4_latency(P9P13|P14P18),'MarkerEdgeColor',[150,150,150]/255,'MarkerFaceColor','none','SizeData',18);
hold on
errorbar((9:18)+0.3,MeanLFPlatency,SEMLFPlatency,'ko-','MarkerSize',6,'LineWidth',0.75)
ax(1).XLim=[8.5,19];
ax(1).XLabel.String='Age (days)';
ax(1).YLabel.String='Peak evoked LFP latency (ms)';
ax(1).FontSize=12;
ax(1).FontName='Arial';
ax(1).Box='off';
ax(1).LineWidth=1;
ax(1).YLim=[0,1000];

swtest(data.LFP_L4_latency(P9P13|P14P18))
% [p,tbl,stats]=kruskalwallis(data.LFP_L4_latency(P9P13|P14P18),data.mouseAge(P9P13|P14P18));
% multcompare(stats,'CType' ,'dunn-sidak')

ax(2)=subplot(3,1,2);
scatter(data.mouseAge(P9P13|P14P18),data.LFP_L4_peak(P9P13|P14P18),'MarkerEdgeColor',[150,150,150]/255,'MarkerFaceColor','none','SizeData',18);
hold on
errorbar((9:18)+0.3,MeanLFPamplitude,SEMLFPamplitude,'ko-','MarkerSize',6,'LineWidth',0.75)
ax(2).XLim=[8.5,19];
ax(2).XLabel.String='Age (days)';
ax(2).YLabel.String='Peak evoked LFP amplitude (\muV)';
ax(2).FontSize=12;
ax(2).FontName='Arial';
ax(2).Box='off';
ax(2).LineWidth=1;

swtest(data.LFP_L4_peak(P9P13|P14P18))
% [p,tbl,stats]=kruskalwallis(data.LFP_L4_peak(P9P13|P14P18),data.mouseAge(P9P13|P14P18));
% multcompare(stats,'CType' ,'dunn-sidak')

ax(3)=subplot(3,1,3);
errorbar((9:18),MeanMUA_L4_peakFast,SEMMUA_L4_peakFast,'ko-','MarkerSize',6,'LineWidth',0.75)
hold on
errorbar((9:18),MeanMUA_L4_peakSlow,SEMMUA_L4_peakSlow,'o-','Color',[150,150,150]/255,'MarkerSize',6,'LineWidth',0.75)
ax(3).XLim=[8.5,19];
ax(3).XLabel.String='Age (days)';
ax(3).YLabel.String='L4 MUA amplitude (spikes/s)';
ax(3).FontSize=12;
ax(3).FontName='Arial';
ax(3).Box='off';
ax(3).LineWidth=1;
legend('Fast','Slow','Location','northwest')
legend('box','off')

% swtest(data.MUA_L4_peakFast(P9P13|P14P18))
% [p,tbl,stats]=kruskalwallis(data.MUA_L4_peakFast(P9P13|P14P18),data.mouseAge(P9P13|P14P18));
% multcompare(stats,'CType' ,'dunn-sidak')

% swtest(data.MUA_L4_peakSlow(P9P13|P14P18))
% [p,tbl,stats]=kruskalwallis(data.MUA_L4_peakSlow(P9P13|P14P18),data.mouseAge(P9P13|P14P18));
% multcompare(stats,'CType' ,'dunn-sidak')

print(gcf, '-dpdf', fullfile(folderFigures,'4.2','MUA_LFP_plots'))
close

%% Chemogenetics V1
ChemoRealV1 = {'K18','K24','K35','K36','K38','K39','K43','K44','K45','K47','K6'};
ChemoRealS1 = {'K51','K52','K54','K55','K56','K58','K60','K62','K64','K71','K72','K73','K74'};
ChemoReal = [ChemoRealV1(:)',ChemoRealS1(:)'];
ChemoCtrlGFP = {'K82','K83','K84','K85','K86','K87','K88'}; 
ChemoCtrlSaline = {'K75','K76'};
ChemoNoExpression = {'K48','K4','K9','K10','K17','K20','K21','K23','K37','K41','K42','K46','K50','K53','K59'};

KORD=ismember(data.MouseID,ChemoReal);
GFP=ismember(data.MouseID,ChemoCtrlGFP);
Saline=ismember(data.MouseID,ChemoCtrlSaline);
data.maxBetaPSD=max(data.PSD(:,(PSD_f>10) & (PSD_f<30)),[],2);
data.maxBetaPSD_K=max(data.PSD_K(:,(PSD_f>10) & (PSD_f<30)),[],2);

%% Chemo SB V1
% figure('units','normalized','outerposition',[0 0 0.2 1]);
% ax(1)=subplot(3,2,1);
% hold on
% plot([-3,-2],[data.maxBetaPSD(KORD&P9P13),data.maxBetaPSD_K(KORD&P9P13)],'-','Color',[150,150,150]/255)
% errorbar([-3,-2],[mean(data.maxBetaPSD(KORD&P9P13),'omitnan'),mean(data.maxBetaPSD_K(KORD&P9P13),'omitnan')],[sem(data.maxBetaPSD(KORD&P9P13)),sem(data.maxBetaPSD_K(KORD&P9P13))],'k-o','Linewidth',1,'CapSize',10)
% plot([-1,0],[data.maxBetaPSD(KORD&P14P18),data.maxBetaPSD_K(KORD&P14P18)],'-','Color',[150,150,150]/255)
% errorbar([-1,0],[mean(data.maxBetaPSD(KORD&P14P18),'omitnan'),mean(data.maxBetaPSD_K(KORD&P14P18),'omitnan')],[sem(data.maxBetaPSD(KORD&P14P18)),sem(data.maxBetaPSD_K(KORD&P14P18))],'k-o','Linewidth',1,'CapSize',10)
% ax(1).YAxis.Label.String='Normalised \beta power';
% ax(1).Title.String = 'KORD';
% 
% % ax(2)=subplot(3,2,2);
% hold on
% plot([-3,-2],[data.maxBetaPSD(GFP&P9P13),data.maxBetaPSD_K(GFP&P9P13)],'-','Color',[150,150,150]/255)
% errorbar([-3,-2],[mean(data.maxBetaPSD(GFP&P9P13),'omitnan'),mean(data.maxBetaPSD_K(GFP&P9P13),'omitnan')],[sem(data.maxBetaPSD(GFP&P9P13)),sem(data.maxBetaPSD_K(GFP&P9P13))],'k-o','Linewidth',1,'CapSize',10)
% plot([-1,0],[data.maxBetaPSD(GFP&P14P18),data.maxBetaPSD_K(GFP&P14P18)],'-','Color',[150,150,150]/255)
% errorbar([-1,0],[mean(data.maxBetaPSD(GFP&P14P18),'omitnan'),mean(data.maxBetaPSD_K(GFP&P14P18),'omitnan')],[sem(data.maxBetaPSD(GFP&P14P18)),sem(data.maxBetaPSD_K(GFP&P14P18))],'k-o','Linewidth',1,'CapSize',10)
% ax(2).YAxis.Label.String='Normalised \beta power';
% ax(2).Title.String = 'GFP';
% 
% ax(3)=subplot(3,2,3);
% hold on
% plot([-3,-2],[data.SB_L4_duration(KORD&P9P13),data.SB_L4_duration_K(KORD&P9P13)],'-','Color',[150,150,150]/255)
% errorbar([-3,-2],[mean(data.SB_L4_duration(KORD&P9P13),'omitnan'),mean(data.SB_L4_duration_K(KORD&P9P13),'omitnan')],[sem(data.SB_L4_duration(KORD&P9P13)),sem(data.SB_L4_duration_K(KORD&P9P13))],'k-o','Linewidth',1,'CapSize',10)
% plot([-1,0],[data.SB_L4_duration(KORD&P14P18),data.SB_L4_duration_K(KORD&P14P18)],'-','Color',[150,150,150]/255)
% errorbar([-1,0],[mean(data.SB_L4_duration(KORD&P14P18),'omitnan'),mean(data.SB_L4_duration_K(KORD&P14P18),'omitnan')],[sem(data.SB_L4_duration(KORD&P14P18)),sem(data.SB_L4_duration_K(KORD&P14P18))],'k-o','Linewidth',1,'CapSize',10)
% ax(3).YAxis.Label.String='SB duration (s)';
% 
% ax(4)=subplot(3,2,4);
% hold on
% plot([-3,-2],[data.SB_L4_duration(GFP&P9P13),data.SB_L4_duration_K(GFP&P9P13)],'-','Color',[150,150,150]/255)
% errorbar([-3,-2],[mean(data.SB_L4_duration(GFP&P9P13),'omitnan'),mean(data.SB_L4_duration_K(GFP&P9P13),'omitnan')],[sem(data.SB_L4_duration(GFP&P9P13)),sem(data.SB_L4_duration_K(GFP&P9P13))],'k-o','Linewidth',1,'CapSize',10)
% plot([-1,0],[data.SB_L4_duration(GFP&P14P18),data.SB_L4_duration_K(GFP&P14P18)],'-','Color',[150,150,150]/255)
% errorbar([-1,0],[mean(data.SB_L4_duration(GFP&P14P18),'omitnan'),mean(data.SB_L4_duration_K(GFP&P14P18),'omitnan')],[sem(data.SB_L4_duration(GFP&P14P18)),sem(data.SB_L4_duration_K(GFP&P14P18))],'k-o','Linewidth',1,'CapSize',10)
% ax(4).YAxis.Label.String='SB duration (s)';
% 
% ax(5)=subplot(3,2,5);
% hold on
% plot([-3,-2],[data.SB_L4_frequency(KORD&P9P13),data.SB_L4_frequency_K(KORD&P9P13)],'-','Color',[150,150,150]/255)
% errorbar([-3,-2],[mean(data.SB_L4_frequency(KORD&P9P13),'omitnan'),mean(data.SB_L4_frequency_K(KORD&P9P13),'omitnan')],[sem(data.SB_L4_frequency(KORD&P9P13)),sem(data.SB_L4_frequency_K(KORD&P9P13))],'k-o','Linewidth',1,'CapSize',10)
% plot([-1,0],[data.SB_L4_frequency(KORD&P14P18),data.SB_L4_frequency_K(KORD&P14P18)],'-','Color',[150,150,150]/255)
% errorbar([-1,0],[mean(data.SB_L4_frequency(KORD&P14P18),'omitnan'),mean(data.SB_L4_frequency_K(KORD&P14P18),'omitnan')],[sem(data.SB_L4_frequency(KORD&P14P18)),sem(data.SB_L4_frequency_K(KORD&P14P18))],'k-o','Linewidth',1,'CapSize',10)
% ax(5).YAxis.Label.String='SB frequency (Hz)';
% 
% ax(6)=subplot(3,2,6);
% hold on
% plot([-3,-2],[data.SB_L4_frequency(GFP&P9P13),data.SB_L4_frequency_K(GFP&P9P13)],'-','Color',[150,150,150]/255)
% errorbar([-3,-2],[mean(data.SB_L4_frequency(GFP&P9P13),'omitnan'),mean(data.SB_L4_frequency_K(GFP&P9P13),'omitnan')],[sem(data.SB_L4_frequency(GFP&P9P13)),sem(data.SB_L4_frequency_K(GFP&P9P13))],'k-o','Linewidth',1,'CapSize',10)
% plot([-1,0],[data.SB_L4_frequency(GFP&P14P18),data.SB_L4_frequency_K(GFP&P14P18)],'-','Color',[150,150,150]/255)
% errorbar([-1,0],[mean(data.SB_L4_frequency(GFP&P14P18),'omitnan'),mean(data.SB_L4_frequency_K(GFP&P14P18),'omitnan')],[sem(data.SB_L4_frequency(GFP&P14P18)),sem(data.SB_L4_frequency_K(GFP&P14P18))],'k-o','Linewidth',1,'CapSize',10)
% ax(6).YAxis.Label.String='SB frequency (Hz)';
% 
% ylims=[1,3; 1,3; 0,1.5; 0,1.5; 0,0.15; 0,0.15];
% for i=1:numel(ax)
%     ax(i).XLim=[-3.5,.5];
%     ax(i).FontSize=10;
%     ax(i).YLim=ylims(i,:);
% end
% print(gcf, '-dpdf', fullfile(folderFigures,'4.17','SB_plots_Raw'))
% close 

figure('units','normalized','outerposition',[0 0 0.15 1]);
ax(1)=subplot(3,1,1);
hold on
plot(ones(height(data(KORD&P9P13,:)),1),log2(data.maxBetaPSD_K(KORD&P9P13)./data.maxBetaPSD(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(log2(data.maxBetaPSD_K(KORD&P9P13)./data.maxBetaPSD(KORD&P9P13)),'omitnan'),sem(log2(data.maxBetaPSD_K(KORD&P9P13)./data.maxBetaPSD(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*4,log2(data.maxBetaPSD_K(KORD&P14P18)./data.maxBetaPSD(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.maxBetaPSD_K(KORD&P14P18)./data.maxBetaPSD(KORD&P14P18)),'omitnan'),sem(log2(data.maxBetaPSD_K(KORD&P14P18)./data.maxBetaPSD(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P9P13,:)),1)*2,log2(data.maxBetaPSD_K(GFP&P9P13)./data.maxBetaPSD(GFP&P9P13)),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(log2(data.maxBetaPSD_K(GFP&P9P13)./data.maxBetaPSD(GFP&P9P13)),'omitnan'),sem(log2(data.maxBetaPSD_K(GFP&P9P13)./data.maxBetaPSD(GFP&P9P13))),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18,:)),1)*5,log2(data.maxBetaPSD_K(GFP&P14P18)./data.maxBetaPSD(GFP&P14P18)),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(log2(data.maxBetaPSD_K(GFP&P14P18)./data.maxBetaPSD(GFP&P14P18)),'omitnan'),sem(log2(data.maxBetaPSD_K(GFP&P14P18)./data.maxBetaPSD(GFP&P14P18))),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(1).YAxis.Label.String='Normalised \beta power change (ratio log2)';

[h,p, stats] = my_ttest2(log2(data.maxBetaPSD_K(KORD&P9P13)./data.maxBetaPSD(KORD&P9P13)),log2(data.maxBetaPSD_K(GFP&P9P13)./data.maxBetaPSD(GFP&P9P13)))
[h,p, stats] = my_ttest2(log2(data.maxBetaPSD_K(KORD&P14P18)./data.maxBetaPSD(KORD&P14P18)),log2(data.maxBetaPSD_K(GFP&P14P18)./data.maxBetaPSD(GFP&P14P18)))

ax(2)=subplot(3,1,2);
hold on
plot(ones(height(data(KORD&P9P13,:)),1),log2(data.SB_L4_duration_K(KORD&P9P13)./data.SB_L4_duration(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(log2(data.SB_L4_duration_K(KORD&P9P13)./data.SB_L4_duration(KORD&P9P13)),'omitnan'),sem(log2(data.SB_L4_duration_K(KORD&P9P13)./data.SB_L4_duration(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*4,log2(data.SB_L4_duration_K(KORD&P14P18)./data.SB_L4_duration(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.SB_L4_duration_K(KORD&P14P18)./data.SB_L4_duration(KORD&P14P18)),'omitnan'),sem(log2(data.SB_L4_duration_K(KORD&P14P18)./data.SB_L4_duration(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P9P13,:)),1)*2,log2(data.SB_L4_duration_K(GFP&P9P13)./data.SB_L4_duration(GFP&P9P13)),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(log2(data.SB_L4_duration_K(GFP&P9P13)./data.SB_L4_duration(GFP&P9P13)),'omitnan'),sem(log2(data.SB_L4_duration_K(GFP&P9P13)./data.SB_L4_duration(GFP&P9P13))),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18,:)),1)*5,log2(data.SB_L4_duration_K(GFP&P14P18)./data.SB_L4_duration(GFP&P14P18)),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(log2(data.SB_L4_duration_K(GFP&P14P18)./data.SB_L4_duration(GFP&P14P18)),'omitnan'),sem(log2(data.SB_L4_duration_K(GFP&P14P18)./data.SB_L4_duration(GFP&P14P18))),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(2).YAxis.Label.String='SB duration change (ratio log2)';

[h,p, stats] = my_ttest2(log2(data.SB_L4_duration_K(KORD&P9P13)./data.SB_L4_duration(KORD&P9P13)),log2(data.SB_L4_duration_K(GFP&P9P13)./data.SB_L4_duration(GFP&P9P13)))

ax(3)=subplot(3,1,3);
hold on
plot(ones(height(data(KORD&P9P13,:)),1),log2(data.SB_L4_frequency_K(KORD&P9P13)./data.SB_L4_frequency(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(log2(data.SB_L4_frequency_K(KORD&P9P13)./data.SB_L4_frequency(KORD&P9P13)),'omitnan'),sem(log2(data.SB_L4_frequency_K(KORD&P9P13)./data.SB_L4_frequency(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*4,log2(data.SB_L4_frequency_K(KORD&P14P18)./data.SB_L4_frequency(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.SB_L4_frequency_K(KORD&P14P18)./data.SB_L4_frequency(KORD&P14P18)),'omitnan'),sem(log2(data.SB_L4_frequency_K(KORD&P14P18)./data.SB_L4_frequency(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P9P13,:)),1)*2,log2(data.SB_L4_frequency_K(GFP&P9P13)./data.SB_L4_frequency(GFP&P9P13)),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(log2(data.SB_L4_frequency_K(GFP&P9P13)./data.SB_L4_frequency(GFP&P9P13)),'omitnan'),sem(log2(data.SB_L4_frequency_K(GFP&P9P13)./data.SB_L4_frequency(GFP&P9P13))),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18,:)),1)*5,log2(data.SB_L4_frequency_K(GFP&P14P18)./data.SB_L4_frequency(GFP&P14P18)),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(log2(data.SB_L4_frequency_K(GFP&P14P18)./data.SB_L4_frequency(GFP&P14P18)),'omitnan'),sem(log2(data.SB_L4_frequency_K(GFP&P14P18)./data.SB_L4_frequency(GFP&P14P18))),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(3).YAxis.Label.String='SB frequency change (ratio log2)';

[h,p, stats] = my_ttest2(log2(data.SB_L4_frequency_K(KORD&P9P13)./data.SB_L4_frequency(KORD&P9P13)),log2(data.SB_L4_frequency_K(GFP&P9P13)./data.SB_L4_frequency(GFP&P9P13)))

for i=1:numel(ax)
    ax(i).XLim=[0.5,5.5];
    ax(i).FontSize=10;
    ax(i).YLim=[-1,1];
end
print(gcf, '-dpdf', fullfile(folderFigures,'4.17','SB_plots_log2Ratio'))
close 

%% Chemo L4 LFP MUA V1
figure('units','normalized','outerposition',[0 0 0.15 1]);
ax(1)=subplot(3,1,1);
hold on
plot(ones(height(data(KORD&P9P13,:)),1),log2(data.MUA_L4_peakFast_K(KORD&P9P13)./data.MUA_L4_peakFast(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(log2(data.MUA_L4_peakFast_K(KORD&P9P13)./data.MUA_L4_peakFast(KORD&P9P13)),'omitnan'),sem(log2(data.MUA_L4_peakFast_K(KORD&P9P13)./data.MUA_L4_peakFast(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*4,log2(data.MUA_L4_peakFast_K(KORD&P14P18)./data.MUA_L4_peakFast(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.MUA_L4_peakFast_K(KORD&P14P18)./data.MUA_L4_peakFast(KORD&P14P18)),'omitnan'),sem(log2(data.MUA_L4_peakFast_K(KORD&P14P18)./data.MUA_L4_peakFast(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P9P13,:)),1)*2,log2(data.MUA_L4_peakFast_K(GFP&P9P13)./data.MUA_L4_peakFast(GFP&P9P13)),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(log2(data.MUA_L4_peakFast_K(GFP&P9P13)./data.MUA_L4_peakFast(GFP&P9P13)),'omitnan'),sem(log2(data.MUA_L4_peakFast_K(GFP&P9P13)./data.MUA_L4_peakFast(GFP&P9P13))),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18,:)),1)*5,log2(data.MUA_L4_peakFast_K(GFP&P14P18)./data.MUA_L4_peakFast(GFP&P14P18)),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(log2(data.MUA_L4_peakFast_K(GFP&P14P18)./data.MUA_L4_peakFast(GFP&P14P18)),'omitnan'),sem(log2(data.MUA_L4_peakFast_K(GFP&P14P18)./data.MUA_L4_peakFast(GFP&P14P18))),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(1).YAxis.Label.String='MUA L4 fast peak change (ratio log2)';

[h,p, stats] = my_ttest2(log2(data.MUA_L4_peakFast_K(KORD&P9P13)./data.MUA_L4_peakFast(KORD&P9P13)),log2(data.MUA_L4_peakFast_K(GFP&P9P13)./data.MUA_L4_peakFast(GFP&P9P13)))
[h,p, stats] = my_ttest2(log2(data.MUA_L4_peakFast_K(KORD&P14P18)./data.MUA_L4_peakFast(KORD&P14P18)),log2(data.MUA_L4_peakFast_K(GFP&P14P18)./data.MUA_L4_peakFast(GFP&P14P18)))

ax(2)=subplot(3,1,2);
hold on
plot(ones(height(data(KORD&P9P13,:)),1),log2(data.MUA_L4_peakSlow_K(KORD&P9P13)./data.MUA_L4_peakSlow(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(log2(data.MUA_L4_peakSlow_K(KORD&P9P13)./data.MUA_L4_peakSlow(KORD&P9P13)),'omitnan'),sem(log2(data.MUA_L4_peakSlow_K(KORD&P9P13)./data.MUA_L4_peakSlow(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*4,log2(data.MUA_L4_peakSlow_K(KORD&P14P18)./data.MUA_L4_peakSlow(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.MUA_L4_peakSlow_K(KORD&P14P18)./data.MUA_L4_peakSlow(KORD&P14P18)),'omitnan'),sem(log2(data.MUA_L4_peakSlow_K(KORD&P14P18)./data.MUA_L4_peakSlow(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P9P13,:)),1)*2,log2(data.MUA_L4_peakSlow_K(GFP&P9P13)./data.MUA_L4_peakSlow(GFP&P9P13)),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(log2(data.MUA_L4_peakSlow_K(GFP&P9P13)./data.MUA_L4_peakSlow(GFP&P9P13)),'omitnan'),sem(log2(data.MUA_L4_peakSlow_K(GFP&P9P13)./data.MUA_L4_peakSlow(GFP&P9P13))),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18,:)),1)*5,log2(data.MUA_L4_peakSlow_K(GFP&P14P18)./data.MUA_L4_peakSlow(GFP&P14P18)),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(log2(data.MUA_L4_peakSlow_K(GFP&P14P18)./data.MUA_L4_peakSlow(GFP&P14P18)),'omitnan'),sem(log2(data.MUA_L4_peakSlow_K(GFP&P14P18)./data.MUA_L4_peakSlow(GFP&P14P18))),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(2).YAxis.Label.String='MUA L4 slow peak change (ratio log2)';

[h,p, stats] = my_ttest2(log2(data.MUA_L4_peakSlow_K(KORD&P9P13)./data.MUA_L4_peakSlow(KORD&P9P13)),log2(data.MUA_L4_peakSlow_K(GFP&P9P13)./data.MUA_L4_peakSlow(GFP&P9P13)))
[h,p, stats] = my_ttest2(log2(data.MUA_L4_peakSlow_K(KORD&P14P18)./data.MUA_L4_peakSlow(KORD&P14P18)),log2(data.MUA_L4_peakSlow_K(GFP&P14P18)./data.MUA_L4_peakSlow(GFP&P14P18)))

ax(3)=subplot(3,1,3);
hold on
plot(ones(height(data(KORD&P9P13,:)),1),log2(data.MUA_L4_baselineFiring_K(KORD&P9P13)./data.MUA_L4_baselineFiring(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(log2(data.MUA_L4_baselineFiring_K(KORD&P9P13)./data.MUA_L4_baselineFiring(KORD&P9P13)),'omitnan'),sem(log2(data.MUA_L4_baselineFiring_K(KORD&P9P13)./data.MUA_L4_baselineFiring(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*4,log2(data.MUA_L4_baselineFiring_K(KORD&P14P18)./data.MUA_L4_baselineFiring(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.MUA_L4_baselineFiring_K(KORD&P14P18)./data.MUA_L4_baselineFiring(KORD&P14P18)),'omitnan'),sem(log2(data.MUA_L4_baselineFiring_K(KORD&P14P18)./data.MUA_L4_baselineFiring(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P9P13,:)),1)*2,log2(data.MUA_L4_baselineFiring_K(GFP&P9P13)./data.MUA_L4_baselineFiring(GFP&P9P13)),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(log2(data.MUA_L4_baselineFiring_K(GFP&P9P13)./data.MUA_L4_baselineFiring(GFP&P9P13)),'omitnan'),sem(log2(data.MUA_L4_baselineFiring_K(GFP&P9P13)./data.MUA_L4_baselineFiring(GFP&P9P13))),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18,:)),1)*5,log2(data.MUA_L4_baselineFiring_K(GFP&P14P18)./data.MUA_L4_baselineFiring(GFP&P14P18)),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(log2(data.MUA_L4_baselineFiring_K(GFP&P14P18)./data.MUA_L4_baselineFiring(GFP&P14P18)),'omitnan'),sem(log2(data.MUA_L4_baselineFiring_K(GFP&P14P18)./data.MUA_L4_baselineFiring(GFP&P14P18))),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(3).YAxis.Label.String='MUA L4 baseline change (ratio log2)';

[h,p, stats] = my_ttest2(log2(data.MUA_L4_baselineFiring_K(KORD&P9P13)./data.MUA_L4_baselineFiring(KORD&P9P13)),log2(data.MUA_L4_baselineFiring_K(GFP&P9P13)./data.MUA_L4_baselineFiring(GFP&P9P13)))
[h,p, stats] = my_ttest2(log2(data.MUA_L4_baselineFiring_K(KORD&P14P18)./data.MUA_L4_baselineFiring(KORD&P14P18)),log2(data.MUA_L4_baselineFiring_K(GFP&P14P18)./data.MUA_L4_baselineFiring(GFP&P14P18)))

for i=1:numel(ax)
    ax(i).XLim=[0.5,5.5];
    ax(i).FontSize=10;
    ax(i).YLim=[-1,1];
end
print(gcf, '-dpdf', fullfile(folderFigures,'4.17','MUA_plots_log2Ratio'))
close 

figure('units','normalized','outerposition',[0 0 0.15 1]);
ax(1)=subplot(3,1,1);
hold on
plot(ones(height(data(KORD&P9P13,:)),1),log2(data.MUA_L4_continuity_K(KORD&P9P13)./data.MUA_L4_continuity(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(log2(data.MUA_L4_continuity_K(KORD&P9P13)./data.MUA_L4_continuity(KORD&P9P13)),'omitnan'),sem(log2(data.MUA_L4_continuity_K(KORD&P9P13)./data.MUA_L4_continuity(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*4,log2(data.MUA_L4_continuity_K(KORD&P14P18)./data.MUA_L4_continuity(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.MUA_L4_continuity_K(KORD&P14P18)./data.MUA_L4_continuity(KORD&P14P18)),'omitnan'),sem(log2(data.MUA_L4_continuity_K(KORD&P14P18)./data.MUA_L4_continuity(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P9P13,:)),1)*2,log2(data.MUA_L4_continuity_K(GFP&P9P13)./data.MUA_L4_continuity(GFP&P9P13)),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(log2(data.MUA_L4_continuity_K(GFP&P9P13)./data.MUA_L4_continuity(GFP&P9P13)),'omitnan'),sem(log2(data.MUA_L4_continuity_K(GFP&P9P13)./data.MUA_L4_continuity(GFP&P9P13))),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18,:)),1)*5,log2(data.MUA_L4_continuity_K(GFP&P14P18)./data.MUA_L4_continuity(GFP&P14P18)),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(log2(data.MUA_L4_continuity_K(GFP&P14P18)./data.MUA_L4_continuity(GFP&P14P18)),'omitnan'),sem(log2(data.MUA_L4_continuity_K(GFP&P14P18)./data.MUA_L4_continuity(GFP&P14P18))),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(1).YAxis.Label.String='MUA L4 continuity change (ratio log2)';

ax(2)=subplot(3,1,2);
hold on
plot(ones(height(data(KORD&P9P13,:)),1),log2(data.LFP_L4_peak_K(KORD&P9P13)./data.LFP_L4_peak(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(log2(data.LFP_L4_peak_K(KORD&P9P13)./data.LFP_L4_peak(KORD&P9P13)),'omitnan'),sem(log2(data.LFP_L4_peak_K(KORD&P9P13)./data.LFP_L4_peak(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*4,log2(data.LFP_L4_peak_K(KORD&P14P18)./data.LFP_L4_peak(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.LFP_L4_peak_K(KORD&P14P18)./data.LFP_L4_peak(KORD&P14P18)),'omitnan'),sem(log2(data.LFP_L4_peak_K(KORD&P14P18)./data.LFP_L4_peak(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P9P13,:)),1)*2,log2(data.LFP_L4_peak_K(GFP&P9P13)./data.LFP_L4_peak(GFP&P9P13)),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(log2(data.LFP_L4_peak_K(GFP&P9P13)./data.LFP_L4_peak(GFP&P9P13)),'omitnan'),sem(log2(data.LFP_L4_peak_K(GFP&P9P13)./data.LFP_L4_peak(GFP&P9P13))),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18,:)),1)*5,log2(data.LFP_L4_peak_K(GFP&P14P18)./data.LFP_L4_peak(GFP&P14P18)),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(log2(data.LFP_L4_peak_K(GFP&P14P18)./data.LFP_L4_peak(GFP&P14P18)),'omitnan'),sem(log2(data.LFP_L4_peak_K(GFP&P14P18)./data.LFP_L4_peak(GFP&P14P18))),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(2).YAxis.Label.String='LFP L4 peak change (ratio log2)';

ax(3)=subplot(3,1,3);
hold on
plot(ones(height(data(KORD&P9P13,:)),1),log2(data.LFP_L4_latency_K(KORD&P9P13)./data.LFP_L4_latency(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(1.2,mean(log2(data.LFP_L4_latency_K(KORD&P9P13)./data.LFP_L4_latency(KORD&P9P13)),'omitnan'),sem(log2(data.LFP_L4_latency_K(KORD&P9P13)./data.LFP_L4_latency(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*4,log2(data.LFP_L4_latency_K(KORD&P14P18)./data.LFP_L4_latency(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.LFP_L4_latency_K(KORD&P14P18)./data.LFP_L4_latency(KORD&P14P18)),'omitnan'),sem(log2(data.LFP_L4_latency_K(KORD&P14P18)./data.LFP_L4_latency(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P9P13,:)),1)*2,log2(data.LFP_L4_latency_K(GFP&P9P13)./data.LFP_L4_latency(GFP&P9P13)),'o','Color',[150,150,150]/255)
errorbar(2.2,mean(log2(data.LFP_L4_latency_K(GFP&P9P13)./data.LFP_L4_latency(GFP&P9P13)),'omitnan'),sem(log2(data.LFP_L4_latency_K(GFP&P9P13)./data.LFP_L4_latency(GFP&P9P13))),'ro','Linewidth',1,'CapSize',10)
plot(ones(height(data(GFP&P14P18,:)),1)*5,log2(data.LFP_L4_latency_K(GFP&P14P18)./data.LFP_L4_latency(GFP&P14P18)),'o','Color',[150,150,150]/255)
errorbar(5.2,mean(log2(data.LFP_L4_latency_K(GFP&P14P18)./data.LFP_L4_latency(GFP&P14P18)),'omitnan'),sem(log2(data.LFP_L4_latency_K(GFP&P14P18)./data.LFP_L4_latency(GFP&P14P18))),'ro','Linewidth',1,'CapSize',10)
plot([0,6],[0,0],'k--','LineWidth',0.75)
ax(3).YAxis.Label.String='LFP L4 latency change(ratio log2)';

for i=1:numel(ax)
    ax(i).XLim=[0.5,5.5];
    ax(i).FontSize=10;
    ax(i).YLim=[-1,1];
end
print(gcf, '-dpdf', fullfile(folderFigures,'4.17','LFP_plots_log2Ratio'))
close 

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(3,2,1)
plot(PSTHbins, data.MUA_L4(data.MouseID=='K6',:),'k')
title('Control - P9-P13')
ylim([0 200])
subplot(3,2,2)
plot(PSTHbins, data.MUA_L4_K(data.MouseID=='K6',:),'k')
title('KORD - P9-P13')
ylim([0 200])

subplot(3,2,3)
plot(PSTHbins, data.MUA_L4(data.MouseID=='K43',:),'k')
title('Control - P14-P18')
ylim([0 200])
subplot(3,2,4)
plot(PSTHbins, data.MUA_L4_K(data.MouseID=='K43',:),'k')
title('KORD - P14-P18')
ylim([0 200])

export_fig(fullfile(folderFigures,'4.17','MUA_L4_traces'),'-pdf','-transparent','-nocrop')
close 



%% Chemo SB S1
% figure('units','normalized','outerposition',[0 0 0.2 1]);
% ax(1)=subplot(3,2,1);
% hold on
% plot([-5,-4],[data.maxBetaPSD(KORD&P5P8),data.maxBetaPSD_K(KORD&P5P8)],'-','Color',[150,150,150]/255)
% errorbar([-5,-4],[mean(data.maxBetaPSD(KORD&P5P8),'omitnan'),mean(data.maxBetaPSD_K(KORD&P5P8),'omitnan')],[sem(data.maxBetaPSD(KORD&P5P8)),sem(data.maxBetaPSD_K(KORD&P5P8))],'k-o','Linewidth',1,'CapSize',10)
% plot([-3,-2],[data.maxBetaPSD(KORD&P9P13),data.maxBetaPSD_K(KORD&P9P13)],'-','Color',[150,150,150]/255)
% errorbar([-3,-2],[mean(data.maxBetaPSD(KORD&P9P13),'omitnan'),mean(data.maxBetaPSD_K(KORD&P9P13),'omitnan')],[sem(data.maxBetaPSD(KORD&P9P13)),sem(data.maxBetaPSD_K(KORD&P9P13))],'k-o','Linewidth',1,'CapSize',10)
% plot([-1,0],[data.maxBetaPSD(KORD&P14P18),data.maxBetaPSD_K(KORD&P14P18)],'-','Color',[150,150,150]/255)
% errorbar([-1,0],[mean(data.maxBetaPSD(KORD&P14P18),'omitnan'),mean(data.maxBetaPSD_K(KORD&P14P18),'omitnan')],[sem(data.maxBetaPSD(KORD&P14P18)),sem(data.maxBetaPSD_K(KORD&P14P18))],'k-o','Linewidth',1,'CapSize',10)
% ax(1).YAxis.Label.String='Normalised \beta power';
% ax(1).Title.String = 'KORD';
% 
% ax(2)=subplot(3,2,2);
% hold on
% plot([-3,-2],[data.maxBetaPSD(Saline&P9P13),data.maxBetaPSD_K(Saline&P9P13)],'-','Color',[150,150,150]/255)
% errorbar([-3,-2],[mean(data.maxBetaPSD(Saline&P9P13),'omitnan'),mean(data.maxBetaPSD_K(Saline&P9P13),'omitnan')],[sem(data.maxBetaPSD(Saline&P9P13)),sem(data.maxBetaPSD_K(Saline&P9P13))],'k-o','Linewidth',1,'CapSize',10)
% % plot([-1,0],[data.maxBetaPSD(Saline&P14P18),data.maxBetaPSD_K(Saline&P14P18)],'-','Color',[150,150,150]/255)
% % errorbar([-1,0],[mean(data.maxBetaPSD(Saline&P14P18),'omitnan'),mean(data.maxBetaPSD_K(Saline&P14P18),'omitnan')],[sem(data.maxBetaPSD(Saline&P14P18)),sem(data.maxBetaPSD_K(Saline&P14P18))],'k-o','Linewidth',1,'CapSize',10)
% ax(2).YAxis.Label.String='Normalised \beta power';
% ax(2).Title.String = 'GFP';
% 
% ax(3)=subplot(3,2,3);
% hold on
% plot([-5,-4],[data.SB_L4_duration(KORD&P5P8),data.SB_L4_duration_K(KORD&P5P8)],'-','Color',[150,150,150]/255)
% errorbar([-5,-4],[mean(data.SB_L4_duration(KORD&P5P8),'omitnan'),mean(data.SB_L4_duration_K(KORD&P5P8),'omitnan')],[sem(data.SB_L4_duration(KORD&P5P8)),sem(data.SB_L4_duration_K(KORD&P5P8))],'k-o','Linewidth',1,'CapSize',10)
% plot([-3,-2],[data.SB_L4_duration(KORD&P9P13),data.SB_L4_duration_K(KORD&P9P13)],'-','Color',[150,150,150]/255)
% errorbar([-3,-2],[mean(data.SB_L4_duration(KORD&P9P13),'omitnan'),mean(data.SB_L4_duration_K(KORD&P9P13),'omitnan')],[sem(data.SB_L4_duration(KORD&P9P13)),sem(data.SB_L4_duration_K(KORD&P9P13))],'k-o','Linewidth',1,'CapSize',10)
% plot([-1,0],[data.SB_L4_duration(KORD&P14P18),data.SB_L4_duration_K(KORD&P14P18)],'-','Color',[150,150,150]/255)
% errorbar([-1,0],[mean(data.SB_L4_duration(KORD&P14P18),'omitnan'),mean(data.SB_L4_duration_K(KORD&P14P18),'omitnan')],[sem(data.SB_L4_duration(KORD&P14P18)),sem(data.SB_L4_duration_K(KORD&P14P18))],'k-o','Linewidth',1,'CapSize',10)
% ax(3).YAxis.Label.String='SB duration (s)';
% 
% ax(4)=subplot(3,2,4);
% hold on
% plot([-3,-2],[data.SB_L4_duration(Saline&P9P13),data.SB_L4_duration_K(Saline&P9P13)],'-','Color',[150,150,150]/255)
% errorbar([-3,-2],[mean(data.SB_L4_duration(Saline&P9P13),'omitnan'),mean(data.SB_L4_duration_K(Saline&P9P13),'omitnan')],[sem(data.SB_L4_duration(Saline&P9P13)),sem(data.SB_L4_duration_K(Saline&P9P13))],'k-o','Linewidth',1,'CapSize',10)
% % plot([-1,0],[data.SB_L4_duration(Saline&P14P18),data.SB_L4_duration_K(Saline&P14P18)],'-','Color',[150,150,150]/255)
% % errorbar([-1,0],[mean(data.SB_L4_duration(Saline&P14P18),'omitnan'),mean(data.SB_L4_duration_K(Saline&P14P18),'omitnan')],[sem(data.SB_L4_duration(Saline&P14P18)),sem(data.SB_L4_duration_K(Saline&P14P18))],'k-o','Linewidth',1,'CapSize',10)
% ax(4).YAxis.Label.String='SB duration (s)';
% 
% ax(5)=subplot(3,2,5);
% hold on
% plot([-3,-2],[data.SB_L4_frequency(KORD&P9P13),data.SB_L4_frequency_K(KORD&P9P13)],'-','Color',[150,150,150]/255)
% errorbar([-3,-2],[mean(data.SB_L4_frequency(KORD&P9P13),'omitnan'),mean(data.SB_L4_frequency_K(KORD&P9P13),'omitnan')],[sem(data.SB_L4_frequency(KORD&P9P13)),sem(data.SB_L4_frequency_K(KORD&P9P13))],'k-o','Linewidth',1,'CapSize',10)
% plot([-1,0],[data.SB_L4_frequency(KORD&P14P18),data.SB_L4_frequency_K(KORD&P14P18)],'-','Color',[150,150,150]/255)
% errorbar([-1,0],[mean(data.SB_L4_frequency(KORD&P14P18),'omitnan'),mean(data.SB_L4_frequency_K(KORD&P14P18),'omitnan')],[sem(data.SB_L4_frequency(KORD&P14P18)),sem(data.SB_L4_frequency_K(KORD&P14P18))],'k-o','Linewidth',1,'CapSize',10)
% ax(5).YAxis.Label.String='SB frequency (Hz)';
% 
% ax(6)=subplot(3,2,6);
% hold on
% plot([-3,-2],[data.SB_L4_frequency(Saline&P9P13),data.SB_L4_frequency_K(Saline&P9P13)],'-','Color',[150,150,150]/255)
% errorbar([-3,-2],[mean(data.SB_L4_frequency(Saline&P9P13),'omitnan'),mean(data.SB_L4_frequency_K(Saline&P9P13),'omitnan')],[sem(data.SB_L4_frequency(Saline&P9P13)),sem(data.SB_L4_frequency_K(Saline&P9P13))],'k-o','Linewidth',1,'CapSize',10)
% % plot([-1,0],[data.SB_L4_frequency(Saline&P14P18),data.SB_L4_frequency_K(Saline&P14P18)],'-','Color',[150,150,150]/255)
% % errorbar([-1,0],[mean(data.SB_L4_frequency(Saline&P14P18),'omitnan'),mean(data.SB_L4_frequency_K(Saline&P14P18),'omitnan')],[sem(data.SB_L4_frequency(Saline&P14P18)),sem(data.SB_L4_frequency_K(Saline&P14P18))],'k-o','Linewidth',1,'CapSize',10)
% ax(6).YAxis.Label.String='SB frequency (Hz)';
% 
% ylims=[1,3; 1,3; 0,1.5; 0,1.5; 0,0.15; 0,0.15];
% for i=1:numel(ax)
%     ax(i).XLim=[-3.5,.5];
%     ax(i).FontSize=10;
%     ax(i).YLim=ylims(i,:);
% end


figure('units','normalized','outerposition',[0 0 0.15 1]);
ax(1)=subplot(3,1,1);
hold on
plot(ones(height(data(KORD&P5P8,:)),1)*0,log2(data.maxBetaPSD(KORD&P5P8)./data.maxBetaPSD_K(KORD&P5P8)),'o','Color',[150,150,150]/255)
errorbar(0.2,mean(log2(data.maxBetaPSD(KORD&P5P8)./data.maxBetaPSD_K(KORD&P5P8)),'omitnan'),sem(log2(data.maxBetaPSD(KORD&P5P8)./data.maxBetaPSD_K(KORD&P5P8))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P9P13,:)),1)*3,log2(data.maxBetaPSD(KORD&P9P13)./data.maxBetaPSD_K(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(3.2,mean(log2(data.maxBetaPSD(KORD&P9P13)./data.maxBetaPSD_K(KORD&P9P13)),'omitnan'),sem(log2(data.maxBetaPSD(KORD&P9P13)./data.maxBetaPSD_K(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*6,log2(data.maxBetaPSD(KORD&P14P18)./data.maxBetaPSD_K(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(6.2,mean(log2(data.maxBetaPSD(KORD&P14P18)./data.maxBetaPSD_K(KORD&P14P18)),'omitnan'),sem(log2(data.maxBetaPSD(KORD&P14P18)./data.maxBetaPSD_K(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P9P13,:)),1)*4,log2(data.maxBetaPSD_K(Saline&P9P13)./data.maxBetaPSD(Saline&P9P13)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.maxBetaPSD_K(Saline&P9P13)./data.maxBetaPSD(Saline&P9P13)),'omitnan'),sem(log2(data.maxBetaPSD_K(Saline&P9P13)./data.maxBetaPSD(Saline&P9P13))),'go','Linewidth',1,'CapSize',10)
plot([-1,8],[0,0],'k--','LineWidth',0.75)
ax(1).YAxis.Label.String='Normalised \beta power change (ratio log2)';

ax(2)=subplot(3,1,2);
hold on
plot(ones(height(data(KORD&P5P8,:)),1)*0,log2(data.SB_L4_duration_K(KORD&P5P8)./data.SB_L4_duration(KORD&P5P8)),'o','Color',[150,150,150]/255)
errorbar(0.2,mean(log2(data.SB_L4_duration_K(KORD&P5P8)./data.SB_L4_duration(KORD&P5P8)),'omitnan'),sem(log2(data.SB_L4_duration_K(KORD&P5P8)./data.SB_L4_duration(KORD&P5P8))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P9P13,:)),1)*3,log2(data.SB_L4_duration_K(KORD&P9P13)./data.SB_L4_duration(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(3.2,mean(log2(data.SB_L4_duration_K(KORD&P9P13)./data.SB_L4_duration(KORD&P9P13)),'omitnan'),sem(log2(data.SB_L4_duration_K(KORD&P9P13)./data.SB_L4_duration(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*6,log2(data.SB_L4_duration_K(KORD&P14P18)./data.SB_L4_duration(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(6.2,mean(log2(data.SB_L4_duration_K(KORD&P14P18)./data.SB_L4_duration(KORD&P14P18)),'omitnan'),sem(log2(data.SB_L4_duration_K(KORD&P14P18)./data.SB_L4_duration(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P9P13,:)),1)*4,log2(data.SB_L4_duration_K(Saline&P9P13)./data.SB_L4_duration(Saline&P9P13)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.SB_L4_duration_K(Saline&P9P13)./data.SB_L4_duration(Saline&P9P13)),'omitnan'),sem(log2(data.SB_L4_duration_K(Saline&P9P13)./data.SB_L4_duration(Saline&P9P13))),'go','Linewidth',1,'CapSize',10)
plot([-1,8],[0,0],'k--','LineWidth',0.75)
ax(2).YAxis.Label.String='SB duration change (ratio log2)';

ax(3)=subplot(3,1,3);
hold on
plot(ones(height(data(KORD&P5P8,:)),1)*0,log2(data.SB_L4_frequency_K(KORD&P5P8)./data.SB_L4_frequency(KORD&P5P8)),'o','Color',[150,150,150]/255)
errorbar(0.2,mean(log2(data.SB_L4_frequency_K(KORD&P5P8)./data.SB_L4_frequency(KORD&P5P8)),'omitnan'),sem(log2(data.SB_L4_frequency_K(KORD&P5P8)./data.SB_L4_frequency(KORD&P5P8))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P9P13,:)),1)*3,log2(data.SB_L4_frequency_K(KORD&P9P13)./data.SB_L4_frequency(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(3.2,mean(log2(data.SB_L4_frequency_K(KORD&P9P13)./data.SB_L4_frequency(KORD&P9P13)),'omitnan'),sem(log2(data.SB_L4_frequency_K(KORD&P9P13)./data.SB_L4_frequency(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*6,log2(data.SB_L4_frequency_K(KORD&P14P18)./data.SB_L4_frequency(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(6.2,mean(log2(data.SB_L4_frequency_K(KORD&P14P18)./data.SB_L4_frequency(KORD&P14P18)),'omitnan'),sem(log2(data.SB_L4_frequency_K(KORD&P14P18)./data.SB_L4_frequency(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P9P13,:)),1)*4,log2(data.SB_L4_frequency_K(Saline&P9P13)./data.SB_L4_frequency(Saline&P9P13)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.SB_L4_frequency_K(Saline&P9P13)./data.SB_L4_frequency(Saline&P9P13)),'omitnan'),sem(log2(data.SB_L4_frequency_K(Saline&P9P13)./data.SB_L4_frequency(Saline&P9P13))),'go','Linewidth',1,'CapSize',10)
plot([-1,8],[0,0],'k--','LineWidth',0.75)
ax(3).YAxis.Label.String='SB frequency change (ratio log2)';

for i=1:numel(ax)
    ax(i).XLim=[-0.5,7.5];
    ax(i).FontSize=10;
    ax(i).YLim=[-3,3];
    ax(i).XTick=[0.5,3.5,6.5];
end
print(gcf, '-dpdf', fullfile(folderFigures,'4.19','SB_log2Ratio'))
close 

%% Chemo L4 LFP MUA S1
figure('units','normalized','outerposition',[0 0 0.15 1]);
ax(1)=subplot(3,1,1);
hold on
plot(ones(height(data(KORD&P5P8,:)),1)*0,log2(data.MUA_L4_peakFast_K(KORD&P5P8)./data.MUA_L4_peakFast(KORD&P5P8)),'o','Color',[150,150,150]/255)
errorbar(0.2,mean(log2(data.MUA_L4_peakFast_K(KORD&P5P8)./data.MUA_L4_peakFast(KORD&P5P8)),'omitnan'),sem(log2(data.MUA_L4_peakFast_K(KORD&P5P8)./data.MUA_L4_peakFast(KORD&P5P8))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P9P13,:)),1)*3,log2(data.MUA_L4_peakFast_K(KORD&P9P13)./data.MUA_L4_peakFast(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(3.2,mean(log2(data.MUA_L4_peakFast_K(KORD&P9P13)./data.MUA_L4_peakFast(KORD&P9P13)),'omitnan'),sem(log2(data.MUA_L4_peakFast_K(KORD&P9P13)./data.MUA_L4_peakFast(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*6,log2(data.MUA_L4_peakFast_K(KORD&P14P18)./data.MUA_L4_peakFast(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(6.2,mean(log2(data.MUA_L4_peakFast_K(KORD&P14P18)./data.MUA_L4_peakFast(KORD&P14P18)),'omitnan'),sem(log2(data.MUA_L4_peakFast_K(KORD&P14P18)./data.MUA_L4_peakFast(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P9P13,:)),1)*4,log2(data.MUA_L4_peakFast_K(Saline&P9P13)./data.MUA_L4_peakFast(Saline&P9P13)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.MUA_L4_peakFast_K(Saline&P9P13)./data.MUA_L4_peakFast(Saline&P9P13)),'omitnan'),sem(log2(data.MUA_L4_peakFast_K(Saline&P9P13)./data.MUA_L4_peakFast(Saline&P9P13))),'go','Linewidth',1,'CapSize',10)
plot([-1,8],[0,0],'k--','LineWidth',0.75)
ax(1).YAxis.Label.String='MUA L4 fast peak change (ratio log2)';

ax(2)=subplot(3,1,2);
hold on
plot(ones(height(data(KORD&P5P8,:)),1)*0,log2(data.MUA_L4_peakSlow_K(KORD&P5P8)./data.MUA_L4_peakSlow(KORD&P5P8)),'o','Color',[150,150,150]/255)
errorbar(0.2,mean(log2(data.MUA_L4_peakSlow_K(KORD&P5P8)./data.MUA_L4_peakSlow(KORD&P5P8)),'omitnan'),sem(log2(data.MUA_L4_peakSlow_K(KORD&P5P8)./data.MUA_L4_peakSlow(KORD&P5P8))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P9P13,:)),1)*3,log2(data.MUA_L4_peakSlow_K(KORD&P9P13)./data.MUA_L4_peakSlow(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(3.2,mean(log2(data.MUA_L4_peakSlow_K(KORD&P9P13)./data.MUA_L4_peakSlow(KORD&P9P13)),'omitnan'),sem(log2(data.MUA_L4_peakSlow_K(KORD&P9P13)./data.MUA_L4_peakSlow(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*6,log2(data.MUA_L4_peakSlow_K(KORD&P14P18)./data.MUA_L4_peakSlow(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(6.2,mean(log2(data.MUA_L4_peakSlow_K(KORD&P14P18)./data.MUA_L4_peakSlow(KORD&P14P18)),'omitnan'),sem(log2(data.MUA_L4_peakSlow_K(KORD&P14P18)./data.MUA_L4_peakSlow(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P9P13,:)),1)*4,log2(data.MUA_L4_peakSlow_K(Saline&P9P13)./data.MUA_L4_peakSlow(Saline&P9P13)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.MUA_L4_peakSlow_K(Saline&P9P13)./data.MUA_L4_peakSlow(Saline&P9P13)),'omitnan'),sem(log2(data.MUA_L4_peakSlow_K(Saline&P9P13)./data.MUA_L4_peakSlow(Saline&P9P13))),'go','Linewidth',1,'CapSize',10)
plot([-1,8],[0,0],'k--','LineWidth',0.75)
ax(2).YAxis.Label.String='MUA L4 slow peak change (ratio log2)';

ax(3)=subplot(3,1,3);
hold on
plot(ones(height(data(KORD&P5P8,:)),1)*0,log2(data.MUA_L4_baselineFiring_K(KORD&P5P8)./data.MUA_L4_baselineFiring(KORD&P5P8)),'o','Color',[150,150,150]/255)
errorbar(0.2,mean(log2(data.MUA_L4_baselineFiring_K(KORD&P5P8)./data.MUA_L4_baselineFiring(KORD&P5P8)),'omitnan'),sem(log2(data.MUA_L4_baselineFiring_K(KORD&P5P8)./data.MUA_L4_baselineFiring(KORD&P5P8))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P9P13,:)),1)*3,log2(data.MUA_L4_baselineFiring_K(KORD&P9P13)./data.MUA_L4_baselineFiring(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(3.2,mean(log2(data.MUA_L4_baselineFiring_K(KORD&P9P13)./data.MUA_L4_baselineFiring(KORD&P9P13)),'omitnan'),sem(log2(data.MUA_L4_baselineFiring_K(KORD&P9P13)./data.MUA_L4_baselineFiring(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*6,log2(data.MUA_L4_baselineFiring_K(KORD&P14P18)./data.MUA_L4_baselineFiring(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(6.2,mean(log2(data.MUA_L4_baselineFiring_K(KORD&P14P18)./data.MUA_L4_baselineFiring(KORD&P14P18)),'omitnan'),sem(log2(data.MUA_L4_baselineFiring_K(KORD&P14P18)./data.MUA_L4_baselineFiring(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P9P13,:)),1)*4,log2(data.MUA_L4_baselineFiring_K(Saline&P9P13)./data.MUA_L4_baselineFiring(Saline&P9P13)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.MUA_L4_baselineFiring_K(Saline&P9P13)./data.MUA_L4_baselineFiring(Saline&P9P13)),'omitnan'),sem(log2(data.MUA_L4_baselineFiring_K(Saline&P9P13)./data.MUA_L4_baselineFiring(Saline&P9P13))),'go','Linewidth',1,'CapSize',10)
plot([-1,8],[0,0],'k--','LineWidth',0.75)
ax(3).YAxis.Label.String='MUA L4 baseline change (ratio log2)';

for i=1:numel(ax)
    ax(i).XLim=[-0.5,7.5];
    ax(i).FontSize=10;
    ax(i).YLim=[-3,3];
    ax(i).XTick=[0.5,3.5,6.5];
end

print(gcf, '-dpdf', fullfile(folderFigures,'4.19','MUA_log2Ratio'))
close 

figure('units','normalized','outerposition',[0 0 0.15 1]);
ax(1)=subplot(3,1,1);
hold on
plot(ones(height(data(KORD&P5P8,:)),1)*0,log2(data.PPR_K(KORD&P5P8)./data.PPR(KORD&P5P8)),'o','Color',[150,150,150]/255)
errorbar(0.2,mean(log2(data.PPC_K(KORD&P5P8)./data.PPR(KORD&P5P8)),'omitnan'),sem(log2(data.PPR_K(KORD&P5P8)./data.PPR(KORD&P5P8))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P9P13,:)),1)*3,log2(data.PPR_K(KORD&P9P13)./data.PPR(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(3.2,mean(log2(data.PPC_K(KORD&P9P13)./data.PPR(KORD&P9P13)),'omitnan'),sem(log2(data.PPR_K(KORD&P9P13)./data.PPR(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*6,log2(data.PPR_K(KORD&P14P18)./data.PPR(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(6.2,mean(log2(data.PPC_K(KORD&P14P18)./data.PPR(KORD&P14P18)),'omitnan'),sem(log2(data.PPR_K(KORD&P14P18)./data.PPR(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P9P13,:)),1)*4,log2(data.PPR_K(Saline&P9P13)./data.PPR(Saline&P9P13)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.PPC_K(Saline&P9P13)./data.PPR(Saline&P9P13)),'omitnan'),sem(log2(data.PPR_K(Saline&P9P13)./data.PPR(Saline&P9P13))),'go','Linewidth',1,'CapSize',10)
plot([-1,8],[0,0],'k--','LineWidth',0.75)
ax(1).YAxis.Label.String='PPR (ratio log2)';

ax(2)=subplot(3,1,2);
hold on
plot(ones(height(data(KORD&P5P8,:)),1)*0,log2(data.LFP_L4_peak_K(KORD&P5P8)./data.LFP_L4_peak(KORD&P5P8)),'o','Color',[150,150,150]/255)
errorbar(0.2,mean(log2(data.LFP_L4_peak_K(KORD&P5P8)./data.LFP_L4_peak(KORD&P5P8)),'omitnan'),sem(log2(data.LFP_L4_peak_K(KORD&P5P8)./data.LFP_L4_peak(KORD&P5P8))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P9P13,:)),1)*3,log2(data.LFP_L4_peak_K(KORD&P9P13)./data.LFP_L4_peak(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(3.2,mean(log2(data.LFP_L4_peak_K(KORD&P9P13)./data.LFP_L4_peak(KORD&P9P13)),'omitnan'),sem(log2(data.LFP_L4_peak_K(KORD&P9P13)./data.LFP_L4_peak(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*6,log2(data.LFP_L4_peak_K(KORD&P14P18)./data.LFP_L4_peak(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(6.2,mean(log2(data.LFP_L4_peak_K(KORD&P14P18)./data.LFP_L4_peak(KORD&P14P18)),'omitnan'),sem(log2(data.LFP_L4_peak_K(KORD&P14P18)./data.LFP_L4_peak(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P9P13,:)),1)*4,log2(data.LFP_L4_peak_K(Saline&P9P13)./data.LFP_L4_peak(Saline&P9P13)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.LFP_L4_peak_K(Saline&P9P13)./data.LFP_L4_peak(Saline&P9P13)),'omitnan'),sem(log2(data.LFP_L4_peak_K(Saline&P9P13)./data.LFP_L4_peak(Saline&P9P13))),'go','Linewidth',1,'CapSize',10)
plot([-1,8],[0,0],'k--','LineWidth',0.75)
ax(2).YAxis.Label.String='LFP L4 peak change (ratio log2)';

ax(3)=subplot(3,1,3);
hold on
plot(ones(height(data(KORD&P5P8,:)),1)*0,log2(data.LFP_L4_latency_K(KORD&P5P8)./data.LFP_L4_latency(KORD&P5P8)),'o','Color',[150,150,150]/255)
errorbar(0.2,mean(log2(data.LFP_L4_latency_K(KORD&P5P8)./data.LFP_L4_latency(KORD&P5P8)),'omitnan'),sem(log2(data.LFP_L4_latency_K(KORD&P5P8)./data.LFP_L4_latency(KORD&P5P8))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P9P13,:)),1)*3,log2(data.LFP_L4_latency_K(KORD&P9P13)./data.LFP_L4_latency(KORD&P9P13)),'o','Color',[150,150,150]/255)
errorbar(3.2,mean(log2(data.LFP_L4_latency_K(KORD&P9P13)./data.LFP_L4_latency(KORD&P9P13)),'omitnan'),sem(log2(data.LFP_L4_latency_K(KORD&P9P13)./data.LFP_L4_latency(KORD&P9P13))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(KORD&P14P18,:)),1)*6,log2(data.LFP_L4_latency_K(KORD&P14P18)./data.LFP_L4_latency(KORD&P14P18)),'o','Color',[150,150,150]/255)
errorbar(6.2,mean(log2(data.LFP_L4_latency_K(KORD&P14P18)./data.LFP_L4_latency(KORD&P14P18)),'omitnan'),sem(log2(data.LFP_L4_latency_K(KORD&P14P18)./data.LFP_L4_latency(KORD&P14P18))),'ko','Linewidth',1,'CapSize',10)
plot(ones(height(data(Saline&P9P13,:)),1)*4,log2(data.LFP_L4_latency_K(Saline&P9P13)./data.LFP_L4_latency(Saline&P9P13)),'o','Color',[150,150,150]/255)
errorbar(4.2,mean(log2(data.LFP_L4_latency_K(Saline&P9P13)./data.LFP_L4_latency(Saline&P9P13)),'omitnan'),sem(log2(data.LFP_L4_latency_K(Saline&P9P13)./data.LFP_L4_latency(Saline&P9P13))),'go','Linewidth',1,'CapSize',10)
plot([-1,8],[0,0],'k--','LineWidth',0.75)
ax(3).YAxis.Label.String='LFP L4 latency  change (ratio log2)';

for i=1:numel(ax)
    ax(i).XLim=[-0.5,7.5];
    ax(i).FontSize=10;
    ax(i).YLim=[-3,3];
    ax(i).XTick=[0.5,3.5,6.5];
end
print(gcf, '-dpdf', fullfile(folderFigures,'4.19','LFP_PPC_log2Ratio'))
close 

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(3,2,1)
hold on
plot(PSTHbins, data.MUA_L4(data.MouseID=='K51',:),'k')
plot([0,0],[0,500],'r')
plot([500 500],[0,500],'r')

title('Control - P5-P8')
ylim([0 200])
subplot(3,2,2)
hold on
plot(PSTHbins, data.MUA_L4_K(data.MouseID=='K51',:),'k')
plot([0,0],[0,500],'r')
plot([500 500],[0,500],'r')
title('KORD - P5-P8')
ylim([0 200])

subplot(3,2,3)
hold on
plot(PSTHbins, data.MUA_L4(data.MouseID=='K55',:),'k')
plot([0,0],[0,500],'r')
plot([500 500],[0,500],'r')
title('Control - P9-P13')
ylim([0 450])
subplot(3,2,4)
hold on
plot(PSTHbins, data.MUA_L4_K(data.MouseID=='K55',:),'k')
plot([0,0],[0,500],'r')
plot([500 500],[0,500],'r')
title('KORD - P9-P13')
ylim([0 450])

export_fig(fullfile(folderFigures,'4.19','MUA_L4_traces'),'-pdf','-transparent','-nocrop')
close 
% figure('units','normalized','outerposition',[0 0 0.2 1]);
% ax(1)=subplot(3,2,1);
% hold on
% plot([-5,-4],[data.PPC(KORD&P5P8),data.PPC_K(KORD&P5P8)],'-','Color',[150,150,150]/255)
% errorbar([-5,-4],[mean(data.PPC(KORD&P5P8),'omitnan'),mean(data.PPC_K(KORD&P5P8),'omitnan')],[sem(data.PPC(KORD&P5P8)),sem(data.PPC_K(KORD&P5P8))],'k-o','Linewidth',1,'CapSize',10)
% plot([-3,-2],[data.PPC(KORD&P9P13),data.PPC_K(KORD&P9P13)],'-','Color',[150,150,150]/255)
% errorbar([-3,-2],[mean(data.PPC(KORD&P9P13),'omitnan'),mean(data.PPC_K(KORD&P9P13),'omitnan')],[sem(data.PPC(KORD&P9P13)),sem(data.PPC_K(KORD&P9P13))],'k-o','Linewidth',1,'CapSize',10)
% plot([-1,0],[data.PPC(KORD&P14P18),data.PPC_K(KORD&P14P18)],'-','Color',[150,150,150]/255)
% errorbar([-1,0],[mean(data.PPC(KORD&P14P18),'omitnan'),mean(data.PPC_K(KORD&P14P18),'omitnan')],[sem(data.PPC(KORD&P14P18)),sem(data.PPC_K(KORD&P14P18))],'k-o','Linewidth',1,'CapSize',10)
% ax(1).YAxis.Label.String='PPC';
% ax(1).Title.String = 'KORD';

%% Awake
P9=data.mouseAge<14;
aw=data.state=='Awake';
ur=data.state=='Urethane';

figure('units','normalized','outerposition',[0 0 0.2 1])

ax(1)=subplot(3,1,1);
v(1,:)=violinplot(maxBetaPSD(P9), data.state(P9));
ax(1).YLabel.String='Max \beta normalised power';
[h,p,stats] = my_ttest2 (maxBetaPSD(P9&aw),maxBetaPSD(P9&ur))


ax(2)=subplot(3,1,2);
v(2,:)=violinplot(data.SB_L4_duration(P9), data.state(P9));
ax(2).YLabel.String='Spindle burst duration (s)';
[h,p,stats] = my_ttest2 (data.SB_L4_duration(P9&aw),data.SB_L4_duration(P9&ur))

ax(3)=subplot(3,1,3);
v(3,:)=violinplot(data.SB_L4_frequency(P9), data.state(P9));
ax(3).YLabel.String='Spindle burst occurrence (Hz)';
[h,p,stats] = my_ttest2 (data.SB_L4_frequency(P9&aw),data.SB_L4_frequency(P9&ur))

color = [200,200,200;20,20,20];

for i=1:numel(ax)
    for j=1:size(v,2)
        v(i,j).ViolinColor=color(j,:)/255;
        v(i,j).ScatterPlot.MarkerFaceColor=[37,37,37]/255;
        v(i,j).ScatterPlot.MarkerFaceAlpha=1;
        v(i,j).ScatterPlot.SizeData=10;
    end

    ax(i).FontSize=12;
    ax(i).LineWidth=1;
    ax(i).XLim=[0.5,2.5];
end

print(gcf, '-dpdf', fullfile(folderFigures,'4.21','PSD_SB_violin'))
% close

figure('units','normalized','outerposition',[0 0 0.2 1])

ax(1)=subplot(3,1,1);
v(1,:)=violinplot(data.MUA_L4_peakFast(P9), data.state(P9));
ax(1).YLabel.String='Max L4 MUA fast (spikes/s)';
[h,p,stats] = my_ttest2 (data.MUA_L4_peakFast(P9&aw),data.MUA_L4_peakFast(P9&ur))


ax(2)=subplot(3,1,2);
v(2,:)=violinplot(data.MUA_L4_baselineFiring(P9), data.state(P9));
ax(2).YLabel.String='Spontaneous L4 MUA (spike/s)';
[h,p,stats] = my_ttest2 (data.MUA_L4_baselineFiring(P9&aw),data.MUA_L4_baselineFiring(P9&ur))

ax(3)=subplot(3,1,3);
v(3,:)=violinplot(data.MUA_L4_peakSlow(P9), data.state(P9));
ax(3).YLabel.String='Max L4 MUA slow (spikes/s)';
[h,p,stats] = my_ttest2 (data.MUA_L4_peakSlow(P9&aw),data.MUA_L4_peakSlow(P9&ur))

color = [200,200,200;20,20,20];

for i=1:numel(ax)
    for j=1:size(v,2)
        v(i,j).ViolinColor=color(j,:)/255;
        v(i,j).ScatterPlot.MarkerFaceColor=[37,37,37]/255;
        v(i,j).ScatterPlot.MarkerFaceAlpha=1;
        v(i,j).ScatterPlot.SizeData=10;
    end

    ax(i).FontSize=12;
    ax(i).LineWidth=1;
    ax(i).XLim=[0.5,2.5];
end

print(gcf, '-dpdf', fullfile(folderFigures,'4.21','MUA_violin'))
% close

figure
hold on
plot(mean(data.MUA_L4(ur,:)),'k')
plot(mean(data.MUA_L4(aw,:)),'r')
