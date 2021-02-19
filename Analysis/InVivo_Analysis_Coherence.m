close all
clear 
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 

recordings=readtable('V1_InVivo_SST;Ai32.csv');
recID=unique(recordings.MouseID);
recFolder='D:\InVivo_SST';
dirOUT='C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures\Coherence and phaseLocking';

mouseID=[];
suid=[];
suage=[];
sulayer=[];
responseTag=[];
troughPeakTime=[];
PPC=[];
vectorLength=[];
vectorAngle=[];
pValuePPC=[];

zscorecoherence=0;
coherence=[];

for i=1:numel(recID)
        dir=fullfile(recFolder,recID{i});
        load(fullfile(dir,'rez.mat'),'rez')
        load(fullfile(dir,'spikes.mat'),'s')
        thisID=cell(numel(s.suid),1);
        
        thisID(:)={recID{i}};
        
        mouseID=[mouseID;thisID];
        suid=[suid;s.suid];
        suage=[suage;s.suage];
        sulayer=[sulayer;s.sulayer];
        responseTag=[responseTag;s.response.optotagging];
        troughPeakTime=[troughPeakTime;s.troughPeakTime];

        PPC=[PPC;rez.ops.phaseLocking.PPC];
        vectorLength=[vectorLength;rez.ops.phaseLocking.vectorLength];
        vectorAngle=[vectorAngle;rez.ops.phaseLocking.vectorAngle];
        pValuePPC=[pValuePPC;rez.ops.phaseLocking.pValue];

        
        [betaCoherence,maxCoherenceIdx]=max(rez.ops.coherence.curve(:,(rez.ops.coherence.freqs>=10) & (rez.ops.coherence.freqs<=30)),[],2);
        [lowCoherence,~]=max(rez.ops.coherence.curve(:,(rez.ops.coherence.freqs>=1) & (rez.ops.coherence.freqs<=5)),[],2);
        [gammaCoherence,~]=max(rez.ops.coherence.curve(:,(rez.ops.coherence.freqs>=50) & (rez.ops.coherence.freqs<=80)),[],2);

        if zscorecoherence
            betaCoherence=(betaCoherence-rez.ops.coherence.curve_shuff_mean(maxCoherenceIdx))./rez.ops.coherence.curve_suff_SD(maxCoherenceIdx);
        end
        coherence=[coherence;lowCoherence,betaCoherence,gammaCoherence];
end
coherence(coherence==0)=NaN;
coherenceFreq=rez.ops.coherence.freqs;

%% Set group logic arrays

%Age
P5P8=suage<9;
P9P13=suage>=9 & suage<14;
P14P18=suage>=14;
SST=responseTag(:,1)==1;
FS=troughPeakTime<=0.7;
FS=(FS&~SST);
RS=(~SST)&(~FS);

dev(P5P8,1)=categorical(cellstr('P5-P8'));
dev(P9P13,1)=categorical(cellstr('P9-P13'));
dev(P14P18,1)=categorical(cellstr('P14-P18'));

frequencyBands={'Low (1-5 Hz)','Beta (10-30 Hz)','Gamma (40-80 Hz)'};

%%
coherencetab=[nanmean(coherence(P5P8,1)),nanstd(coherence(P5P8,1)),numel(coherence(P5P8,1)),nanmean(coherence(P9P13,1)),nanstd(coherence(P9P13,1)),numel(coherence(P9P13,1)),nanmean(coherence(P14P18,1)),nanstd(coherence(P14P18,1)),numel(coherence(P14P18,1));...
    nanmean(coherence(P5P8,2)),nanstd(coherence(P5P8,2)),numel(coherence(P5P8,2)),nanmean(coherence(P9P13,2)),nanstd(coherence(P9P13,2)),numel(coherence(P9P13,2)),nanmean(coherence(P14P18,2)),nanstd(coherence(P14P18,2)),numel(coherence(P14P18,2));...
    nanmean(coherence(P5P8,3)),nanstd(coherence(P5P8,3)),numel(coherence(P5P8,3)),nanmean(coherence(P9P13,3)),nanstd(coherence(P9P13,3)),numel(coherence(P9P13,3)),nanmean(coherence(P14P18,3)),nanstd(coherence(P14P18,3)),numel(coherence(P14P18,3))];


%% Standard coherence per development plot
figure('units','normalized','outerposition',[0 0 .5 .75]);
boxPlotY=[{coherence(:,1)},{coherence(:,2)},{coherence(:,3)}];
boxPlotX=[{dev},{dev},{dev}];
groupedBoxPlot(boxPlotY,boxPlotX,'ColorScale','Black','Orientation','vertical','Jitter',0.07);
ax=gca;
ax.FontName='Arial';
ax.YLabel.String='Coherence';
ax.YLim=[0 1];
text([1 2 3],[0 0 0]-0.05,{'P5-P8','P9-P13','P14-P18'},'FontSize',20)
legend(frequencyBands)
legend('boxoff')
export_fig(fullfile(dirOUT,'CoherenceAll'),'-pdf','-transparent','-nocrop')
close
%% SST vs pyr vs fs Coherence
figure('units','normalized','outerposition',[0 0 1 .5]);
for i=1:numel(frequencyBands)
    subplot(1,3,i)
    boxPlotY=[{coherence(SST,i)},{coherence(RS,i)},{coherence(FS,i)}];
    boxPlotX=[{dev(SST)},{dev(RS)},{dev(FS)}];
    groupedBoxPlot(boxPlotY,boxPlotX,'ColorScale','Rainbow','Orientation','vertical','Jitter',0.07);
    ax=gca;
    ax.Title.String=frequencyBands{i};
    ax.FontName='Arial';
    if i==1; ax.YLabel.String='Coherence'; end
    ax.YLim=[0 1];
    text([1 2 3],[0 0 0]-0.05,{'P5-P8','P9-P13','P14-P18'},'FontSize',20)
    if i==numel(frequencyBands)
        legend('SST','RS','FS')
        legend('boxoff')    
    end
end
export_fig(fullfile(dirOUT,'CoherenceSSTvsRSvsFS'),'-pdf','-transparent','-nocrop')
close
%% PPC all
figure('units','normalized','outerposition',[0 0 1 .5]);
for i=1:numel(frequencyBands)
    subplot(1,3,i)
    boxPlotY=[{PPC(SST,1)},{PPC(RS,1)},{PPC(FS,1)}];
    boxPlotX=[{dev(SST)},{dev(RS)},{dev(FS)}];
    groupedBoxPlot(boxPlotY,boxPlotX,'ColorScale','Rainbow','Orientation','vertical','Jitter',0.07);
    ax=gca;
    ax.Title.String=frequencyBands{i};
    ax.FontName='Arial';
    if i==1; ax.YLabel.String='PPC'; end
    ax.YLim=[-.1 .3];
    text([1 2 3],[0 0 0]-0.12,{'P5-P8','P9-P13','P14-P18'},'FontSize',20)
    if i==numel(frequencyBands)
        legend('SST','RS','FS')
        legend('boxoff')    
    end
end
export_fig(fullfile(dirOUT,'PPC_All_SSTvsRSvsFS'),'-pdf','-transparent','-nocrop')
close
%% Ratio of entrained
entrained=pValuePPC<0.05;
figure('units','normalized','outerposition',[0 0 1 .5]);
colors={[215,48,39]/255,[254,224,144]/255,[69,117,180]/255};

for i =1:numel(frequencyBands)
    subplot(1,3,i)
    entrainedPercentage=[nnz(entrained(P5P8&SST,i))/numel(entrained(P5P8&SST,i)), nnz(entrained(P5P8&RS,i))/numel(entrained(P5P8&RS,i)), nnz(entrained(P5P8&FS,i))/numel(entrained(P5P8&FS,i));...
        nnz(entrained(P9P13&SST,i))/numel(entrained(P9P13&SST,i)), nnz(entrained(P9P13&RS,i))/numel(entrained(P9P13&RS,i)), nnz(entrained(P9P13&FS,i))/numel(entrained(P9P13&FS,i));...
        nnz(entrained(P14P18&SST,i))/numel(entrained(P14P18&SST,i)), nnz(entrained(P14P18&RS,i))/numel(entrained(P14P18&RS,i)), nnz(entrained(P14P18&FS,i))/numel(entrained(P14P18&FS,i))];
    b=bar(entrainedPercentage);
    for j=1:size(b,2); b(j).FaceColor=colors{j}; end
    ax=gca;
    ax.FontName='Arial';
    ax.YLim=[0 1];
    ax.Title.String=frequencyBands{i};
    ax.Box='off';
    ax.FontSize=20;
    ax.LineWidth=1;
    ax.XAxis.Visible='off';
    if i==1; ax.YLabel.String='Entrained single units ratio'; end

    
    text([1 2 3]-0.3,[0 0 0]-0.05,{'P5-P8','P9-P13','P14-P18'},'FontSize',20)
    if i==numel(frequencyBands)
        legend('SST','RS','FS','Location','northwest')
        legend('boxoff')    
    end
end
export_fig(fullfile(dirOUT,'PPC_RatioEntrained'),'-pdf','-transparent','-nocrop')
close
%% PPC entrained only
figure('units','normalized','outerposition',[0 0 1 .5]);
for i =1:numel(frequencyBands)
    subplot(1,3,i)
    boxPlotY=[{PPC(entrained(:,i)&SST,i)},{PPC(entrained(:,i)&RS,i)},{PPC(entrained(:,i)&FS,i)}];
    boxPlotX=[{dev(entrained(:,i)&SST)},{dev(entrained(:,i)&RS)},{dev(entrained(:,i)&FS)}];
    groupedBoxPlot(boxPlotY,boxPlotX,'ColorScale','Rainbow','Orientation','vertical','Jitter',0.05);
    ax=gca;
    ax.FontName='Arial';
    ax.Title.String=frequencyBands{i};
    ax.XLabel.String='Total IPSC charge (pC)';
    ax.YLim=[-.01 .55];
    if i==1; ax.YLabel.String='PPC'; end
    text([1 2 3],[0 0 0]-0.02,{'P5-P8','P9-P13','P14-P18'},'FontSize',20)
    if i==numel(frequencyBands)
        legend('SST','RS','FS','Location','northeast')
        legend('boxoff')    
    end
end
export_fig(fullfile(dirOUT,'PPC_Entrained_SSTvsRSvsFS'),'-pdf','-transparent','-nocrop')
close