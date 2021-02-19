close all
clear 
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 

recordings=readtable('V1_InVivo_SST;Ai32.csv');
recID=unique(recordings.MouseID);
recFolder='D:\InVivo_SST';

mouseID=[];
suid=[];
suage=[];
sulayer=[];
wf=[];
halfWidth=[];
troughPeakTime=[];
peakTroughRatio=[];
endSlope=[];
PSTHvisual=[];
PSTHoptotagging=[];
PSTHvisualOpto=[];
PSTHlaser=[];
responseVisual=[];
responseTag=[];
responseVisualOpto=[];
responseLaser=[];

PPC=[];

for i=1:numel(recID)
    if i~=7
        dir=fullfile(recFolder,recID{i});

        load(fullfile(dir,'spikes.mat'),'s')
        load(fullfile(dir,'rez.mat'),'rez')

        thisID=cell(numel(s.suid),1);
        thisID(:)={recID{i}};

        mouseID=[mouseID;thisID];
        suid=[suid;s.suid];
        suage=[suage;s.suage];
        sulayer=[sulayer;s.sulayer];
        halfWidth=[halfWidth;s.halfWidth];
        troughPeakTime=[troughPeakTime;s.troughPeakTime];
        peakTroughRatio=[peakTroughRatio;s.peakTroughRatio];
        endSlope=[endSlope;s.endSlope];
        wf=[wf;s.suWf];
        PSTHvisual=[PSTHvisual;s.PSTHvisual];
        PSTHoptotagging=[PSTHoptotagging;s.PSTHoptotagging];
        PSTHvisualOpto=[PSTHvisualOpto;s.PSTHvisualOpto];
        if isempty(s.PSTHlaser)
            PSTHlaser=[PSTHlaser;nan(numel(s.suid),numel(s.PSTHbins))];
            responseLaser=[responseLaser;nan(numel(s.suid),size(responseLaser,2))];
        else
            PSTHlaser=[PSTHlaser;s.PSTHlaser];
            responseLaser=[responseLaser;s.response.laser];
        end
        responseVisual=[responseVisual;s.response.visual];
        responseTag=[responseTag;s.response.optotagging];
        responseVisualOpto=[responseVisualOpto;s.response.visualOpto];
    end
    PPC=[PPC;rez.ops.PPC];
end

%remove units for which don't have laser only condition
PSTHvisual=PSTHvisual(~isnan(PSTHlaser(:,1)),:);
PSTHvisualOpto=PSTHvisualOpto(~isnan(PSTHlaser(:,1)),:);
PSTHlaser=PSTHlaser(~isnan(PSTHlaser(:,1)),:);

%Z-scoring SU
Z_PSTHvisual=zscoreBaseline(PSTHvisual);
Z_PSTHvisualOpto=zscoreBaseline(PSTHvisualOpto);
Z_PSTHlaser=zscoreBaseline(PSTHlaser);
Z_PSTHoptotagging=zscoreBaseline(PSTHoptotagging);

% Make table
t=table;
t.MouseID=mouseID;
t.suid=suid;
t.SST=responseTag;

PSTHbins=s.PSTHbins;

%% Set group logic arrays

%Age
P5P8=suage<9;
P9P13=suage>=9 & suage<14;
P9P11=suage>=9 & suage<12;
P12P13=suage>=12 & suage<14;
P14P18=suage>=14;

%Layer
L23=(sulayer==1);
L4=(sulayer==2);
L56=(sulayer==3);

%Visual response
ON=responseVisual(:,1)==1 & responseVisual(:,2)==0 & responseVisual(:,3)==0;
OFF=responseVisual(:,2)==1 & responseVisual(:,1)==0 & responseVisual(:,3)==0;
late=responseVisual(:,3)==1 & responseVisual(:,2)==0 & responseVisual(:,1)==0;
nonResponsive=responseVisual(:,3)==0 & responseVisual(:,2)==0 & responseVisual(:,1)==0;
ONOFF=responseVisual(:,1)==1 & responseVisual(:,2)==1 & responseVisual(:,3)==0;
ONlate=responseVisual(:,1)==1 & responseVisual(:,2)==0 & responseVisual(:,3)==1;
OFFlate=responseVisual(:,1)==0 & responseVisual(:,2)==1 & responseVisual(:,3)==1;
ONOFFlate=responseVisual(:,1)==1 & responseVisual(:,2)==1 & responseVisual(:,3)==1;

responsive=~nonResponsive;

SST=responseTag(:,1)==1;

SST_E=responseLaser(:,1)==1;
SST_I=responseLaser(:,1)==2;
SST_Nonresponsive=responseLaser(:,1)==0;


FS=troughPeakTime<=0.7;
SST=(SST & ~FS);
RS=((~FS)&(~SST));

%% PSTH plotting - responsive cells
figure('units','normalized','outerposition',[0 0 1 1]);

ax=subplot(3,3,1);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P5P8&L23,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P5P8&L23,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P5P8&L23,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
% ax.YLim=[0 max(mean(PSTHvisual(responsive&FS&P5P8&L23,:),1))+1];
ax.XLim=[-200,5000];
ax.Title.String='P5-P8';

ax=subplot(3,3,2);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P9P13&L23,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P9P13&L23,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P9P13&L23,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(PSTHvisual(responsive&FS&P9P13&L23,:),1))+1];
ax.XLim=[-200,5000];
ax.Title.String='P9-P13';

ax=subplot(3,3,3);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P14P18&L23,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P14P18&L23,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P14P18&L23,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(PSTHvisual(responsive&FS&P14P18&L23,:),1))+1];
ax.XLim=[-200,5000];
ax.Title.String='P14-P18';


ax=subplot(3,3,4);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P5P8&L4,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P5P8&L4,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P5P8&L4,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
% ax.YLim=[0 max(mean(PSTHvisual(responsive&FS&P5P8&L4,:),1))+1];
ax.XLim=[-200,5000];

ax=subplot(3,3,5);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P9P13&L4,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P9P13&L4,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P9P13&L4,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(PSTHvisual(responsive&FS&P9P13&L4,:),1))+1];
ax.XLim=[-200,5000];

ax=subplot(3,3,6);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P14P18&L4,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P14P18&L4,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P14P18&L4,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(PSTHvisual(responsive&FS&P14P18&L4,:),1))+1];
ax.XLim=[-200,5000];

ax=subplot(3,3,7);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P5P8&L56,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P5P8&L56,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P5P8&L56,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
% ax.YLim=[0 max(mean(PSTHvisual(responsive&FS&P5P8&L56,:),1))+1];
ax.XLim=[-200,5000];

ax=subplot(3,3,8);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P9P13&L56,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P9P13&L56,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P9P13&L56,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(PSTHvisual(responsive&FS&P9P13&L56,:),1))+1];
ax.XLim=[-200,5000];

ax=subplot(3,3,9);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P14P18&L56,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P14P18&L56,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P14P18&L56,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(PSTHvisual(responsive&FS&P14P18&L56,:),1))+1];
ax.XLim=[-200,5000];

legend('RS','FS','SST')
sgtitle('Responsive cells')

%% PSTH plotting - all cells
figure('units','normalized','outerposition',[0 0 1 1]);

ax=subplot(3,3,1);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(Z_PSTHvisual(RS&P5P8&L23,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(FS&P5P8&L23,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(SST&P5P8&L23,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(Z_PSTHvisual(RS&P5P8&L23,:),1))+1];
ax.XLim=[-200,5000];
ax.Title.String='P5-P8';

ax=subplot(3,3,2);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(Z_PSTHvisual(RS&P9P13&L23,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(FS&P9P13&L23,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(SST&P9P13&L23,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(Z_PSTHvisual(FS&P9P13&L23,:),1))+1];
ax.XLim=[-200,5000];
ax.Title.String='P9-P13';


ax=subplot(3,3,3);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(Z_PSTHvisual(RS&P14P18&L23,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(FS&P14P18&L23,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(SST&P14P18&L23,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(Z_PSTHvisual(FS&P14P18&L23,:),1))+1];
ax.XLim=[-200,5000];
ax.Title.String='P14-P18';


ax=subplot(3,3,4);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(Z_PSTHvisual(RS&P5P8&L4,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(FS&P5P8&L4,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(SST&P5P8&L4,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(Z_PSTHvisual(RS&P5P8&L4,:),1))+1];
ax.XLim=[-200,5000];

ax=subplot(3,3,5);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(Z_PSTHvisual(RS&P9P13&L4,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(FS&P9P13&L4,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(SST&P9P13&L4,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(Z_PSTHvisual(FS&P9P13&L4,:),1))+1];
ax.XLim=[-200,5000];

ax=subplot(3,3,6);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(Z_PSTHvisual(RS&P14P18&L56,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(FS&P14P18&L56,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(SST&P14P18&L56,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(Z_PSTHvisual(FS&P14P18&L56,:),1))+1];
ax.XLim=[-200,5000];

ax=subplot(3,3,7);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(Z_PSTHvisual(RS&P5P8&L56,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(FS&P5P8&L56,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(SST&P5P8&L56,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(Z_PSTHvisual(RS&P5P8&L56,:),1))+1];
ax.XLim=[-200,5000];

ax=subplot(3,3,8);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(Z_PSTHvisual(RS&P9P13&L56,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(FS&P9P13&L56,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(SST&P9P13&L56,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(Z_PSTHvisual(FS&P9P13&L56,:),1))+1];
ax.XLim=[-200,5000];

ax=subplot(3,3,9);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(Z_PSTHvisual(RS&P14P18&L56,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(FS&P14P18&L56,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(Z_PSTHvisual(SST&P14P18&L56,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(Z_PSTHvisual(FS&P14P18&L56,:),1))+1];
ax.XLim=[-200,5000];

legend('RS','FS','SST')
sgtitle('All cells')


%% PSTH plotting - responsive cells - 4 Dev windows
figure('units','normalized','outerposition',[0 0 1 1]);

ax=subplot(3,4,1);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P5P8&L23,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P5P8&L23,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P5P8&L23,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 12];
ax.XLim=[-200,5000];
ax.Title.String='P5-P8';

ax=subplot(3,4,2);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P9P11&L23,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P9P11&L23,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P9P11&L23,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 12];
ax.XLim=[-200,5000];
ax.Title.String='P9-P11';

ax=subplot(3,4,3);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P12P13&L23,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P12P13&L23,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P12P13&L23,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 15];
ax.XLim=[-200,5000];
ax.Title.String='P12-P13';

ax=subplot(3,4,4);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P14P18&L23,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P14P18&L23,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P14P18&L23,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(PSTHvisual(responsive&FS&P14P18&L23,:),1))+1];
ax.XLim=[-200,5000];
ax.Title.String='P14-P18';

ax=subplot(3,4,5);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P5P8&L4,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P5P8&L4,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P5P8&L4,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 12];
ax.XLim=[-200,5000];

ax=subplot(3,4,6);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P9P11&L4,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P9P11&L4,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P9P11&L4,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 12];
ax.XLim=[-200,5000];

ax=subplot(3,4,7);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P12P13&L4,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P12P13&L4,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P12P13&L4,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 25];
ax.XLim=[-200,5000];

ax=subplot(3,4,8);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P14P18&L4,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P14P18&L4,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P14P18&L4,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(PSTHvisual(responsive&FS&P14P18&L4,:),1))+1];
ax.XLim=[-200,5000];

ax=subplot(3,4,9);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P5P8&L56,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P5P8&L56,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P5P8&L56,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 12];
ax.XLim=[-200,5000];

ax=subplot(3,4,10);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P9P11&L56,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P9P11&L56,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P9P11&L56,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 12];
ax.XLim=[-200,5000];

ax=subplot(3,4,11);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P12P13&L56,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P12P13&L56,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P12P13&L56,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 35];
ax.XLim=[-200,5000];

ax=subplot(3,4,12);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(responsive&RS&P14P18&L56,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&FS&P14P18&L56,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(responsive&SST&P14P18&L56,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(PSTHvisual(responsive&FS&P14P18&L56,:),1))+1];
ax.XLim=[-200,5000];

legend('RS','FS','SST')
sgtitle('Responsive cells')

%% PSTH plotting - responsive cells - 4 Dev windows
figure('units','normalized','outerposition',[0 0 1 1]);

ax=subplot(3,4,1);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(RS&P5P8&L23,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(FS&P5P8&L23,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(SST&P5P8&L23,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 12];
ax.XLim=[-200,5000];
ax.Title.String='P5-P8';

ax=subplot(3,4,2);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(RS&P9P11&L23,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(FS&P9P11&L23,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(SST&P9P11&L23,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 12];
ax.XLim=[-200,5000];
ax.Title.String='P9-P11';

ax=subplot(3,4,3);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(RS&P12P13&L23,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(FS&P12P13&L23,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(SST&P12P13&L23,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 15];
ax.XLim=[-200,5000];
ax.Title.String='P12-P13';

ax=subplot(3,4,4);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(RS&P14P18&L23,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(FS&P14P18&L23,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(SST&P14P18&L23,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(PSTHvisual(FS&P14P18&L23,:),1))+1];
ax.XLim=[-200,5000];
ax.Title.String='P14-P18';

ax=subplot(3,4,5);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(RS&P5P8&L4,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(FS&P5P8&L4,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(SST&P5P8&L4,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 12];
ax.XLim=[-200,5000];

ax=subplot(3,4,6);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(RS&P9P11&L4,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(FS&P9P11&L4,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(SST&P9P11&L4,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 12];
ax.XLim=[-200,5000];

ax=subplot(3,4,7);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(RS&P12P13&L4,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(FS&P12P13&L4,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(SST&P12P13&L4,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 25];
ax.XLim=[-200,5000];

ax=subplot(3,4,8);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(RS&P14P18&L4,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(FS&P14P18&L4,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(SST&P14P18&L4,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(PSTHvisual(FS&P14P18&L4,:),1))+1];
ax.XLim=[-200,5000];

ax=subplot(3,4,9);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(RS&P5P8&L56,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(FS&P5P8&L56,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(SST&P5P8&L56,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 12];
ax.XLim=[-200,5000];

ax=subplot(3,4,10);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(RS&P9P11&L56,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(FS&P9P11&L56,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(SST&P9P11&L56,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 12];
ax.XLim=[-200,5000];

ax=subplot(3,4,11);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(RS&P12P13&L56,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(FS&P12P13&L56,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(SST&P12P13&L56,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 35];
ax.XLim=[-200,5000];

ax=subplot(3,4,12);
hold on;
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
plot(PSTHbins,mean(PSTHvisual(RS&P14P18&L56,:),1),'k','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(FS&P14P18&L56,:),1),'r','LineWidth',2)
plot(PSTHbins,mean(PSTHvisual(SST&P14P18&L56,:),1),'b','LineWidth',2)
ax=applyFont(ax,0);
ax.YLim=[0 max(mean(PSTHvisual(FS&P14P18&L56,:),1))+1];
ax.XLim=[-200,5000];

legend('RS','FS','SST')
sgtitle('Responsive cells')