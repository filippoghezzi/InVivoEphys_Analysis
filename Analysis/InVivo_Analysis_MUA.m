clear 
close all
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 
recordings=readtable('V1_InVivo_SST;Ai32.csv');
recID=unique(recordings.MouseID);
recFolder='D:\InVivo_V1';

smoothingIndex=5;

for i=1:numel(recID)
    dir=fullfile(recFolder,recID{i});

    %%
    if i~=7
        verbose=1;
        load(fullfile(dir,'rez.mat'),'rez')
        load(fullfile(dir,'stimuli.mat'),'stim')
        ops=rez.ops;

        s=loadSpikes(dir);
        MUA_L4=getLayerMUA(s,stim.ledR,ops.fs,ops.LFPwindow,ops.PSTHbinSize,ops.L4);
        MUAvisualL4(i,:)=smooth(MUA_L4.raw,smoothingIndex);

        MUA_L2=getLayerMUA(s,stim.ledR,ops.fs,ops.LFPwindow,ops.PSTHbinSize,1:min(ops.L4)-1);
        MUAvisualL2(i,:)=smooth(MUA_L2.raw,smoothingIndex);

        MUA_L5=getLayerMUA(s,stim.ledR,ops.fs,ops.LFPwindow,ops.PSTHbinSize,max(ops.L4)+1:ops.NchanTOT);
        MUAvisualL5(i,:)=smooth(MUA_L5.raw,smoothingIndex);
        
        
        artefactRemoval=[-51,49;149,153]*10^-3; %s

        MUA_L4opto=getLayerMUA(s,stim.laserAndLedR,ops.fs,ops.LFPwindow,ops.PSTHbinSize,ops.L4,artefactRemoval);
        MUAoptoL4(i,:)=smooth(MUA_L4opto.raw,smoothingIndex);

        MUA_L2opto=getLayerMUA(s,stim.laserAndLedR,ops.fs,ops.LFPwindow,ops.PSTHbinSize,1:min(ops.L4)-1,artefactRemoval);
        MUAoptoL2(i,:)=smooth(MUA_L2opto.raw,smoothingIndex);

        MUA_L5opto=getLayerMUA(s,stim.laserAndLedR,ops.fs,ops.LFPwindow,ops.PSTHbinSize,max(ops.L4)+1:ops.NchanTOT,artefactRemoval);
        MUAoptoL5(i,:)=smooth(MUA_L5opto.raw,smoothingIndex);

        MUA_L4laser=getLayerMUA(s,stim.laser+0.05*ops.fs,ops.fs,ops.LFPwindow,ops.PSTHbinSize,ops.L4,artefactRemoval);
        MUAlaserL4(i,:)=smooth(MUA_L4laser.raw,smoothingIndex);

        MUA_L2laser=getLayerMUA(s,stim.laser+0.05*ops.fs,ops.fs,ops.LFPwindow,ops.PSTHbinSize,1:min(ops.L4)-1,artefactRemoval);
        MUAlaserL2(i,:)=smooth(MUA_L2laser.raw,smoothingIndex);

        MUA_L5laser=getLayerMUA(s,stim.laser+0.05*ops.fs,ops.fs,ops.LFPwindow,ops.PSTHbinSize,max(ops.L4)+1:ops.NchanTOT,artefactRemoval);
        MUAlaserL5(i,:)=smooth(MUA_L5laser.raw,smoothingIndex);

        
        
        
        

        if ops.age <9
            dev(i,1)=1;
        elseif ops.age>=9 && ops.age <14
            dev(i,1)=2;
        elseif ops.age>=14
            dev(i,1)=3;
        end
    end
end
bins=MUA_L4.bins;

%% Z-scoring by the baseline
Z_MUAvisualL2=zscoreBaseline(MUAvisualL2);
Z_MUAvisualL4=zscoreBaseline(MUAvisualL4);
Z_MUAvisualL5=zscoreBaseline(MUAvisualL5);
Z_MUAoptoL2=zscoreBaseline(MUAoptoL2);
Z_MUAoptoL4=zscoreBaseline(MUAoptoL4);
Z_MUAoptoL5=zscoreBaseline(MUAoptoL5);
Z_MUAlaserL2=zscoreBaseline(MUAlaserL2);
Z_MUAlaserL4=zscoreBaseline(MUAlaserL4);
Z_MUAlaserL5=zscoreBaseline(MUAlaserL5);

%% Compare development
figure('units','normalized','outerposition',[0 0 1 1]);
ax1=subplot(3,2,1);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none')
hold on; 
plot(bins,mean(MUAvisualL2(dev==1,:),1),'LineWidth',2)
plot(bins,mean(MUAvisualL2(dev==2,:),1),'LineWidth',2)
plot(bins,mean(MUAvisualL2(dev==3,:),1),'LineWidth',2)
legend('Light Stim','P5-P8','P9-P13','P14-P18')
ax1=applyFont(ax1,0);


ax2=subplot(3,2,3);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none') 
hold on; 
plot(bins,mean(MUAvisualL4(dev==1,:),1),'LineWidth',2)
plot(bins,mean(MUAvisualL4(dev==2,:),1),'LineWidth',2)
plot(bins,mean(MUAvisualL4(dev==3,:),1),'LineWidth',2)
legend('Light Stim','P5-P8','P9-P13','P14-P18')
ax2=applyFont(ax2,0);

ax3=subplot(3,2,5);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none') 
hold on; 
plot(bins,mean(MUAvisualL5(dev==1,:),1),'LineWidth',2)
plot(bins,mean(MUAvisualL5(dev==2,:),1),'LineWidth',2)
plot(bins,mean(MUAvisualL5(dev==3,:),1),'LineWidth',2)
legend('Light Stim','P5-P8','P9-P13','P14-P18')
ax3=applyFont(ax3,0);
ax3.XLabel.String='Time (ms)';


ax4=subplot(3,2,2);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none') 
hold on; 
plot(bins,mean(Z_MUAvisualL2(dev==1,:),1),'LineWidth',2)
plot(bins,mean(Z_MUAvisualL2(dev==2,:),1),'LineWidth',2)
plot(bins,mean(Z_MUAvisualL2(dev==3,:),1),'LineWidth',2)
legend('Light Stim','P5-P8','P9-P13','P14-P18')
ax4=applyFont(ax4,1);

ax5=subplot(3,2,4);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none') 
hold on; 
plot(bins,mean(Z_MUAvisualL4(dev==1,:),1),'LineWidth',2)
plot(bins,mean(Z_MUAvisualL4(dev==2,:),1),'LineWidth',2)
plot(bins,mean(Z_MUAvisualL4(dev==3,:),1),'LineWidth',2)
legend('Light Stim','P5-P8','P9-P13','P14-P18')
ax5=applyFont(ax5,1);

ax6=subplot(3,2,6);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none') 
hold on; 
plot(bins,mean(Z_MUAvisualL5(dev==1,:),1),'LineWidth',2)
plot(bins,mean(Z_MUAvisualL5(dev==2,:),1),'LineWidth',2)
plot(bins,mean(Z_MUAvisualL5(dev==3,:),1),'LineWidth',2)
legend('Light Stim','P5-P8','P9-P13','P14-P18')
ax6=applyFont(ax6,1);
ax6.XLabel.String='Time (ms)';

%save figure
figname='MUA_Layer_Dev';
% sg=sgtitle(strcat('P',int2str(ops.age)));
% sg.FontSize=30;
export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')




%% Compare visual vs. opto
figure('units','normalized','outerposition',[0 0 1 1]);
plotLaserOnly=1;

%%P5-P8
ax=subplot(3,3,1);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL2(dev==1,:),1),'k','LineWidth',2)
plot(bins,mean(Z_MUAoptoL2(dev==1,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(bins,mean(Z_MUAlaserL2(dev==1,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.XLim=[-200 4000];


ax=subplot(3,3,4);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL4(dev==1,:),1),'k','LineWidth',2) 
plot(bins,mean(Z_MUAoptoL4(dev==1,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(bins,mean(Z_MUAlaserL4(dev==1,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.XLim=[-200 4000];

ax=subplot(3,3,7);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL5(dev==1,:),1),'k','LineWidth',2)
plot(bins,mean(Z_MUAoptoL5(dev==1,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(bins,mean(Z_MUAlaserL5(dev==1,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.XLim=[-200 4000];
ax.XLabel.String='Time (ms)';


%%P9-P13
ax=subplot(3,3,2);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL2(dev==2,:),1),'k','LineWidth',2)
plot(bins,mean(Z_MUAoptoL2(dev==2,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(bins,mean(Z_MUAlaserL2(dev==2,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.XLim=[-200 4000];

ax=subplot(3,3,5);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL4(dev==2,:),1),'k','LineWidth',2)
plot(bins,mean(Z_MUAoptoL4(dev==2,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(bins,mean(Z_MUAlaserL4(dev==2,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.XLim=[-200 4000];

ax=subplot(3,3,8);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL5(dev==2,:),1),'k','LineWidth',2)
plot(bins,mean(Z_MUAoptoL5(dev==2,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(bins,mean(Z_MUAlaserL5(dev==2,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.XLim=[-200 4000];
ax.XLabel.String='Time (ms)';

%%P14-P18
ax=subplot(3,3,3);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL2(dev==3,:),1),'k','LineWidth',2)
plot(bins,mean(Z_MUAoptoL2(dev==3,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(bins,mean(Z_MUAlaserL2(dev==3,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.XLim=[-200 4000];
if plotLaserOnly; legend('Control','+ SST opto','Laser only'); else; legend('Control','+ SST opto'); end


ax=subplot(3,3,6);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL4(dev==3,:),1),'k','LineWidth',2)
plot(bins,mean(Z_MUAoptoL4(dev==3,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(bins,mean(Z_MUAlaserL4(dev==3,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.XLim=[-200 4000];

ax=subplot(3,3,9);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL5(dev==3,:),1),'k','LineWidth',2)
plot(bins,mean(Z_MUAoptoL5(dev==3,:),1),'r','LineWidth',2)
if plotLaserOnly; plot(bins,mean(Z_MUAlaserL5(dev==3,:),1),'c','LineWidth',2); end
ax=applyFont(ax,1);
ax.XLim=[-200 4000];
ax.XLabel.String='Time (ms)';

if plotLaserOnly; figname='SST_MUA_layers_dev_withLaser'; else; figname='SST_MUA_layers_dev';end
% sg=sgtitle(strcat('P',int2str(ops.age)));
% sg.FontSize=30;
export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')


%% Compre visual vs opto with subtraction of laser only condition


MUAsubL2=Z_MUAoptoL2-Z_MUAlaserL2;
MUAsubL4=Z_MUAoptoL4-Z_MUAlaserL4;
MUAsubL5=Z_MUAoptoL5-Z_MUAlaserL5;

% MUAsubL2=zscoreBaseline(MUAsubL2);
% MUAsubL4=zscoreBaseline(MUAsubL4);
% MUAsubL5=zscoreBaseline(MUAsubL5);


figure('units','normalized','outerposition',[0 0 1 1]);

%%P5-P8
ax=subplot(3,3,1);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL2(dev==1,:),1),'k','LineWidth',2)
plot(bins,mean(MUAsubL2(dev==1,:),1),'r','LineWidth',2)
ax=applyFont(ax,1);
ax.XLim=[-200 1000];


ax=subplot(3,3,4);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL4(dev==1,:),1),'k','LineWidth',2) 
plot(bins,mean(MUAsubL4(dev==1,:),1),'r','LineWidth',2)
ax=applyFont(ax,1);
ax.XLim=[-200 1000];

ax=subplot(3,3,7);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL5(dev==1,:),1),'k','LineWidth',2)
plot(bins,mean(MUAsubL5(dev==1,:),1),'r','LineWidth',2)
ax=applyFont(ax,1);
ax.XLim=[-200 1000];
ax.XLabel.String='Time (ms)';


%%P9-P13
ax=subplot(3,3,2);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL2(dev==2,:),1),'k','LineWidth',2)
plot(bins,mean(MUAsubL2(dev==2,:),1),'r','LineWidth',2)
ax=applyFont(ax,1);
ax.XLim=[-200 1000];

ax=subplot(3,3,5);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL4(dev==2,:),1),'k','LineWidth',2)
plot(bins,mean(MUAsubL4(dev==2,:),1),'r','LineWidth',2)
ax=applyFont(ax,1);
ax.XLim=[-200 1000];

ax=subplot(3,3,8);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL5(dev==2,:),1),'k','LineWidth',2)
plot(bins,mean(MUAsubL5(dev==2,:),1),'r','LineWidth',2)
ax=applyFont(ax,1);
ax.XLim=[-200 1000];
ax.XLabel.String='Time (ms)';

%%P14-P18
ax=subplot(3,3,3);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL2(dev==3,:),1),'k','LineWidth',2)
plot(bins,mean(MUAsubL2(dev==3,:),1),'r','LineWidth',2)
ax=applyFont(ax,1);
ax.XLim=[-200 1000];
legend('Control','+ SST opto - laser only')


ax=subplot(3,3,6);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL4(dev==3,:),1),'k','LineWidth',2)
plot(bins,mean(MUAsubL4(dev==3,:),1),'r','LineWidth',2)
ax=applyFont(ax,1);
ax.XLim=[-200 1000];

ax=subplot(3,3,9);
patch([0 100 100 0],[-10 -10 1000 1000],'y','EdgeColor','none','HandleVisibility','off') 
hold on;
plot([-50,-50,150,150],[-5,75,75,-5],'b--','LineWidth',2,'HandleVisibility','off')
plot(bins,mean(Z_MUAvisualL5(dev==3,:),1),'k','LineWidth',2)
plot(bins,mean(MUAsubL5(dev==3,:),1),'r','LineWidth',2)
ax=applyFont(ax,1);
ax.XLim=[-200 1000];
ax.XLabel.String='Time (ms)';

if plotLaserOnly; figname='SST_MUA_layers_dev_withSubtraction_Zoom'; else; figname='SST_MUA_layers_dev';end
% sg=sgtitle(strcat('P',int2str(ops.age)));
% sg.FontSize=30;
export_fig(fullfile('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SummaryFigures',figname),'-tiff','-transparent')


function ax=applyFont(ax,normalization)
    ax.XLim=[-50 4000];

    if normalization
        ax.YLabel.String='Z-score';
        ax.YLim=[-5 50];
    else
        ax.YLabel.String='Spike/s'; 
        ax.YLim=[-20 700];
    end
    ax.XAxis.Visible='on';
    ax.FontSize=20;
    ax.Box='off';
    ax.LineWidth = 1.5;


end