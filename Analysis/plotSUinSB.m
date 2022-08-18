close all
clear
clc
addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 

%% Load
dir=fullfile('G:\InVivo_SpikeSorting_V1\SC38');
[ops,s,stim] = InVivo_dataProcessing_loadData(dir);
LFP = loadLFP_baseline(ops.fbinary,ops.fs,ops.NchanTOT,0,20*60*ops.fs,'LFP');  

%% Plot
figure
hold on
plot(LFP(ops.L4best,ops.fs*60*12:ops.fs*60*17),'k')
plot(LFP(ops.L4best-6,ops.fs*60*12:ops.fs*60*17)+4000,'k')
plot(LFP(ops.L4best+8,ops.fs*60*12:ops.fs*60*17)-4000,'k')
plot(LFP(ops.L4best+16,ops.fs*60*12:ops.fs*60*17)-8000,'k')

U78=s.st(s.sclu==78);
U78=U78((U78>ops.fs*60*12) & (U78<ops.fs*60*17))-ops.fs*60*12;
plot([U78,U78],[1500,2000],'r-')

U97=s.st(s.sclu==97);
U97=U97((U97>ops.fs*60*12) & (U97<ops.fs*60*17))-ops.fs*60*12;
plot([U97,U97],[-2700,-2200],'r-')

U59=s.st(s.sclu==59);
U59=U59((U59>ops.fs*60*12) & (U59<ops.fs*60*17))-ops.fs*60*12;
plot([U59,U59],[-6300,-5800],'r-')

U83=s.st(s.sclu==83);
U83=U83((U83>ops.fs*60*12) & (U83<ops.fs*60*17))-ops.fs*60*12;
plot([U83,U83],[-10000,-9500],'r-')

plot([ops.fs*60*2,ops.fs*60*2.5],[-11,-11]*10^3,'k','LineWidth',1)% 30 s
plot([ops.fs*60*2,ops.fs*60*2],[-11,-10]*10^3,'k','LineWidth',1)% 1 mV

%% Save
ax=gca;
ax.XLim=[3,9]*10^6;
ax.Box='off';
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';
print(gcf,'-dpdf','C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\WT Doctoral Programme in Neuroscience\DPhil in Neuroscience\Manuscripts\V1 S1\Figures\Fig. 2\SB and units')

ax.XLim=[4.45,4.65]*10^6;
plot([4.5,4.55]*10^6,[-10.5,-10.5]*10^3,'k','LineWidth',1) % 1.5 s
print(gcf,'-dpdf','C:\Users\Butt Lab\OneDrive - OnTheHub - The University of Oxford\University of Oxford\WT Doctoral Programme in Neuroscience\DPhil in Neuroscience\Manuscripts\V1 S1\Figures\Fig. 2\SB and units_Zoom')
