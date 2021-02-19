close all
clear 
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 

recordings=readtable('V1_InVivo_SST;Ai32.csv');
recID=unique(recordings.MouseID);
recFolder='D:\InVivo_SST';

mouseID=[];
age=[];

corticalResponse=[];
VEP=[];
VEPtime=[];
PSD=[];



for i=1:numel(recID)
        dir=fullfile(recFolder,recID{i});

        load(fullfile(dir,'rez.mat'),'rez')
        ops=rez.ops;
        
        mouseID=[mouseID;recID(i)];
        age=[age;ops.age];
        corticalResponse=[corticalResponse,ops.corticalResponseOnset];
        VEP=[VEP,ops.peakVEP];
        VEPtime=[VEPtime,ops.timeVEP];
        
        PSD=[PSD;mean(ops.LFP_PSD(ops.L4,:))];

end
%% Set group logic arrays

%Age
P5P8=age<9;
P9P13=age>=9 & age<14;
P9P11=age>=9 & age<12;
P12P13=age>=12 & age<14;
P14P18=age>=14;

dev(P5P8,1)=categorical(cellstr('P5-P8'));
dev(P9P13,1)=categorical(cellstr('P9-P13'));
dev(P14P18,1)=categorical(cellstr('P14-P18'));
dev=reordercats(dev,{'P5-P8','P9-P13','P14-P18'});

dev2(P5P8,1)=categorical(cellstr('P5-P8'));
dev2(P9P11,1)=categorical(cellstr('P9-P11'));
dev2(P12P13,1)=categorical(cellstr('P12-P13'));
dev2(P14P18,1)=categorical(cellstr('P14-P18'));
dev2=reordercats(dev2,{'P5-P8','P9-P11','P12-P13','P14-P18'});

figure
boxchartWithScatter(VEPtime,dev2);
figure
boxchartWithScatter(VEP,dev2);
figure
boxchartWithScatter(corticalResponse,dev2);