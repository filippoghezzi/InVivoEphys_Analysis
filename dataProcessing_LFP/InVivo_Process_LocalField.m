clear 
clc
close all

addpath(genpath('  C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 

recordings=readtable('V1_InVivo_SST;Ai32.csv');
recordings=recordings(~isnan(recordings.Sorting),:);
recID=unique(recordings.MouseID);
recID={'SC33'};%,'NK7'};
recFolder='D:\InVivo_V1';

for i=1:numel(recID)
    %% Load data
    dir=fullfile(recFolder,recID{i});
    load(fullfile(dir,'stimuli.mat'),'stim')
    load(fullfile(dir,'rez.mat'),'rez')
    s=loadSpikes(dir);
    
    ops=rez.ops;
    ops.binaryRoot=dir;
    ops.fbinary=fullfile(ops.binaryRoot,'RawData.dat');
    ops.dirOUT=fullfile(dir,'Figures');
    ops.recID=recID{i};
    if ~exist(ops.dirOUT,'dir'); mkdir(ops.dirOUT); end
    
    %% Analyse Stimulus Evoked Field Potential
    ops=analyseStimulusEvokedFieldPotential(ops,s,stim);
    
    %% Analyse baseline
    ops=analyseBaseline(ops,s);
    
    %% Save outputs
    rez.ops=ops;
    save(fullfile(dir,'rez.mat'),'rez')

end

