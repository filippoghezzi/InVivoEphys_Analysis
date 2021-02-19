clear 
clc
close all

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 

recordings=readtable('V1_InVivo_SST;Ai32.csv');
recordings=recordings(~isnan(recordings.Sorting),:);
recID=unique(recordings.MouseID);
% recID={'SC28'};%,'NK7'};

recFolder='D:\InVivo_V1';

for i=1:numel(recID)
    %% Load data
    fprintf('analysing %s',recID{i})
    dir=fullfile(recFolder,recID{i});
    verbose=1;
    load(fullfile(dir,'rez.mat'),'rez')
    load(fullfile(dir,'stimuli.mat'),'stim')
    ops=rez.ops;
    ops.dirOUT=fullfile(dir,'Figures');


    s=loadSpikes(dir);


    %% 
    if ~isempty(s.suid)
        s.suage=ones(numel(s.suid),1)*ops.age;
        s.sulayer=findSingleUnitLayer(s,ops);        
        s=getTemplateFeatures(s,ops,0,ops.dirOUT);

        s = PSTHbyCondition (s,stim,ops,0 );

    end

    save(fullfile(dir,'spikes.mat'),'s', '-v7.3')
end