clear 
clc
close all

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 
addpath('C:\Users\Butt Lab\Documents\MATLAB\Add-Ons\CircStat2012a')

recordings=readtable('V1_InVivo.csv');
recordings=recordings(~isnan(recordings.Sorting),:);

recID=unique(recordings.MouseID);
% recID={'NK32'};

for i=1:numel(recID)
    fprintf(strcat('Looking into ...',recID{i},'\n'))
    mouseData=recordings(strcmp(recordings.MouseID,recID{i}),:);
    recFolder=mouseData.Folder{1};
    dir=fullfile(recFolder,recID{i});
    outputFileName=fullfile(dir,'processData.mat');

    %% Load data
    if exist(outputFileName)==2
        load(outputFileName,'results')
    else
        results=struct;
        results.baseline=struct;
        results.SalB=struct;
        results.SalB.baseline=struct;
    end
    [ops,s,stim] = InVivo_dataProcessing_loadData(dir);
    
    %% Select stimuli
    if strcmp(mouseData.BrainArea{1},'S1BF')
        stimControl=stim.whiskerStim;
        if ops.SalB; stimChemo=stim.whiskerStim_SalB; end
        ops.brainArea='S1BF';
    elseif strcmp(mouseData.BrainArea{1},'V1')
        stimControl=stim.ledR;
        if ops.SalB; stimChemo=stim.ledR_SalB; end
        ops.brainArea='V1';
    end

%     %% Analyse Stimulus Evoked Field Potential  
%     [results.evokedLFP, ops]=InVivo_dataProcessing_evokedLFP(ops,s,stimControl,'Control');
%     rez.ops=ops;
%     save(fullfile(dir,'rez.mat'),'rez')

    %% Analyse baseline
%     results.baseline=InVivo_dataProcessing_baseline(ops,s,results.baseline,'Condition','Control');
 
    %% Analyse single units
    results.s=InVivo_dataProcessing_singleUnits(s,ops,stim);
    
    %% Chemogenetic
    
%     if ops.SalB
%         fprintf(strcat('Analysing chemogenetic condition...','\n'))
%         ops_SalB=ops;   
%         ops_SalB=rmfield(ops_SalB,'L4');

%         [results.SalB.evokedLFP, ops_SalB]=InVivo_dataProcessing_evokedLFP(ops_SalB,s,stimChemo,'SalB');
%         results.SalB.baseline=InVivo_dataProcessing_baseline(ops_SalB,s,results.SalB.baseline,'Condition','SalB');
%         results.SalB.ops=ops_SalB;
%     end
    
    %% Save outputs
    save(fullfile(dir,'processData.mat'),'results')

end

