clear 
clc
close all

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 
addpath('C:\Users\Butt Lab\Documents\MATLAB\Add-Ons\CircStat2012a')

recordings=readtable('V1_InVivo_SST;Ai32.csv');
recordings=recordings(~isnan(recordings.Sorting),:);
recID=unique(recordings.MouseID);
recFolder='Q:\';
% recID={'NK58','NSC1','NSC2','NSC4','NSC5','SC55','SC56','SC91','SC92b','SC93b','SC107l','SC107u'};
% recFolder='C:\Users\Butt Lab\Documents\SpikeSorting';

for i=1:numel(recID)
    fprintf(strcat('Looking into ...',recID{i},'\n'))
    dir=fullfile(recFolder,recID{i});
    results=struct;

    %% Load data
    [ops,s,stim] = InVivo_dataProcessing_loadData(dir);
    
    %% Analyse Stimulus Evoked Field Potential
    [results.evokedLFP, ops]=InVivo_dataProcessing_evokedLFP(ops,s,stim.ledR,'Control');
    rez.ops=ops;
    save(fullfile(dir,'rez.mat'),'rez')

    %% Analyse baseline
    results.baseline=InVivo_dataProcessing_baseline(ops,s,'Condition','Control');
 
    %% Analyse single units
    results.s=InVivo_dataProcessing_singleUnits(s,ops,stim);
    
    %% Chemogenetic
    if ops.SalB
        fprintf(strcat('Analysing chemogenetic condition...','\n'))
%         ops=rmfield(ops,'L4');
        [results.SalB.evokedLFP, ops_SalB]=InVivo_dataProcessing_evokedLFP(ops,s,stim.ledR_SalB,'SalB');
        results.SalB.baseline=InVivo_dataProcessing_baseline(ops_SalB,s,'Condition','SalB');
        results.SalB.ops=ops_SalB;
    end
    
    %% Save outputs
    save(fullfile(dir,'processData.mat'),'results')

end

