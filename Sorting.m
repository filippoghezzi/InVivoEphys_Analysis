clear 
close all
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\Kilosort2')) % path to kilosort folder

mainSortingDir='C:\Users\Butt Lab\Documents\SpikeSorting';
load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\InVivo_V1_SST;Ai32.mat'); 
tab=tab(tab.Use~=0,:);
recToAnalyse=unique(tab.RecID);

for recording=1:length(recToAnalyse)
    subtab=tab(tab.RecID==recToAnalyse(recording),:);
    dataFolders=fullfile(subtab.Folder{1},subtab.Experiment(1:height(subtab)));
    sortingFolder=fullfile(mainSortingDir,subtab.MouseID{1});
    if ~exist(sortingFolder,'dir')
        mkdir(sortingFolder)
        master_kilosort_Fun(dataFolders,sortingFolder)   
    end
end
