clear 
close all
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\Kilosort2')) % path to kilosort folder

mainSortingDir='C:\Users\Butt Lab\Documents\SpikeSorting';
load('InVivo_V1_SST;Ai32.mat'); 
tab=tab(tab.Use==1,:);
recToAnalyse=unique(tab.RecID);

for recording=1:length(recToAnalyse)
    subtab=tab(tab.RecID==recToAnalyse(recording),:);
    
    dataFolders=fullfile(subtab.Folder{1},subtab.Experiment(1:height(subtab)));
%     sortingFolder=strcat(subtab.Folder{1},'\Sorting');
    sortingFolder=fullfile(mainSortingDir,subtab.MouseID{1});
    if ~exist(sortingFolder,'dir')
        mkdir(sortingFolder)
    end
    
    master_kilosort_Fun(dataFolders,sortingFolder)   
    
end
