clear 
close all
clc

mainSortingDir='C:\Users\Butt Lab\Documents\SpikeSorting';

data=readtable('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\V1_InVivo_SST;Ai32.csv'); 
data=data(data.Use~=0,:);
data.Experiment=fullfile(data.Folder,data.Experiment);

recToAnalyse=unique(data.MouseID);
for recording=1:length(recToAnalyse)    
    fprintf('************************************************************************************************************************ \n')
    fprintf(strcat('Sorting ...',recToAnalyse{recording},'\n'))
    dataFolders=data(strcmp(data.MouseID,recToAnalyse{recording}),:).Experiment;    
    sortingFolder=fullfile(mainSortingDir,recToAnalyse{recording});
    if ~exist(sortingFolder,'dir')
        mkdir(sortingFolder)
        master_kilosort_Fun(recToAnalyse{recording},dataFolders,sortingFolder)   
    end
end
fprintf('************************************************************************************************************************ \n')