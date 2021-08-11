clear 
close all
clc

mainSortingDir='C:\Users\Butt Lab\Documents\SpikeSorting';
data=readtable('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\V1_InVivo_SST;KORD.csv');
data=data((data.Sorting~=0) & (~isnan(data.Sorting)),:);
data.Experiment=fullfile(data.Folder,data.MouseID,data.Experiment);
IDs=unique(data.MouseID);

for recording=1:length(IDs)    
    fprintf('************************************************************************************************************************ \n')
    fprintf(strcat('Looking into ...',IDs{recording},'\n'))
    mouseData=data(strcmp(data.MouseID,IDs{recording}),:);
    
    clear ops stim
    dirIN=mouseData.Experiment;    
    dirOUT=fullfile(mainSortingDir,IDs{recording});
    if ~exist(dirOUT,'dir'); mkdir(dirOUT); end
        
    if ~exist(fullfile(dirOUT,'RawData.dat'),'file')
        %% Load electrode map
        if strcmp(mouseData.Electrode{1},'A1x32')
            load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\A1x32_Map.mat');
        elseif strcmp(mouseData.Electrode{1},'OA1x32')
            load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\OA1x32_Map.mat');
        elseif strcmp(mouseData.Electrode{1},'A1x32_Poly3')
            load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\A1x32_Poly3_Map.mat');
        end
        ops.ElectrodeMap=ElectrodeMap; 

        %% Build binary file
        ops.Protocol=mouseData.Protocol;
        ops.age = mouseData.Age(1);
        fprintf(strcat('Build binary *.dat file for ...',IDs{recording},'\n'))
        tic; 
        
        ops = makeBinary(ops,dirIN,dirOUT);
        save(fullfile(dirOUT, 'ops.mat'), 'ops', '-v7.3');
        fprintf('Time %3.0fs. Finished building binary file... \n', toc);
    end
    
    %% Build stimuli structure
    if ~exist(fullfile(dirOUT,'stimuli.mat'),'file')
        fprintf(strcat('Finding stimuli for ...',IDs{recording},'\n'))
        load(fullfile(dirOUT, 'ops.mat'), 'ops');
        stim = getStimuli(ops,mouseData); 
        save(fullfile(dirOUT, 'stimuli.mat'), 'stim', '-v7.3')
        fprintf('Time %3.0fs. Saved stimuli... \n', toc);
    end

    
end