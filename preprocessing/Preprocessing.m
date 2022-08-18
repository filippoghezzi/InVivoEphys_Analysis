clear 
close all
clc

mainSortingDir='C:\Users\Butt Lab\Documents\SpikeSorting';
data=readtable('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\V1_InVivo.csv');
data=data(data.Preprocessing==1,:);
data.Experiment=fullfile(data.Folder,data.MouseID,data.Experiment);
IDs=unique(data.MouseID);

for recording=1:length(IDs)    
    fprintf('************************************************************************************************************************ \n')
    fprintf(strcat('Looking into ...',IDs{recording},'\n'))
    mouseData=data(strcmp(data.MouseID,IDs{recording}),:);
    tic; 

    clear ops stim
    dirIN=mouseData.Experiment;    
    dirOUT=fullfile(mainSortingDir,IDs{recording});
    if ~exist(dirOUT,'dir'); mkdir(dirOUT); end
        
%     if ~exist(fullfile(dirOUT,'RawData.dat'),'file')
%         %% Load electrode map
%         ops.reorderChannels = 1;
%         if strcmp(mouseData.Electrode{1},'A1x32')
%             load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\A1x32_Map.mat');
%             ops.ElectrodeStart = 16;
%         elseif strcmp(mouseData.Electrode{1},'OA1x32')
%             load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\OA1x32_Map.mat');
%             ops.ElectrodeStart = 16;
%         elseif strcmp(mouseData.Electrode{1},'A1x32_Poly3')
%             load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\A1x32_Poly3_Map.mat');
%             ops.ElectrodeStart = 16;
%         elseif strcmp(mouseData.Electrode{1},'A1x32_Intan32')
%             load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\A1x32_Intan32_Map.mat');
%             ops.ElectrodeStart = 64;
%         elseif strcmp(mouseData.Electrode{1},'A1x32_NH')
%             load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\A1x32_Map.mat');
%             ops.ElectrodeStart = 0;
%             ElectrodeMap=ElectrodeMap-16;
%             ops.reorderChannels = 0;
%         end
%         ops.ElectrodeMap=ElectrodeMap; 
% 
%         %% Build binary file
%         ops.Protocol=mouseData.Protocol;
%         ops.age = mouseData.Age(1);
%         fprintf(strcat('Build binary *.dat file for ...',IDs{recording},'\n'))        
%         ops = makeBinary(ops,dirIN,dirOUT);
%         save(fullfile(dirOUT, 'ops.mat'), 'ops', '-v7.3');
%         fprintf('Time %3.0fs. Finished building binary file... \n', toc);
%     end
    
    %% Build stimuli structure
    if ~exist(fullfile(dirOUT,'stimuli.mat'),'file')
        fprintf(strcat('Finding stimuli for ...',IDs{recording},'\n'))
        load(fullfile(dirOUT, 'ops.mat'), 'ops');
        stim = getStimuli(ops,mouseData); 
        save(fullfile(dirOUT, 'stimuli.mat'), 'stim', '-v7.3')
        fprintf('Time %3.0fs. Saved stimuli... \n', toc);
    end   
end