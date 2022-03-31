clear 
% close all
addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\spikes'))
addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\npy-matlab'))
mainSortingDir='C:\Users\Butt Lab\Documents\SpikeSorting';
data=readtable('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\V1_InVivo.csv');
data=data((data.Sorting~=0) & (~isnan(data.Sorting)),:);
ID=unique(data.MouseID);

for recording=1:length(ID)    
    fprintf('************************************************************************************************************************ \n')
    fprintf(strcat('Sorting with KS2 ...',ID{recording},'\n'))
    dirOUT=fullfile(mainSortingDir,ID{recording});   
    if ~exist(fullfile(dirOUT,'RawData.dat'),'file'); error('RawData.dat file does not exist! - Run preprocessing first'); end
        
    %% Do Sorting
    if ~exist(fullfile(dirOUT,'rez2.mat'),'file')
        load(fullfile(dirOUT,'ops.mat'))
        master_kilosort_Fun(ID{recording},dirOUT,ops)  
        
        %% Plot drift map
        [spikeTimes, spikeAmps, spikeDepths, ~] = ksDriftmap(dirOUT);
        figure; plotDriftmap(spikeTimes, spikeAmps, spikeDepths);
        title(ID{recording},'FontSize',20)
    end    
end
fprintf('************************************************************************************************************************ \n')