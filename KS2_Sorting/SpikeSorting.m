clear 
close all

mainSortingDir='C:\Users\Butt Lab\Documents\SpikeSorting';

data=readtable('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\V1_InVivo_SST;Ai32.csv');
data=data((data.Sorting~=0) & (~isnan(data.Sorting)),:);
ID=unique(data.MouseID);

for recording=1:length(ID)    
    fprintf('************************************************************************************************************************ \n')
    fprintf(strcat('Sorting with KS2 ...',ID{recording},'\n'))

    dirOUT=fullfile(mainSortingDir,ID{recording});
        
    if ~exist(fullfile(dirOUT,'RawData.dat'),'file')
        error('RawData.dat file does not exist! - Run preprocessing first')
    end
    
        
    %% Do Sorting
    if ~exist(fullfile(dirOUT,'rez2.mat'),'file')
        load(fullfile(dirOUT,'ops.mat'))
        master_kilosort_Fun(ID{recording},dirOUT,ops)   
    end
    
end
fprintf('************************************************************************************************************************ \n')