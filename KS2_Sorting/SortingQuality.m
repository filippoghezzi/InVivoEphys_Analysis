clear 
close all
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\sortingQuality'))
addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis'))


mainSortingDir='C:\Users\Butt Lab\Documents\SpikeSorting';
recName={'SC1','SC2','SC3','SC4','SC5','SC6','SC7','SC8','SC9','SC10','SC11','SC12','SC13','SC14','SC15','SC18','SC19','SC20','SC21','SC22','SC23','SC24','SC25','SC26','SC27','SC28','SC29','SC30','SC31','SC32','SC33','SC34','SC35','SC36','SC37'};

for rec=1:size(recName,2)
    results=table;
    fprintf('************************************************************************************************************************ \n')
    fprintf(strcat('SortingQuality ...',recName{rec},'\n'))
    recFolder=fullfile(mainSortingDir,recName{rec});
    [cids, cgs, uQ, cR, isiV] = sqKilosort.computeAllMeasures(recFolder);
    
    %Make output
    results.ClusterID=cids';
    results.ClusterGroup=cgs';
    results.ClusterQuality=uQ;
    results.ClusterContamination=cR;
    results.ClusterViolationRate=isiV';
    writetable(results,fullfile(recFolder,'SortingQuality.csv'))
    fprintf('************************************************************************************************************************ \n')
end
