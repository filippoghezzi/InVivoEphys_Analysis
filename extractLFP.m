clear
close all
clc

origFolder=cd;
load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\A1x32_Map.mat');

dataFolder='C:\Users\Butt Lab\Desktop\In vivo data';
cd(dataFolder);
dataFolderContent=dir;

for i=1:length(dataFolderContent)
    if contains(dataFolderContent(i).name,'2019')
        recFolder=strcat(dataFolder,'\',dataFolderContent(i).name);
        cd(recFolder)
        recContent=dir;
        dataPresence=[];
        for j=1:length(recContent)
            dataPresence=[dataPresence;strcmp(recContent(j).name,'Data.mat')];
        end
        if ~any(dataPresence)
            [data,dataTime,eventArray]=loadData(recFolder,ElectrodeMap);
            save('Data.mat','data','dataTime','eventArray')
        end
    cd(dataFolder);
    end
end

cd(origFolder)


