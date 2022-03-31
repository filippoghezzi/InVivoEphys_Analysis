function [ops,s,stim]=InVivo_dataProcessing_loadData(folder)


    %% Stimuli
    load(fullfile(folder,'stimuli.mat'),'stim')
    
    %% Ops
    load(fullfile(folder,'rez.mat'),'rez') 
    ops=rez.ops;
    if ~isfield(ops,'Protocol')
        ops.Protocol={'Control'};%%%%%%%%%%%%%%%%%%%
    end
    ops.verbose=0;  
    ops.chanMap=(1:ops.NchanTOT)';
    ops.recID=split(folder,'\');
    ops.recID=ops.recID{end};
    ops.binaryRoot=folder;
    ops.fbinary=fullfile(ops.binaryRoot,'RawData.dat');
    ops.dirOUT=fullfile(folder,'Figures');
    if ~exist(ops.dirOUT,'dir'); mkdir(ops.dirOUT); end
    
    %% Stimulus-evoked LFP 
    ops.LFPwindow=[1 5]; %s
    ops.PSTHbinSize=0.003;
    
    %% Set electrode map and features
    if ops.ElectrodeMap(1)==32
        ops.electrodeType = 'Poly3';
        ops.electrodeLength = 300; %um
        ops.electrodeSpacing = 25*10^-6; %um
        ops.idxCentralLine=[1,2,5,8,11,14,17,20,23,26,29,32];
        ops.electrodeChannelsForCSD = logical(zeros(ops.NchanTOT,1));
        ops.electrodeChannelsForCSD(ops.idxCentralLine) = true;
    else 
        ops.electrodeType = 'Poly2';
        ops.electrodeLength = 800; %um
        ops.electrodeSpacing = 25*10^-6; %um
        ops.electrodeChannelsForCSD = logical(ones(ops.NchanTOT,1));
    end
    
    %% Chemogenetics
    if any(cellfun(@(x) x(end), ops.Protocol)=='K') 
        ops.SalB = 1;
    else 
        ops.SalB = 0;
    end
    
    %% Spikes
    s=loadSpikes(folder);
    
    if strcmp(ops.recID,'SC107l')
        s.suDepth=[0;-50;-100;-100;-150;-150;-200;-200;-200;-200;-200;-200;-250;-300;-350;-300;-300;-400;-400;-400;-100;-200;-400;-350;-350;-50];
    elseif strcmp(ops.recID,'SC107u')
        s.suDepth=[0;-50;-50;-50;-50;-100;-100;-100;-150;-150;-150;-200;-150;-200;-300;-300;-350;-300;-50;-200;-50;-200];
    end
    
    
end