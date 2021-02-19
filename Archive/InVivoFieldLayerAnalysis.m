clear 
clc
close all

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 


% rec='SC35';
% 
% dataFolder='C:\Users\Butt Lab\Documents\SpikeSorting';
% dataFolder='D:\InVivo';
% dir=fullfile(dataFolder,rec);



tab=readtable('V1_InVivo_SST;Ai32_reduced.csv');
% rec={'SC3','SC4','SC9','SC18','SC19'};

% spikeFolder='C:\Users\Butt Lab\Documents\SpikeSorting';

for i=1:height(tab)
    dir=fullfile(tab.Folder{i},tab.ID{i});
%%
    load(fullfile(dir,'stimuli.mat'),'stim')
    load(fullfile(dir,'rez.mat'),'rez')
    s=loadSpikesNew(dir);
    
    ops=rez.ops;
    ops.binaryRoot=dir;
    ops.fbinary=fullfile(ops.binaryRoot,'RawData.dat');
    ops.dirOUT=fullfile(dir,'Figures');
    if ~exist(ops.dirOUT,'dir')
        mkdir(ops.dirOUT) 
    end
    ops.chanMap=(1:ops.NchanTOT)';
    ops.LFPwindow=[1 5]; %s
    ops.PSTHbinSize=0.003;
    ops.electrodeLength = 800; %um
    ops.electrodeSpacing = 25*10^-6; %um
    ops.verbose=1;
    
    ops.L2=1:min(ops.L4)-1;
    ops.L5=max(ops.L4)+1:ops.NchanTOT;
    
    %Select index for plotting a single trace per layer
    L4Idx=ops.L4(2);
    L2Idx=min(ops.L4)-4;
    while L2Idx<1
        L2Idx=L2Idx+1;
        if L2Idx==L4Idx
            error('L4 index is too low')
        end
    end
    L5Idx=max(ops.L4)+6;
    while L5Idx>ops.NchanTOT
        L5Idx=L5Idx-1;
        if L5Idx==L4Idx
            error('L4 index is too high')
        end
    end
    
    
    eLFP = get_eLFP(ops.fbinary,stim.ledR,ops.fs,ops.NchanTOT,ops.LFPwindow,'LFP');
    eMUA_trace=get_eLFP(ops.fbinary,stim.ledR(10),ops.fs,ops.NchanTOT,ops.LFPwindow,'MUA');

    MUA(:,1)=getLayerMUA(s,stim.ledR,ops.fs,ops.LFPwindow,ops.PSTHbinSize,ops.L2);
    MUA(:,2)=getLayerMUA(s,stim.ledR,ops.fs,ops.LFPwindow,ops.PSTHbinSize,ops.L4);
    MUA(:,3)=getLayerMUA(s,stim.ledR,ops.fs,ops.LFPwindow,ops.PSTHbinSize,ops.L5);
    
    plotLFPbyLayer(eLFP([L2Idx,L4Idx,L5Idx],:,:),eMUA_trace([L2Idx,L4Idx,L5Idx],:),MUA,ops); 
    
    %% Time frequency analysis
    getSpectrogram(squeeze(eLFP(L4Idx,:,:)),ops.fs,ops.LFPwindow)

end

