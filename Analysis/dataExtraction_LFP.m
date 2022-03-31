close all
clear 
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 

recordings=readtable('V1_InVivo.csv');
recordings=recordings(~isnan(recordings.Sorting),:);
recID=unique(recordings.MouseID);
% recID={'K6','K24','K35','K36','K38','K39'};
% recFolder='D:\InVivo_V1';

mouseID=[];
optotagging=[];
mouseAge=[];
brainArea=[];
state=[];

MUA_L4=[];
MUA_L4_peakFast=[];
MUA_L4_peakSlow=[];
MUA_L4_baselineFiring=[];
MUA_L4_continuity=[];
LFP_L4_onset=[];
LFP_L4_peak=[];
LFP_L4_latency=[];
SB_L4_power=[];
SB_L4_duration=[];
SB_L4_frequency=[];
PSD=[];

MUA_L4_K=[];
MUA_L4_peakFast_K=[];
MUA_L4_peakSlow_K=[];
MUA_L4_baselineFiring_K=[];
MUA_L4_continuity_K=[];
LFP_L4_onset_K=[];
LFP_L4_peak_K=[];
LFP_L4_latency_K=[];
SB_L4_power_K=[];
SB_L4_duration_K=[];
SB_L4_frequency_K=[];
PSD_K=[];

for i=1:numel(recID)
    fprintf(strcat('Looking into ...',recID{i},'\n'))

    thisTagging=recordings.Tagging(strcmp(recordings.MouseID,recID{i}),:);
    thisTagging=thisTagging(1);
    recFolder=recordings.Folder(strcmp(recordings.MouseID,recID{i}),:);
    recFolder=recFolder{1};
    Age=recordings.Age(strcmp(recordings.MouseID,recID{i}),:);
    Age=Age(1);
    Area=recordings.BrainArea(strcmp(recordings.MouseID,recID{i}),:);
    Area=categorical(cellstr((Area{1})));
    State=recordings.State(strcmp(recordings.MouseID,recID{i}),:);
    State=categorical(cellstr((State{1})));
    
    dir=fullfile(recFolder,recID{i});
    load(fullfile(dir,'processData.mat'),'results')
    load(fullfile(dir,'rez.mat'),'rez')
    
    mouseID=[mouseID;{recID{i}}];
    optotagging=[optotagging;thisTagging];
    mouseAge=[mouseAge;Age];
    brainArea=[brainArea;Area];
    state=[state;State];
    
    selectedMUA=results.evokedLFP.MUA.raw(rez.ops.L4best,:);
    MUA_L4=[MUA_L4;selectedMUA];
    MUA_L4_peakFast=[MUA_L4_peakFast;max(selectedMUA(results.evokedLFP.MUA.bins>0 & results.evokedLFP.MUA.bins<200))]; 
    MUA_L4_peakSlow=[MUA_L4_peakSlow;max(selectedMUA(results.evokedLFP.MUA.bins>200 & results.evokedLFP.MUA.bins<4000))]; 

    MUA_L4_baselineFiring=[MUA_L4_baselineFiring;results.baseline.spectral.MUA_spikeRate(rez.ops.L4best)];
    MUA_L4_continuity=[MUA_L4_continuity;results.baseline.spectral.MUA_continuity(rez.ops.L4best)];

    LFP_L4_onset=[LFP_L4_onset;results.evokedLFP.LFP.corticalResponseOnset];
    LFP_L4_peak=[LFP_L4_peak;results.evokedLFP.LFP.peakVEP];
    LFP_L4_latency=[LFP_L4_latency;results.evokedLFP.LFP.timeVEP];

    SB_L4_power=[SB_L4_power;mean(results.baseline.spindleBurst.amplitude)];
    SB_L4_duration=[SB_L4_duration;mean(results.baseline.spindleBurst.duration)];
    SB_L4_frequency=[SB_L4_frequency;numel(results.baseline.spindleBurst.amplitude)/results.baseline.durationBaseline];
    
    PSD = [PSD; results.baseline.spectral.LFP_PSD]; 
    
    if isfield(results,'SalB')
    if isfield(results.SalB.baseline,'durationBaseline')

        selectedMUA_K=results.SalB.evokedLFP.MUA.raw(rez.ops.L4best,:);
        MUA_L4_K=[MUA_L4_K;selectedMUA_K];
        MUA_L4_peakFast_K=[MUA_L4_peakFast_K;max(selectedMUA_K(results.SalB.evokedLFP.MUA.bins>0 & results.SalB.evokedLFP.MUA.bins<200))]; 
        MUA_L4_peakSlow_K=[MUA_L4_peakSlow_K;max(selectedMUA_K(results.SalB.evokedLFP.MUA.bins>200 & results.SalB.evokedLFP.MUA.bins<4000))]; 

        MUA_L4_baselineFiring_K=[MUA_L4_baselineFiring_K;results.SalB.baseline.spectral.MUA_spikeRate(rez.ops.L4best)];
        MUA_L4_continuity_K=[MUA_L4_continuity_K;results.SalB.baseline.spectral.MUA_continuity(rez.ops.L4best)];

        LFP_L4_onset_K=[LFP_L4_onset_K;results.SalB.evokedLFP.LFP.corticalResponseOnset];
        LFP_L4_peak_K=[LFP_L4_peak_K;results.SalB.evokedLFP.LFP.peakVEP];
        LFP_L4_latency_K=[LFP_L4_latency_K;results.SalB.evokedLFP.LFP.timeVEP];

        SB_L4_power_K=[SB_L4_power_K;mean(results.SalB.baseline.spindleBurst.amplitude)];
        SB_L4_duration_K=[SB_L4_duration_K;mean(results.SalB.baseline.spindleBurst.duration)];
        SB_L4_frequency_K=[SB_L4_frequency_K;numel(results.SalB.baseline.spindleBurst.amplitude)/results.baseline.durationBaseline];
        
        PSD_K=[PSD_K; results.SalB.baseline.spectral.LFP_PSD];

    else
        MUA_L4_K=[MUA_L4_K;nan(1,2000)];
        MUA_L4_peakFast_K=[MUA_L4_peakFast_K;nan(1,1)]; 
        MUA_L4_peakSlow_K=[MUA_L4_peakSlow_K;nan(1,1)]; 

        MUA_L4_baselineFiring_K=[MUA_L4_baselineFiring_K;nan(1,1)];
        MUA_L4_continuity_K=[MUA_L4_continuity_K;nan(1,1)];

        LFP_L4_onset_K=[LFP_L4_onset_K;nan(1,1)];
        LFP_L4_peak_K=[LFP_L4_peak_K;nan(1,1)];
        LFP_L4_latency_K=[LFP_L4_latency_K;nan(1,1)];

        SB_L4_power_K=[SB_L4_power_K;nan(1,1)];
        SB_L4_duration_K=[SB_L4_duration_K;nan(1,1)];
        SB_L4_frequency_K=[SB_L4_frequency_K;nan(1,1)];
        
        PSD_K=[PSD_K; nan(size(results.baseline.spectral.LFP_PSD))];

    end
    else
        MUA_L4_K=[MUA_L4_K;nan(1,2000)];
        MUA_L4_peakFast_K=[MUA_L4_peakFast_K;nan(1,1)]; 
        MUA_L4_peakSlow_K=[MUA_L4_peakSlow_K;nan(1,1)]; 

        MUA_L4_baselineFiring_K=[MUA_L4_baselineFiring_K;nan(1,1)];
        MUA_L4_continuity_K=[MUA_L4_continuity_K;nan(1,1)];

        LFP_L4_onset_K=[LFP_L4_onset_K;nan(1,1)];
        LFP_L4_peak_K=[LFP_L4_peak_K;nan(1,1)];
        LFP_L4_latency_K=[LFP_L4_latency_K;nan(1,1)];

        SB_L4_power_K=[SB_L4_power_K;nan(1,1)];
        SB_L4_duration_K=[SB_L4_duration_K;nan(1,1)];
        SB_L4_frequency_K=[SB_L4_frequency_K;nan(1,1)];
        
        PSD_K=[PSD_K; nan(size(results.baseline.spectral.LFP_PSD))];
    end
end

%% Make table
data=table;
data.MouseID=categorical(cellstr(mouseID));
data.Tagging=categorical(cellstr(optotagging));
data.brainArea=categorical(cellstr(brainArea));
data.state=categorical(cellstr(state));

data.mouseAge=mouseAge;

data.MUA_L4_peakFast=MUA_L4_peakFast;
data.MUA_L4_peakSlow=MUA_L4_peakSlow;
data.MUA_L4_baselineFiring=MUA_L4_baselineFiring;
data.MUA_L4_continuity=MUA_L4_continuity;
data.LFP_L4_onset=LFP_L4_onset;
data.LFP_L4_peak=LFP_L4_peak;
data.LFP_L4_latency=LFP_L4_latency;
data.SB_L4_power=SB_L4_power;
data.SB_L4_duration=SB_L4_duration;
data.SB_L4_frequency=SB_L4_frequency;
data.PSD=PSD;

data.MUA_L4_peakFast_K=MUA_L4_peakFast_K;
data.MUA_L4_peakSlow_K=MUA_L4_peakSlow_K;
data.MUA_L4_baselineFiring_K=MUA_L4_baselineFiring_K;
data.MUA_L4_continuity_K=MUA_L4_continuity_K;
data.LFP_L4_onset_K=LFP_L4_onset_K;
data.LFP_L4_peak_K=LFP_L4_peak_K;
data.LFP_L4_latency_K=LFP_L4_latency_K;
data.SB_L4_power_K=SB_L4_power_K;
data.SB_L4_duration_K=SB_L4_duration_K;
data.SB_L4_frequency_K=SB_L4_frequency_K;
data.PSD_K=PSD_K;

data.MUA_L4=MUA_L4;
data.MUA_L4_K=MUA_L4_K;

PSTHbins=results.evokedLFP.MUA.bins;
PSD_f=results.baseline.spectral.LFP_PSD_f;
%% Save
save('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\LFPData.mat', 'data', 'PSTHbins','PSD_f')