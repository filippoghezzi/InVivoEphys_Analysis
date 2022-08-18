close all
clear 
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 
load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\NH_SingleUnitData.mat')
NH=data;


data = table;
%% Liad's data
load('LJB_PPCData.mat')
load('D:\InVivo_S1\Liad - S1Data\AllData\NewClassificationData.mat')
LJB=struct2table(unitT);

nLJB=height(LJB);
nanLJB=nan(nLJB,1);

LJBtagging=repmat(categorical(cellstr('SST')),nLJB,1);
LJBtagging(logical(LJB.Mut))=categorical(cellstr('SST;NrgKO'));
LJBresponsive=LJB.Responsive;
LJBresponsive(LJBresponsive==-1)=2;
LJBtag=LJB.Tagged;
LJBtag(LJBtag==-1)=2;

data.MouseID = [NH.MouseID;categorical(cellstr(LJB.AnimalId))];
data.Tagging = [NH.Tagging;LJBtagging];
data.brainArea = [NH.brainArea; repmat(categorical(cellstr('S1BF')),nLJB,1)];
data.state = [NH.state; repmat(categorical(cellstr('Urethane')),nLJB,1)];
data.suid = [NH.suid; LJB.UnitId];
data.Age = [NH.Age; LJB.Age];
data.Layer = [NH.Layer; LJB.Layer];
data.Depth = [NH.Depth; LJB.Depth];
data.wf = [NH.wf; nan(nLJB, size(NH.wf,2))];
data.filt_wf = [NH.filt_wf; nan(nLJB, size(NH.filt_wf,2))];
data.halfWidth = [NH.halfWidth; nanLJB];
data.troughPeakTime = [NH.troughPeakTime; nanLJB];
data.peakTroughRatio = [NH.peakTroughRatio; nanLJB];
data.endSlope = [NH.endSlope; nanLJB];
data.PSTHvisual = [NH.PSTHvisual; nan(nLJB,size(NH.PSTHvisual,2))];
data.PSTHwhisker =  [NH.PSTHvisual; nan(nLJB,size(NH.PSTHvisual,2))];
data.PSTHoptotagging = [NH.PSTHoptotagging; nan(nLJB,size(NH.PSTHvisual,2))];
data.PSTHvisualOpto = [NH.PSTHvisualOpto; nan(nLJB,size(NH.PSTHvisual,2))];
data.PSTHlaser = [NH.PSTHlaser; nan(nLJB,size(NH.PSTHvisual,2))];
data.responseVisual = [NH.responseVisual; nanLJB];
data.responseWhisker = [NH.responseWhisker; LJBresponsive];
data.responseTag = [NH.responseTag; LJBtag];
data.responseVisualOpto = [NH.responseVisualOpto; nanLJB];
data.responseLaser = [NH.responseLaser; nan(nLJB,2)];
data.reliabilityVisual = [NH.reliabilityVisual; nanLJB];
data.reliabilityVisual_K = [NH.reliabilityVisual_K; nanLJB];
data.fanoVisual = [NH.fanoVisual; nanLJB];
data.fanoVisualOpto = [NH.fanoVisualOpto; nanLJB];
data.fanoVisual_K = [NH.fanoVisual_K; nanLJB];
data.fanoWhisker = [NH.fanoWhisker; nanLJB];
data.fanoWhisker_K = [NH.fanoWhisker_K; nanLJB];
data.PPC = [NH.PPC;PPC_all]; 
data.vectorLength = [NH.vectorLength;vectorLength_all]; 
data.vectorAngle = [NH.vectorAngle;vectorAngle_all]; 
data.pValuePPC = [NH.pValuePPC;pValue_all]; 
data.rw_responsive = [NH.rw_responsive; nanLJB];
data.rw_firing = [NH.rw_firing; nanLJB];
data.rw_baselineFiring = [NH.rw_baselineFiring; nanLJB];
data.rw_rho = [NH.rw_rho; nan(nLJB,size(NH.rw_rho,2))];
data.baseline_firing = [NH.baseline_firing; LJB.MeanBaseline];
data.sb_spikeProb = [NH.sb_spikeProb; nanLJB];
data.sb_entrainedP = [NH.sb_entrainedP; nanLJB];
data.sb_firing = [NH.sb_firing; LJB.SpindleFiringRate];
data.PSTHvisual_K = [NH.PSTHvisual_K; nan(nLJB,size(NH.PSTHvisual_K,2))];
data.responseVisual_K = [NH.responseVisual_K; nanLJB];
data.PSTHwhisker_K = [NH.PSTHwhisker_K; nan(nLJB,size(NH.PSTHwhisker_K,2))];
data.responseWhisker_K = [NH.responseLaser; nan(nLJB,2)];
data.PPC_K = [NH.PPC_K;nan(nLJB,size(NH.PPC_K,2))]; 
data.vectorLength_K = [NH.vectorLength_K;nan(nLJB,size(NH.vectorLength_K,2))]; 
data.vectorAngle_K = [NH.vectorAngle_K;nan(nLJB,size(NH.vectorAngle_K,2))]; 
data.pValuePPC_K = [NH.pValuePPC_K;nan(nLJB,size(NH.pValuePPC_K,2))]; 
data.rw_responsive_K = [NH.rw_responsive_K; nanLJB];
data.rw_firing_K = [NH.rw_firing_K; nanLJB];
data.rw_baselineFiring_K = [NH.rw_baselineFiring_K; nanLJB];
data.rw_rho_K = [NH.rw_rho_K; nan(nLJB,size(NH.rw_rho_K,2))];
data.baseline_firing_K = [NH.baseline_firing_K; nanLJB];
data.sb_spikeProb_K = [NH.sb_spikeProb_K; nanLJB];
data.sb_entrainedP_K = [NH.sb_entrainedP_K; nanLJB];
data.sb_firing_K = [NH.sb_firing_K; nanLJB];

save('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\NH_LJB_SingleUnitData.mat', 'data', 'PSTHbins', 'rwFreq')


