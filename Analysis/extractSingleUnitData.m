close all
clear 
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 

recordings=readtable('V1_InVivo_SST;Ai32.csv');
recordings=recordings(~isnan(recordings.Sorting),:);
recID=unique(recordings.MouseID);

recFolder='D:\InVivo_V1';

mouseID=[];
optotagging=[];
suid=[];
suage=[];
sulayer=[];
wf=[];
filt_wf=[];
halfWidth=[];
troughPeakTime=[];
peakTroughRatio=[];
endSlope=[];
PSTHvisual=[];
PSTHoptotagging=[];
PSTHvisualOpto=[];
PSTHlaser=[];
responseVisual=[];
responseTag=[];
responseVisualOpto=[];
responseLaser=[];


PPC=[];
vectorLength=[];
vectorAngle=[];
pValuePPC=[];

% zscorecoherence=0;
coherence=[];

for i=1:numel(recID)
        thisTagging=recordings.Tagging(strcmp(recordings.MouseID,recID{i}),:);
        thisTagging=thisTagging{1};
        dir=fullfile(recFolder,recID{i});
        load(fullfile(dir,'rez.mat'),'rez')
        load(fullfile(dir,'spikes.mat'),'s')
        
        thisID=cell(numel(s.suid),1);
        thisID(:)={recID{i}};
        tagging=cell(numel(s.suid),1);
        tagging(:)={thisTagging};
        % Spike features
        mouseID=[mouseID;thisID];
        optotagging=[optotagging;tagging];
        suid=[suid;s.suid];
        suage=[suage;s.suage];
        sulayer=[sulayer;s.sulayer];
        halfWidth=[halfWidth;s.halfWidth];
        troughPeakTime=[troughPeakTime;s.troughPeakTime];
        peakTroughRatio=[peakTroughRatio;s.peakTroughRatio];
        endSlope=[endSlope;s.endSlope];
        wf=[wf;s.suWf];
        filt_wf=[filt_wf;s.filt_suWf];

        % PSTH
        PSTHvisual=[PSTHvisual;s.PSTHvisual];
        responseVisual=[responseVisual;s.response.visual];
        PSTHoptotagging=[PSTHoptotagging;s.PSTHoptotagging];
        responseTag=[responseTag;s.response.optotagging];
        PSTHvisualOpto=[PSTHvisualOpto;s.PSTHvisualOpto];
        responseVisualOpto=[responseVisualOpto;s.response.visualOpto];
        PSTHlaser=[PSTHlaser;s.PSTHlaser];
        responseLaser=[responseLaser;s.response.laser];

        % PPC
        PPC=[PPC;rez.ops.phaseLocking.PPC];
        vectorLength=[vectorLength;rez.ops.phaseLocking.vectorLength];
        vectorAngle=[vectorAngle;rez.ops.phaseLocking.vectorAngle];
        pValuePPC=[pValuePPC;rez.ops.phaseLocking.pValue];
        
        % Coherence
        coherenceFreq=rez.ops.coherence.freqs;
        [betaCoherence,maxCoherenceIdx]=max(rez.ops.coherence.curve(:,(coherenceFreq>10) & (coherenceFreq<30)),[],2);
        [lowCoherence,~]=max(rez.ops.coherence.curve(:,(coherenceFreq>0) & (coherenceFreq<5)),[],2);
        [gammaCoherence,~]=max(rez.ops.coherence.curve(:,(coherenceFreq>40) & (coherenceFreq<80)),[],2);
%         if zscorecoherence
%             betaCoherence=(betaCoherence-rez.ops.coherence.curve_shuff_mean(maxCoherenceIdx))./rez.ops.coherence.curve_suff_SD(maxCoherenceIdx);
%         end
        coherence=[coherence;lowCoherence,betaCoherence,gammaCoherence];
end

coherence(coherence==0)=NaN;
PSTHbins=s.PSTHbins;

% %Z-scoring SU
% Z_PSTHvisual=zscoreBaseline(PSTHvisual);
% Z_PSTHvisualOpto=zscoreBaseline(PSTHvisualOpto);
% Z_PSTHlaser=zscoreBaseline(PSTHlaser);
% Z_PSTHoptotagging=zscoreBaseline(PSTHoptotagging);

%% Make table
data=table;
data.MouseID=categorical(mouseID);
data.Tagging=categorical(optotagging);
data.suid=suid;
data.Age=suage;
data.Layer=sulayer;
data.wf=wf;
data.filt_wf=filt_wf;
data.halfWidth=halfWidth;
data.troughPeakTime=troughPeakTime;
data.peakTroughRatio=peakTroughRatio;
data.endSlope=endSlope;
data.PSTHvisual=PSTHvisual;
data.PSTHoptotagging=PSTHoptotagging;
data.PSTHvisualOpto=PSTHvisualOpto;
data.PSTHlaser=PSTHlaser;
data.responseVisual=responseVisual;
data.responseTag=responseTag;
data.responseVisualOpto=responseVisualOpto;
data.responseLaser=responseLaser;
data.PPC=PPC;
data.vectorLength=vectorLength;
data.vectorAngle=vectorAngle;
data.pValuePPC=pValuePPC;
data.coherence=coherence;

%% Save
save('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\SingleUnitData.mat', 'data', 'PSTHbins')