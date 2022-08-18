close all
clear 
clc

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis')) 

recordings=readtable('G:\Nesta_Spikesorting\NH_InVivo_RED.csv');
% recordings=readtable('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\V1_InVivo.csv');
recordings=recordings(~isnan(recordings.Sorting),:);
recID=unique(recordings.MouseID);
% recID={'NK58','NSC1','NSC2','NSC4','NSC5','SC55','SC56','SC91','SC92b','SC93b','SC107l','SC107u'};
% recFolder='D:\InVivo_V1';
 
mouseID=[];
optotagging=[];
brainArea=[];
state=[];

suid=[];
suage=[];
sulayer=[];
suDepth=[];
wf=[];
filt_wf=[];
halfWidth=[];
troughPeakTime=[];
peakTroughRatio=[];
endSlope=[];
PSTHvisual=[];
PSTHwhisker=[];
PSTHoptotagging=[];
PSTHvisualOpto=[];
PSTHlaser=[];
responseVisual=[];
responseWhisker=[];
responseTag=[];
responseVisualOpto=[];
responseLaser=[];
PSTHvisualK=[];
responseVisualK=[];
PSTHwhiskerK=[];
responseWhiskerK=[];
reliabilityVisual=[];
reliabilityVisual_K=[];
fanoVisual=[];
fanoVisual_K=[];
fanoVisualOpto=[];
fanoWhisker=[];
fanoWhisker_K=[];

baseline_PPC=[];
baseline_vectorLength=[];
baseline_vectorAngle=[];
baseline_pValuePPC=[];

rw_responsive=[];
rw_firing=[];
rw_baselineFiring=[];
rw_rho=[];
baseline_firing=[];

coherence=[];

baseline_PPC_K=[];
baseline_vectorLength_K=[];
baseline_vectorAngle_K=[];
baseline_pValuePPC_K=[];

rw_responsive_K=[];
rw_firing_K=[];
rw_baselineFiring_K=[];
rw_rho_K=[];
baseline_firing_K=[];

sb_spikeProb=[];
sb_entrainedP=[];
sb_firing=[];
sb_spikeProb_K=[];
sb_entrainedP_K=[];
sb_firing_K=[];

coherence_K=[];

for i=1:numel(recID)
        fprintf(strcat('Looking into ...',recID{i},'\n'))

        thisTagging=recordings.Tagging(strcmp(recordings.MouseID,recID{i}),:);
        thisTagging=thisTagging{1};
        recFolder=recordings.Folder(strcmp(recordings.MouseID,recID{i}),:);
        recFolder=recFolder{1};
        Area=recordings.BrainArea(strcmp(recordings.MouseID,recID{i}),:);
        Area=Area{1};
        State=recordings.State(strcmp(recordings.MouseID,recID{i}),:);
        State=State{1};
        
        dir=fullfile(recFolder,recID{i});
        load(fullfile(dir,'processData.mat'),'results')
%         load(fullfile(dir,'spikes.mat'),'s')
        
        thisID=cell(numel(results.s.suid),1);
        thisID(:)={recID{i}};
        tagging=cell(numel(results.s.suid),1);
        tagging(:)={thisTagging};
        thisBrainArea=cell(numel(results.s.suid),1);
        thisBrainArea(:)={Area};
        thisState=cell(numel(results.s.suid),1);
        thisState(:)={State};  
        
        %% Spike features
        mouseID=[mouseID;thisID];
        optotagging=[optotagging;tagging];
        brainArea=[brainArea;thisBrainArea];
        state=[state;thisState];

        suid=[suid;results.s.suid];
        suage=[suage;results.s.suage];
        sulayer=[sulayer;results.s.sulayer];
        suDepth=[suDepth;results.s.suDepth];
        halfWidth=[halfWidth;results.s.halfWidth];
        troughPeakTime=[troughPeakTime;results.s.troughPeakTime];
        peakTroughRatio=[peakTroughRatio;results.s.peakTroughRatio];
        endSlope=[endSlope;results.s.endSlope];
        wf=[wf;results.s.suWf];
        filt_wf=[filt_wf;results.s.filt_suWf];

        %% PSTH        
        tmpPSTHvisual=nan(numel(results.s.suid),2000);
        tmpresponseVisual=nan(numel(results.s.suid),1);
        tmpPSTHwhisker=nan(numel(results.s.suid),2000);
        tmpresponseWhisker=nan(numel(results.s.suid),1);        
        tmpPSTHoptotagging=nan(numel(results.s.suid),2000);
        tmpresponseTag=nan(numel(results.s.suid),1);
        tmpPSTHvisualOpto=nan(numel(results.s.suid),2000);
        tmpresponseVisualOpto=nan(numel(results.s.suid),1);
        tmpPSTHlaser=nan(numel(results.s.suid),2000);
        tmpresponseLaser=nan(numel(results.s.suid),2);
        tmpPSTHvisual_K=nan(numel(results.s.suid),2000);
        tmpresponseVisualK=nan(numel(results.s.suid),1);
        tmpPSTHwhisker_K=nan(numel(results.s.suid),2000);
        tmpresponseWhisker_K=nan(numel(results.s.suid),1);
        tmpreliabilityVisual=nan(numel(results.s.suid),1);
        tmpreliabilityVisual_K=nan(numel(results.s.suid),1);
        tmpFanoVisual=nan(numel(results.s.suid),1);
        tmpFanoVisual_K=nan(numel(results.s.suid),1);
        tmpFanoVisualOpto=nan(numel(results.s.suid),1);
        tmpFanoWhisker=nan(numel(results.s.suid),1);
        tmpFanoWhisker_K=nan(numel(results.s.suid),1);

        
        if isfield(results.s,'PSTHvisual') && ~isempty(results.s.PSTHvisual)
            tmpPSTHvisual=results.s.PSTHvisual;
            tmpresponseVisual=results.s.response.visual;
            tmpreliabilityVisual=results.s.reliability.visual';
            tmpFanoVisual=results.s.fano.visual';
        end        
        if isfield(results.s,'PSTHwhisker') && ~isempty(results.s.PSTHwhisker)
            tmpPSTHwhisker=results.s.PSTHwhisker;
            tmpresponseWhisker=results.s.response.whisker;
            tmpFanoWhisker=results.s.fano.whisker';
        end      
        if isfield(results.s,'PSTHoptotagging') && ~isempty(results.s.PSTHoptotagging)
            tmpPSTHoptotagging=results.s.PSTHoptotagging;
            tmpresponseTag=results.s.response.optotagging;
        end    
        if isfield(results.s,'PSTHvisualOpto') && ~isempty(results.s.PSTHvisualOpto)
            tmpPSTHvisualOpto=results.s.PSTHvisualOpto;
            tmpresponseVisualOpto=results.s.response.visualOpto;
            tmpFanoVisualOpto=results.s.fano.visualOpto';
        end
        if isfield(results.s,'PSTHlaser') && ~isempty(results.s.PSTHlaser)
            tmpPSTHlaser=results.s.PSTHlaser;
            tmpresponseLaser=results.s.response.laser;
        end
        if isfield(results.s,'PSTHvisual_K') && ~isempty(results.s.PSTHvisual_K)
            tmpPSTHvisual_K=results.s.PSTHvisual_K;
            tmpresponseVisualK=results.s.response.visual_K;
            tmpreliabilityVisual_K=results.s.reliability.visual_K';
            tmpFanoVisual_K=results.s.fano.visual_K';
        end
        if isfield(results.s,'PSTHwhisker_K') && ~isempty(results.s.PSTHwhisker_K)
            tmpPSTHwhisker_K=results.s.PSTHwhisker_K;
            tmpresponseWhisker_K=results.s.response.whisker_K;
            tmpFanoWhisker_K=results.s.fano.whisker_K';
        end
        
        PSTHvisual=[PSTHvisual;tmpPSTHvisual];
        responseVisual=[responseVisual;tmpresponseVisual];        
        PSTHwhisker=[PSTHwhisker;tmpPSTHwhisker];
        responseWhisker=[responseWhisker;tmpresponseWhisker];        
        PSTHoptotagging=[PSTHoptotagging;tmpPSTHoptotagging];
        responseTag=[responseTag;tmpresponseTag];
        PSTHvisualOpto=[PSTHvisualOpto;tmpPSTHvisualOpto];
        responseVisualOpto=[responseVisualOpto;tmpresponseVisualOpto];
        PSTHlaser=[PSTHlaser;tmpPSTHlaser];
        responseLaser=[responseLaser;tmpresponseLaser];
        PSTHvisualK=[PSTHvisualK;tmpPSTHvisual_K];
        responseVisualK=[responseVisualK;tmpresponseVisualK];
        PSTHwhiskerK=[PSTHwhiskerK;tmpPSTHwhisker_K];
        responseWhiskerK=[responseWhiskerK;tmpresponseWhisker_K];  
        reliabilityVisual=[reliabilityVisual;tmpreliabilityVisual];
        reliabilityVisual_K=[reliabilityVisual_K;tmpreliabilityVisual_K];
        fanoVisual=[fanoVisual;tmpFanoVisual];
        fanoVisualOpto=[fanoVisualOpto;tmpFanoVisualOpto];
        fanoVisual_K=[fanoVisual_K;tmpFanoVisual_K];
        fanoWhisker=[fanoWhisker;tmpFanoWhisker];
        fanoWhisker_K=[fanoWhisker_K;tmpFanoWhisker_K];

        %% Control    
        % PPC
        baseline_PPC=[baseline_PPC;results.baseline.phaseLocking.PPC];
        baseline_vectorLength=[baseline_vectorLength;results.baseline.phaseLocking.vectorLength];
        baseline_vectorAngle=[baseline_vectorAngle;results.baseline.phaseLocking.vectorAngle];
        baseline_pValuePPC=[baseline_pValuePPC;results.baseline.phaseLocking.pValue];
        
%         % Coherence
%         coherenceFreq=results.baseline.coherence.freqs;
%         [betaCoherence,maxCoherenceIdx]=max(results.baseline.coherence.curve(:,(coherenceFreq>10) & (coherenceFreq<30)),[],2);
%         [lowCoherence,~]=max(results.baseline.coherence.curve(:,(coherenceFreq>0) & (coherenceFreq<5)),[],2);
%         [gammaCoherence,~]=max(results.baseline.coherence.curve(:,(coherenceFreq>40) & (coherenceFreq<80)),[],2);
%         coherence=[coherence;lowCoherence,betaCoherence,gammaCoherence];
        
        % RW
        if isfield(results.evokedLFP,'rw')
            rw_p=[];
            responsive=[];
            for unit=1:numel(results.s.suid)
                rw_p=signrank(results.evokedLFP.rw.baselineFiringFreq(unit,:)',results.evokedLFP.rw.firingFreq(unit,:)');
                meanRWfiring=mean(results.evokedLFP.rw.firingFreq(unit,:));
                meanBaselineFiring=mean(results.evokedLFP.rw.baselineFiringFreq(unit,:),2);
                
                if (rw_p<=0.05) && (meanRWfiring>meanBaselineFiring)
                    responsive(unit,1) = 1;
                else
                    responsive(unit,1) = 0;
                end
            end
%             responsive=zeros(numel(results.s.suid),1);
%             responsive(rw_p<=0.05 & mean(results.evokedLFP.rw.firingFreq,2)>mean(results.evokedLFP.rw.baselineFiringFreq,2))=1;

            rw_responsive=[rw_responsive;responsive];
            rw_firing=[rw_firing;mean(results.evokedLFP.rw.firingFreq,2)];
            rw_baselineFiring=[rw_baselineFiring;mean(results.evokedLFP.rw.baselineFiringFreq,2)];
            rw_rho=[rw_rho;results.evokedLFP.rw.spikeLFP_rho];
        else
            rw_responsive=[rw_responsive;zeros(numel(results.s.suid),1)];
            rw_firing=[rw_firing;zeros(numel(results.s.suid),1)];
            rw_baselineFiring=[rw_baselineFiring;zeros(numel(results.s.suid),1)];
            rw_rho=[rw_rho;zeros(numel(results.s.suid),16)];
        end
        
        baseline_firing=[baseline_firing;results.baseline.singleUnitFiringFrequency];
        
        % SB
        if size(results.baseline.spindleBurst.spikeProbability,1)==1
            sb_spikeProb=[sb_spikeProb;results.baseline.spindleBurst.spikeProbability'];
            sb_entrainedP=[sb_entrainedP;results.baseline.spindleBurst.unitEntrained_p'];
            sb_firing=[sb_firing;results.baseline.spindleBurst.firingFrequency'];
        else
            sb_spikeProb=[sb_spikeProb;results.baseline.spindleBurst.spikeProbability];
            sb_entrainedP=[sb_entrainedP;results.baseline.spindleBurst.unitEntrained_p];
            sb_firing=[sb_firing;results.baseline.spindleBurst.firingFrequency];
        end
        
        %% Chemogenetic
    if isfield(results,'SalB')
    if isfield(results.SalB.baseline,'durationBaseline')
       % PPC
            baseline_PPC_K=[baseline_PPC_K;results.SalB.baseline.phaseLocking.PPC];
            baseline_vectorLength_K=[baseline_vectorLength_K;results.SalB.baseline.phaseLocking.vectorLength];
            baseline_vectorAngle_K=[baseline_vectorAngle_K;results.SalB.baseline.phaseLocking.vectorAngle];
            baseline_pValuePPC_K=[baseline_pValuePPC_K;results.SalB.baseline.phaseLocking.pValue];

    %         % Coherence
    %         [betaCoherence_K,maxCoherenceIdx_K]=max(results.SalB.baseline.coherence.curve(:,(coherenceFreq>10) & (coherenceFreq<30)),[],2);
    %         [lowCoherence_K,~]=max(results.SalB.baseline.coherence.curve(:,(coherenceFreq>0) & (coherenceFreq<5)),[],2);
    %         [gammaCoherence_K,~]=max(results.SalB.baseline.coherence.curve(:,(coherenceFreq>40) & (coherenceFreq<80)),[],2);
    %         coherence_K=[coherence_K;lowCoherence_K,betaCoherence_K,gammaCoherence_K];

            % RW
            if isfield(results.SalB.evokedLFP,'rw')
                rw_p_K=[];
                responsive_K=[];
                for unit=1:numel(results.s.suid)
                    rw_p_K=signrank(results.SalB.evokedLFP.rw.baselineFiringFreq(unit,:)',results.SalB.evokedLFP.rw.firingFreq(unit,:)');
                    meanRWfiring_K=mean(results.SalB.evokedLFP.rw.firingFreq(unit,:));
                    meanBaselineFiring_K=mean(results.SalB.evokedLFP.rw.baselineFiringFreq(unit,:),2);

                    if (rw_p_K<=0.05) && (meanRWfiring_K>meanBaselineFiring_K)
                        responsive_K(unit,1) = 1;
                    else
                        responsive_K(unit,1) = 0;
                    end
                end
                rw_responsive_K=[rw_responsive_K;responsive_K];
                rw_firing_K=[rw_firing_K;mean(results.SalB.evokedLFP.rw.firingFreq,2)];
                rw_baselineFiring_K=[rw_baselineFiring_K;mean(results.SalB.evokedLFP.rw.baselineFiringFreq,2)];
                rw_rho_K=[rw_rho_K;results.SalB.evokedLFP.rw.spikeLFP_rho];
            else
                rw_responsive_K=[rw_responsive_K;zeros(numel(results.s.suid),1)];
                rw_firing_K=[rw_firing_K;zeros(numel(results.s.suid),1)];
                rw_baselineFiring_K=[rw_baselineFiring_K;zeros(numel(results.s.suid),1)];
                rw_rho_K=[rw_rho_K;zeros(numel(results.s.suid),16)];
            end

            baseline_firing_K=[baseline_firing_K;results.SalB.baseline.singleUnitFiringFrequency];
            if size(results.SalB.baseline.spindleBurst.spikeProbability,1)==1
                sb_spikeProb_K=[sb_spikeProb_K;results.SalB.baseline.spindleBurst.spikeProbability'];
                sb_entrainedP_K=[sb_entrainedP_K;results.SalB.baseline.spindleBurst.unitEntrained_p'];
                sb_firing_K=[sb_firing_K;results.SalB.baseline.spindleBurst.firingFrequency'];
            else
                sb_spikeProb_K=[sb_spikeProb_K;results.SalB.baseline.spindleBurst.spikeProbability];
                sb_entrainedP_K=[sb_entrainedP_K;results.SalB.baseline.spindleBurst.unitEntrained_p];
                sb_firing_K=[sb_firing_K;results.SalB.baseline.spindleBurst.firingFrequency];
            end
            
            
    else
        baseline_PPC_K=[baseline_PPC_K;nan(numel(results.s.suid),3)];
        baseline_vectorLength_K=[baseline_vectorLength_K;nan(numel(results.s.suid),3)];
        baseline_vectorAngle_K=[baseline_vectorAngle_K;nan(numel(results.s.suid),3)];
        baseline_pValuePPC_K=[baseline_pValuePPC_K;nan(numel(results.s.suid),3)];
        rw_responsive_K=[rw_responsive_K;nan(numel(results.s.suid),1)];
        rw_firing_K=[rw_firing_K;nan(numel(results.s.suid),1)];
        rw_baselineFiring_K=[rw_baselineFiring_K;nan(numel(results.s.suid),1)];
        rw_rho_K=[rw_rho_K;nan(numel(results.s.suid),16)];
        baseline_firing_K=[baseline_firing_K;nan(numel(results.s.suid),1)];
        sb_spikeProb_K=[sb_spikeProb_K;nan(numel(results.s.suid),1)];
        sb_entrainedP_K=[sb_entrainedP_K;nan(numel(results.s.suid),1)];
        sb_firing_K=[sb_firing_K;nan(numel(results.s.suid),1)];
    end
        else
        baseline_PPC_K=[baseline_PPC_K;nan(numel(results.s.suid),3)];
        baseline_vectorLength_K=[baseline_vectorLength_K;nan(numel(results.s.suid),3)];
        baseline_vectorAngle_K=[baseline_vectorAngle_K;nan(numel(results.s.suid),3)];
        baseline_pValuePPC_K=[baseline_pValuePPC_K;nan(numel(results.s.suid),3)];
        rw_responsive_K=[rw_responsive_K;nan(numel(results.s.suid),1)];
        rw_firing_K=[rw_firing_K;nan(numel(results.s.suid),1)];
        rw_baselineFiring_K=[rw_baselineFiring_K;nan(numel(results.s.suid),1)];
        rw_rho_K=[rw_rho_K;nan(numel(results.s.suid),16)];
        baseline_firing_K=[baseline_firing_K;nan(numel(results.s.suid),1)];
        sb_spikeProb_K=[sb_spikeProb_K;nan(numel(results.s.suid),1)];
        sb_entrainedP_K=[sb_entrainedP_K;nan(numel(results.s.suid),1)];
        sb_firing_K=[sb_firing_K;nan(numel(results.s.suid),1)];
    end
end

% coherence(coherence==0)=NaN;
% coherence_K(coherence_K==0)=NaN;

PSTHbins=results.s.PSTHbins;
rwFreq=results.evokedLFP.rw.freq;
% %Z-scoring SU
% Z_PSTHvisual=zscoreBaseline(PSTHvisual);
% Z_PSTHvisualOpto=zscoreBaseline(PSTHvisualOpto);
% Z_PSTHlaser=zscoreBaseline(PSTHlaser);
% Z_PSTHoptotagging=zscoreBaseline(PSTHoptotagging);

%% Make table
data=table;
data.MouseID=categorical(mouseID);
data.Tagging=categorical(optotagging);
data.brainArea=categorical(cellstr(brainArea));
data.state=categorical(cellstr(state));
data.suid=suid;
data.Age=suage;
data.Layer=sulayer;
data.Depth=suDepth;
data.wf=wf;
data.filt_wf=filt_wf;
data.halfWidth=halfWidth;
data.troughPeakTime=troughPeakTime;
data.peakTroughRatio=peakTroughRatio;
data.endSlope=endSlope;
data.PSTHvisual=PSTHvisual;
data.PSTHwhisker=PSTHwhisker;
data.PSTHoptotagging=PSTHoptotagging;
data.PSTHvisualOpto=PSTHvisualOpto;
data.PSTHlaser=PSTHlaser;
data.responseVisual=responseVisual;
data.responseWhisker=responseWhisker;
data.responseTag=responseTag;
data.responseVisualOpto=responseVisualOpto;
data.responseLaser=responseLaser;
data.reliabilityVisual=reliabilityVisual;
data.reliabilityVisual_K=reliabilityVisual_K;
data.fanoVisual=fanoVisual;
data.fanoVisualOpto=fanoVisualOpto;
data.fanoVisual_K=fanoVisual_K;
data.fanoWhisker=fanoWhisker;
data.fanoWhisker_K=fanoWhisker_K;

data.PPC=baseline_PPC;
data.vectorLength=baseline_vectorLength;
data.vectorAngle=baseline_vectorAngle;
data.pValuePPC=baseline_pValuePPC;
% data.coherence=coherence;
data.rw_responsive=rw_responsive;
data.rw_firing=rw_firing;
data.rw_baselineFiring=rw_baselineFiring;
data.rw_rho=rw_rho;
data.baseline_firing=baseline_firing;
data.sb_spikeProb=sb_spikeProb;
data.sb_entrainedP=sb_entrainedP;
data.sb_firing=sb_firing;

data.PSTHvisual_K=PSTHvisualK;
data.responseVisual_K=responseVisualK;
data.PSTHwhisker_K=PSTHwhiskerK;
data.responseWhisker_K=responseWhiskerK;
data.PPC_K=baseline_PPC_K;
data.vectorLength_K=baseline_vectorLength_K;
data.vectorAngle_K=baseline_vectorAngle_K;
data.pValuePPC_K=baseline_pValuePPC_K;
% data.coherence_K=coherence_K;
data.rw_responsive_K=rw_responsive_K;
data.rw_firing_K=rw_firing_K;
data.rw_baselineFiring_K=rw_baselineFiring_K;
data.rw_rho_K=rw_rho_K;
data.baseline_firing_K=baseline_firing_K;
data.sb_spikeProb_K=sb_spikeProb_K;
data.sb_entrainedP_K=sb_entrainedP_K;
data.sb_firing_K=sb_firing_K;

%% Save
save('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\NH_SingleUnitData.mat', 'data', 'PSTHbins', 'rwFreq')