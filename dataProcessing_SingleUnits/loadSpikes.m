function spikeStruct = loadSpikes(dir)
% function spikeStucture = loadSpikesNew(dir)
%
% Load spike data obtained from KS2 and phy2 in directory dir.
% Cluster groups: 0 -> noise; 1 -> mua; 2 -> single unit; 3 -> unsorted.

addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\npy-matlab')) 

%% Load spike times and clusters
    spikeTimes = readNPY(fullfile(dir,'spike_times.npy'));
    spikeClusters = readNPY(fullfile(dir,'spike_clusters.npy'));
%     spikeTemplates = readNPY(fullfile(dir,'spike_templates.npy'));

    clusterinfo=tdfread(fullfile(dir,'cluster_info.tsv'));
    clusterChannel=clusterinfo.ch';
    [clusterID,clusterGroup]=readClusterGroupsTSV(fullfile(dir,'cluster_group.tsv'));
%     clusterTemplate=readNPY(fullfile(dir,'templates.npy'));

    suid=clusterID(clusterGroup==2);

%% Remove noise
    noiseID=clusterID(clusterGroup==0);

    clusterChannel=clusterChannel(~ismember(clusterID,noiseID));
    clusterGroup=clusterGroup(~ismember(clusterID,noiseID));
    clusterID=clusterID(~ismember(clusterID,noiseID));

    spikeTimes=spikeTimes(~ismember(spikeClusters,noiseID));
%     spikeTemplates=spikeTemplates(~ismember(spikeClusters,noiseID));
    spikeClusters=spikeClusters(~ismember(spikeClusters,noiseID));

%% Templates
%     suWf=zeros(numel(suid),size(clusterTemplate,2));
%     for i=1:numel(suid)
%         templateIdx=unique(spikeTemplates(spikeClusters==suid(i)));
%         channel=clusterChannel(ismember(clusterID,suid(i)));
%         if numel(templateIdx) ~= 1
%             error('More than one average template associated with one cluster')
%         end
%         suWf(i,:)=clusterTemplate(templateIdx,:,channel);
%     end

%% Waveforms: only average waveform obtained from phy
    if exist(fullfile(dir,'su_Waveforms.npy'),'file')
        suWf=readNPY(fullfile(dir,'su_Waveforms.npy'));
        suWf_ids=readNPY(fullfile(dir,'su_Waveforms_ids.npy'));
        [suWf, suWf_ids] = findBestChannelFromWaveforms(suWf,suWf_ids);
    end
    clusterChannel(ismember(clusterID,suid))=suWf_ids;
    
%% Build output structure
    spikeStruct.st=spikeTimes;
    spikeStruct.sclu=spikeClusters;
    spikeStruct.cids=clusterID;
    spikeStruct.cgs=clusterGroup;
    spikeStruct.cch=clusterChannel+1; % Original clusterChannel starts  0:NchanTOT-1
    spikeStruct.suid=suid';
    spikeStruct.suWf=suWf;

end

function [bestWf,bestIds] = findBestChannelFromWaveforms(wf,ids)
    
    bestWf=zeros(size(wf,3),size(wf,1));
    bestIds=zeros(size(wf,3),1);
    
    for suIdx=1:size(wf,3)
        thisWf=wf(:,:,suIdx);
        thisWf=detrend(thisWf);
        thisWf=thisWf-mean(thisWf(1:40,:),1);
        [~,minimumIndex]=min(min(thisWf,[],1));
        bestWf(suIdx,:)=wf(:,minimumIndex,suIdx);
        bestIds(suIdx,1)=ids(minimumIndex,suIdx);
    end
end
