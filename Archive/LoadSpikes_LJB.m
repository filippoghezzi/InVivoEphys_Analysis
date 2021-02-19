function [spikes, templates, suid] = LoadSpikes(directory, ElectrodeMap)
% Load spikes in directory
% Inputs: 
%    directory -> where kilosort and phy data are stored;
%    ElectrodeMap -> map of electrode used.
% Outputs:
%   spikes -> matrix of the form [spikeCluster clusterGroup channel spikeTime];
%   clusterGroups is 0=noise, 1=MUA, 2=Good, 3=unsorted
%   templates ->
%   suid -> single unit IDs.

%Spikes and clusters
spikeTimes = readNPY(fullfile(directory,'spike_times.npy'));
spikeClusters = readNPY(fullfile(directory,'spike_clusters.npy'));
if isfile(fullfile(directory,'cluster_groups.csv'))
    [clustersIDs, clustersGroups]  = readClusterGroupsCSV(fullfile(directory,'cluster_groups.csv')); %clusterIDs corresponds to unique(spikeClusters)
elseif isfile(fullfile(directory,'cluster_group.tsv'))
    [clustersIDs, clustersGroups] = readClusterGroupsTSV(fullfile(directory,'cluster_group.tsv'));
end
[~,Quality] = ismember(spikeClusters,clustersIDs); %Take index of clusterId when spikeClusters is represented in clustersIDs; Quality is a vector (nspike,1)
clusterGroups=clustersGroups(Quality)';

%Templates
spikeTemplates = readNPY(fullfile(directory,'spike_templates.npy'));
templates = readNPY(fullfile(directory,'templates.npy'));
[templateChannel] = findBiggestTemplate(templates); % get biggest chan per template
[~,actChan] = ismember(templateChannel+16,ElectrodeMap); %act stands for actual channels with 1 top channel and 32 bottom channel
channel=actChan((spikeTemplates+1)');


spikes = [spikeClusters clusterGroups channel spikeTimes]; %OUTPUT


% Get sorted unit (su) ids
su = find(spikes(:,2)==2);
suid = unique(spikes(su,1));

templatesT = [];
for i = 1:length(suid)
    
    a = find(spikes(:,1)==suid(i),1,'first');
    templatesT(i,:) = templates(spikeTemplates(a)+1,:,templateChannel(spikeTemplates(a)+1));    
end

templates = templatesT;

end


function [templateChan] = findBiggestTemplate(templates)
%find the largest amplitude

% run through templates and pick peaks
for t = 1:size(templates,1)
    temp = squeeze(templates(t,:,:));
    
    pMax = 0;
    chMax = 0;
    for ch = 1:size(templates,3)
        [p] = findpeaks(-temp(:,ch),'MinPeakProminence',0.05);
        p = max(p);
        if ~(isempty(p))&& p>pMax
            chMax = ch;
            pMax = p;
        end
    end
    
 templateChan(t) = chMax;   
end

templateChan = templateChan';
end