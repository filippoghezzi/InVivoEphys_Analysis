function [spikes, templates, suid] = LoadSpikes(directory, ElectrodeMap)
% Load spikes in directory
% return a [spikeIds quality chan time] structure
% 0=noise, 1=MUA, 2=Good, 3=unsorted)
spikeTimes = readNPY(fullfile(directory,'spike_times.npy'));
spikeClusters = readNPY(fullfile(directory,'spike_clusters.npy'));
if isfile(fullfile(directory,'cluster_groups.csv'))
    [clustersIDs, clustersGroups]  = readClusterGroupsCSV(fullfile(directory,'cluster_groups.csv')); %clusterIDs corresponds to unique(spikeClusters)
elseif isfile(fullfile(directory,'cluster_group.tsv'))
    [clustersIDs, clustersGroups] = readClusterGroupsTSV(fullfile(directory,'cluster_group.tsv'));
end
templateId = readNPY(fullfile(directory,'spike_templates.npy'));
templates = readNPY(fullfile(directory,'templates.npy'));
[templateChan] = findBiggestTemplate(templates); % get biggest chan per template

[~,actChan] = ismember(templateChan+16,ElectrodeMap); %act stands for actual channels with 1 top channel and 32 bottom channel

% templateId(templateId==0) = 1; % remove zero ones hopefully they coincide with noise templates


[~,Quality] = ismember(spikeClusters,clustersIDs);
% [~,ChIndex] = ismember(templateId+1,1:length(actChan));
spikes = [spikeClusters clustersGroups(Quality)' actChan((templateId+1)') spikeTimes];


%Remove too big spikes taller 1% to remove noise
% amplitude=(readNPY(fullfile(directory,'amplitudes.npy')));
% [cdf,bin] = histcounts(amplitude,'Normalization','cdf');
% cutoff=bin(find(cdf>0.999,1,'first')); %in uV
% spikes(amplitude>cutoff,2)=0;


% Get sorted unit (su) ids
su = find(spikes(:,2)==2);
suid = unique(spikes(su,1));

templatesT = [];
for i = 1:length(suid)
    
    a = find(spikes(:,1)==suid(i),1,'first');
    templatesT(i,:) = templates(templateId(a)+1,:,templateChan(templateId(a)+1));    
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