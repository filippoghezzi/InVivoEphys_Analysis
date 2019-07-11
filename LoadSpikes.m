function [spikes, templates, suid] = LoadSpikes(directory, ElectrodeMap)
% Load spikes in directory
% return a [spikeIds quality chan time] structure
% 0=noise, 1=MUA, 2=Good, 3=unsorted)
spikeTimes = readNPY(fullfile(directory,'spike_times.npy'));
spikeClusters = readNPY(fullfile(directory,'spike_clusters.npy'));
[cids, cgs]  = readClusterGroupsCSV(fullfile(directory,'cluster_groups.csv'));
templateId = readNPY(fullfile(directory,'spike_templates.npy'));
templates = readNPY(fullfile(directory,'templates.npy'));
[templateChan] = findBiggestTemplate(templates); % get biggest chan per template

[~,actChan] = ismember(templateChan+16,ElectrodeMap);

% templateId(templateId==0) = 1; % remove zero ones hopefully they coincide with noise templates


[~,Quality] = ismember(spikeClusters,cids);
[~,ChIndex] = ismember(templateId+1,1:length(actChan));
spikes = [spikeClusters cgs(Quality)' actChan((templateId+1)') spikeTimes];

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
        if ~(isempty(p))& p>pMax
            chMax = ch;
            pMax = p;
        end
    end
    
 templateChan(t) = chMax;   
end

templateChan = templateChan';
end