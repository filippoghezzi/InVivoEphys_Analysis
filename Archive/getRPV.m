clear 
close all
clc

load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\A1x32_Map.mat');
[spikeTimes,templates,suid]=LoadSpikes('C:\Users\Butt Lab\Documents\SpikeSorting\SC28',ElectrodeMap);

cluster_id=unique(spikeTimes(:,1));
allSpikeTimes=double(spikeTimes(:,4))/30000;

for i=1:length(cluster_id)
    cluster_spike_times = allSpikeTimes(spikeTimes(:,1)==cluster_id(i));
%     ax1=subplot(2,1,1);
%     ax2=subplot(2,1,2);
    [x,y,~,~]=myACG(cluster_spike_times,[],[]);
    RPV(i,1)=sum(y(x<0.002).*0.0005);
    
end
