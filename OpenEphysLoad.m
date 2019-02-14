
% Electrode Map File Load
PolytrodeMapPath  = 'C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\A1x32_Map.mat';
% FourShankMapPath = 'C:\Electrophys\A32_4X8ElectrodeEphysMap.mat';

load(PolytrodeMapPath);
%  load(FourShankMapPath);

ElectrodeMap = reshape(ElectrodeMap,32,1);

data  = [];
j = 1;
for i = 17:48
   [tmp, timestamps, info] = load_open_ephys_data_faster(strcat('100_CH',int2str(i),'.continuous'));
   tmp = tmp(1:20:end);
   timestamps = timestamps(1:20:end);
   data(j,:) =  tmp';
   j = j+1; 
end

% data = data - median(data,1);
% tmp = [];
data = data(ElectrodeMap-16,:);
% j = 17;
% for i = 49:64
%    [tmp, timestamps, info] = load_open_ephys_data_faster(strcat('100_CH',int2str(i),'_2.continuous'));
%    data(j,:) =  tmp';
%    j = j+1;
% end


% 



%%
[events, timestamps2, info] = load_open_ephys_data_faster('all_channels.events');


% events = events(52:end);
% timestamps2 = timestamps2(52:end);

timestamps2 = timestamps2(events==2);

st = timestamps2(1:2:end);
  
et = timestamps2(2:2:end);
% 
% 


factor = 1; % the time factor (20 or 1 per msecond depending on sampling

stimT= [];
stims=[];
stimTraces = [];
stimTrials = [];
for ch = 1:size(data,1)
%     d = filtfilt(b,a,data(ch,:));
       d = data(ch,:);
    for i = 1:length(st)
        
        [~,s] = min(abs(timestamps-st(i)));
        [~,e] = min(abs(timestamps-et(i)));   
        if (ch ==1)
            stimT = [stimT;s e];
        end
%         if (e-s)>2000
%             continuei
%         end
        t = (max(1,s-500*factor):e+500);    
    %     stims(i,:) =  data(15,t);
%     bl = mean(d(:,s-500*factor:s),2) ;
%     stims(i,1:length(t)) =  d(t')-bl;
      stims(i,1:length(t)) =  d(t');
%     
    end
    stimTrials(ch,1:size(stims,1),1:size(stims,2)) = stims;
    stimTraces(ch,1:size(stims,2)) = mean(stims);
end 

save('Data.mat','stimTrials','stimTraces','stimT','data','ElectrodeMap','-v7.3');
%% MMN
% [events, timestamps2, info] = load_open_ephys_data_faster('all_channels.events');
%   st = timestamps2(1:2:end);
%   et = timestamps2(2:2:end);
% stimsGen = [stims1 stims2];
% st40 = st(stimsGen==40);
% et40 = et(stimsGen==40);
% 
% st10 = st(stimsGen==10);
% et10 = et(stimsGen==10);
% 
% 
% %     d = data(20,:);
% d = mean(data);
%     mx = max([length(st40) length(st10)]);
%     
%     for i = 1:mx
%     if (i<=length(st40))
%         [~,s] = min(abs(timestamps-st40(i)));
%         [~,e] = min(abs(timestamps-et40(i)));        
%         t = (max(1,s-500):e); 
%         bl = mean(d((max(1,s-500):s-1)));
%         stims40(i,1:length(t)) =  d(t')-bl;
%     end
%    d 
%      if (i<=length(st10))
%        [~,s] = min(abs(timestamps-st10(i)));
%         [~,e] = min(abs(timestamps-et10(i))); 
%         t = (max(1,s-500):e); 
%         bl = mean(d((max(1,s-500):s-1)));
%         stims10(i,1:length(t)) =  d(t')-bl;    
%      end
%         
%            
% 
%     
%     end
% 
