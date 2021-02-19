function MUA = getDepthMUA(s,stimulus,fs,window,binSize,Nchan)
% function MUA = getMUA(s,stimulus,ops)
%
% Build depth profile PSTH of MUA through all channels in Nchan. MUA is
% aligned to stimulus and build using makePSTHandRaster function (requires
% window and binSize)
        
    for channel=1:Nchan
        clusterID=s.cids(ismember(s.cch,channel));
        spiketimes=s.st(ismember(s.sclu,clusterID));
        
        [PSTH(channel,:),PSTHbins,~]=makePSTHandRaster(spiketimes,stimulus,fs,window,binSize,[]);

    end

        
    % Output
    MUA.raw = PSTH;
    MUA.bins=PSTHbins;
end

