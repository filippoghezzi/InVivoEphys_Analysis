function MUA = getLayerMUA(s,stimulus,fs,window,binSize,channels,varargin)
% function MUA = getMUA(s,stimulus,ops)
%
% Build depth profile PSTH of MUA through all channels in Nchan. MUA is
% aligned to stimulus and build using makePSTHandRaster function (requires
% window and binSize)
    
    if ~isempty(varargin)
        artefactRemoval=varargin{1,1};
    else 
        artefactRemoval=[];
    end
        
    clusterID=s.cids(ismember(s.cch,channels));
    spiketimes=s.st(ismember(s.sclu,clusterID));

    [PSTH,PSTHbins,~]=makePSTHandRaster(spiketimes,stimulus,fs,window,binSize,artefactRemoval);


    % Output
    MUA.raw = PSTH;
    MUA.bins=PSTHbins;
end

