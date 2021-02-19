function ops=analyseStimulusEvokedFieldPotential(ops,s,stim)
%
% Inputs: ops -> structure of recording info and results;
%         s -> structure of spike data;    


	ops.chanMap=(1:ops.NchanTOT)';
    ops.LFPwindow=[1 5]; %s
    ops.PSTHbinSize=0.003;
    ops.electrodeLength = 800; %um
    ops.electrodeSpacing = 25*10^-6; %um
    ops.verbose=1;  
    
    %% Process LFP, CSD and MUA; find L4
    eLFP = get_eLFP(ops.fbinary,stim.ledR,ops.fs,ops.NchanTOT,ops.LFPwindow,'LFP');
    CSD = getCSD(eLFP,ops.fs,ops.electrodeSpacing,ops.electrodeLength);
    MUA = getDepthMUA(s,stim.ledR,ops.fs,ops.LFPwindow,ops.PSTHbinSize,ops.NchanTOT);
    
%     if ~isfield(ops,'L4')
        ops=findL4(eLFP,CSD,MUA,ops);
%     end    
    
    %% Analyse LFP
    ops=get_VEPfeatures(eLFP,ops);

    
    %% MUA by cortical layer
    ops.L2=1:min(ops.L4)-1;
    ops.L5=max(ops.L4)+1:ops.NchanTOT;
    
    %Select index for plotting a single trace per layer
    L4Idx=ops.L4(2);
    L2Idx=min(ops.L4)-4;
    while L2Idx<1
        L2Idx=L2Idx+1;
        if L2Idx==L4Idx
            error('L4 index is too low')
        end
    end
    L5Idx=max(ops.L4)+6;
    while L5Idx>ops.NchanTOT
        L5Idx=L5Idx-1;
        if L5Idx==L4Idx
            error('L4 index is too high')
        end
    end
    
    eMUA_trace=get_eLFP(ops.fbinary,stim.ledR(10),ops.fs,ops.NchanTOT,ops.LFPwindow,'MUA');

    MUA(:,1)=getLayerMUA(s,stim.ledR,ops.fs,ops.LFPwindow,ops.PSTHbinSize,ops.L2);
    MUA(:,2)=getLayerMUA(s,stim.ledR,ops.fs,ops.LFPwindow,ops.PSTHbinSize,ops.L4);
    MUA(:,3)=getLayerMUA(s,stim.ledR,ops.fs,ops.LFPwindow,ops.PSTHbinSize,ops.L5);
    
    plotLFPbyLayer(eLFP([L2Idx,L4Idx,L5Idx],:,:),eMUA_trace([L2Idx,L4Idx,L5Idx],:),MUA,ops); 
    
    %% Time frequency analysis L4
    getSpectrogram(squeeze(eLFP(L4Idx,:,:)),ops)
    
    
end