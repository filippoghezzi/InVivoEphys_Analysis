function [LFPfeatures,ops]=InVivo_dataProcessing_evokedLFP_getFeaturesLFP(eLFP,ops)


    meanLFP=mean(eLFP,3);
    time=(-ops.fs*ops.LFPwindow(1):ops.fs*ops.LFPwindow(2));
    time=time/ops.fs*1000;

    %Calculate cortical response
    baselineWindow=time<0 & time>=-100;
    conditionWindow=time>5 & time<=1000;
    
    for channel=1:size(meanLFP,1)         
        meanBaseline=mean(meanLFP(channel,baselineWindow),2);
        stdBaseline=std(meanLFP(channel,baselineWindow),[],2);
        cutoff=[meanBaseline-3*stdBaseline,meanBaseline+3*stdBaseline];
        
        responseOnsetIdx = find(meanLFP(channel,conditionWindow)<cutoff(1),1,'first')+find(time==5);
        if isempty(responseOnsetIdx)
            t_response(channel,1)=NaN;
        else
            t_response(channel,1)=time(responseOnsetIdx);
        end
    end
    
    
   [peakVEP_V,peakVEP_Idx]=min(meanLFP(:,conditionWindow),[],2);
   peakVEP_t=time(peakVEP_Idx+find(time==5))';
    
    
   [~,L4idx]=min(peakVEP_t(ops.L4));
   L4idx = L4idx+min(ops.L4)-1;
    
    
    
    %% Output
    ops.L4best = L4idx;
    LFPfeatures.corticalResponseOnset=t_response(L4idx);
    LFPfeatures.peakVEP=peakVEP_V(L4idx);
    LFPfeatures.timeVEP=peakVEP_t(L4idx);

