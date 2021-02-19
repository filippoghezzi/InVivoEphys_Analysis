function ops=get_VEPfeatures(eLFP,ops)


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
    
    
    
    
    
    
    
    %% Output
    ops.corticalResponseOnset=mean(t_response(ops.L4));
    ops.peakVEP=mean(peakVEP_V(ops.L4));
    ops.timeVEP=mean(peakVEP_t(ops.L4));

