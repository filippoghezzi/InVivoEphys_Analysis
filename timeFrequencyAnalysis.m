function timeFrequencyAnalysis(evokedLFP)

fs=1000;
t=(-1000:1:5100);


    for trial=1:size(evokedLFP,1)
        filteredLFP=BandPassButerworthFilter(evokedLFP(trial,:),1000);
        [cfs(trial,:,:),f] = cwt(filteredLFP,'morse',fs,'FrequencyLimits',[0,50],'VoicesPerOctave',32);
    end
    
    
    %%%%%%%% STILL TO OPTIMIZE THE NORMALIZATION OF THE
    %%%%%%%% CSD%%%%%%%%%%%%%%%%%%%%%%%%
    cfs=abs(cfs).^2;
    ps_w=squeeze(mean(cfs,1));
    ps_w=ps_w./mean(ps_w(:,1:100),2);
    
    imagesc(t,f,ps_w)
    set(gca,'Ydir','normal')
    xlim([-1000 5000]);
    ylim([0 50])
%     helperCWTTimeFreqPlot(ps_w,t,f,'image')
 
   







end
