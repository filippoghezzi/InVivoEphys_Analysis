function timeFrequencyAnalysis(evokedLFP)

fs=1000;
% f=(0:0.5:50);
t=(-1000:1:5100);
% figure
% for channel=1:size(evokedLFP,1)
%     for trial=1:size(evokedLFP,2)
%         p(trial,:,:)=pspectrum(reshape(evokedLFP(channel,trial,:),1,size(evokedLFP,3)),fs,'spectrogram','FrequencyLimits',[0,50],'TimeResolution',0.2);    
%     end
%     ps=reshape(mean(p,1),size(p,2),size(p,3));
%     subplot(4,8,channel)
%     imagesc(ps)
% 
% end

figure

for channel=1:size(evokedLFP,1)
    for trial=1:size(evokedLFP,2)
        filteredLFP=BandPassButerworthFilter(evokedLFP(channel,trial,:),1000);
        [cfs(trial,:,:),f] = cwt(reshape(filteredLFP,1,size(evokedLFP,3)),'morse',fs,'FrequencyLimits',[0,50],'VoicesPerOctave',32);

%         p_w(trial,:,:)=pspectrum(reshape(evokedLFP(channel,trial,:),1,size(evokedLFP,3)),fs,'spectrogram','FrequencyLimits',[0,50],'TimeResolution',0.2);    
    end
    ps_w=reshape(mean(cfs,1),size(cfs,2),size(cfs,3));
    subplot(4,1,channel)
    helperCWTTimeFreqPlot(ps_w,t,f,'surf')
 
    
    
    
end







end
