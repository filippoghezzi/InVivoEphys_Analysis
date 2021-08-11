function rw=InVivo_dataProcessing_evokedLFP_retinalWave(ops,s,LFP,stim)
    
    fprintf(strcat('Processing light-evoked retinal waves...','\n'))

    LFP=LFP(:,ops.LFPwindow(1)*ops.fs:end,:);
    LFP=permute(LFP,[2,1,3]);
%     time=(0:ops.fs*ops.LFPwindow(2));
%     time=time(1:end)/ops.fs*1000;    
    
    %% Initialize variable  

    frequencyBins=(0:5:80)';
    PPC=nan(numel(s.suid),size(frequencyBins,1)-1,numel(stim));
    vectorLength=nan(numel(s.suid),size(frequencyBins,1)-1,numel(stim));
    vectorAngle=nan(numel(s.suid),size(frequencyBins,1)-1,numel(stim));
    pValue=nan(numel(s.suid),size(frequencyBins,1)-1,numel(stim));
    powerLFP=nan(numel(s.suid),size(frequencyBins,1)-1,numel(stim));
    singleUnitFiringFrequency_retinalWave=zeros(numel(s.suid),numel(stim));
    singleUnitFiringFrequency_baseline=zeros(numel(s.suid),numel(stim));
    spikeLFP_rho=zeros(numel(s.suid),size(frequencyBins,1)-1);
    
    %% Analyse light-induced retinal wave
    for frequency=1:size(frequencyBins,1)-1
        filterFrequency=[frequencyBins(frequency)+1, frequencyBins(frequency+1)];
        [b,a]=butter(2,filterFrequency/(ops.fs/2),'bandpass');
        filteredLFP=filtfilt(b,a,LFP);

        for singleUnit=1:numel(s.suid)
            singleUnitChannel=s.cch(s.cids==s.suid(singleUnit));
            selectedChannel=singleUnitChannel+2;
            if selectedChannel>32; selectedChannel=singleUnitChannel-2; end

            for trial=1:numel(stim)
                singleUnitSpikeTimes=double(s.st(s.sclu==s.suid(singleUnit)))-stim(trial);
                singleUnitSpikeTimes_retinalWave=singleUnitSpikeTimes(singleUnitSpikeTimes>0.4*ops.fs & singleUnitSpikeTimes<4.4*ops.fs);

                singleUnitFiringFrequency_retinalWave(singleUnit,trial)=numel(singleUnitSpikeTimes_retinalWave)/4;
                [PPC(singleUnit,frequency,trial), vectorAngle(singleUnit,frequency,trial), vectorLength(singleUnit,frequency,trial), pValue(singleUnit,frequency,trial), powerLFP(singleUnit,frequency,trial)]...
                    =getPPC(singleUnitSpikeTimes_retinalWave,filteredLFP(:,selectedChannel,trial)', 'Threshold', 1);
                
                singleUnitSpikeTimes_baseline=singleUnitSpikeTimes(singleUnitSpikeTimes>-4.4*ops.fs & singleUnitSpikeTimes<-0.4*ops.fs);
                singleUnitFiringFrequency_baseline(singleUnit,trial)=numel(singleUnitSpikeTimes_baseline)/4;
            end
        end
    end
    
    %% Get LFP power/single unit firing correlation
    for singleUnit=1:numel(s.suid)
        fr=singleUnitFiringFrequency_retinalWave(singleUnit,:)';
        p=squeeze(powerLFP(singleUnit,:,:))';
        spikeLFP_rho(singleUnit,:)=corr(fr,p,'Type','Spearman');
    end
    
    %% Output
    rw.freq=frequencyBins(2:end);
    rw.PPC=PPC;
    rw.vectorAngle=vectorAngle;
    rw.vectorLength=vectorLength;
    rw.pValue=pValue;
    rw.powerLFP=powerLFP;
    rw.spikeLFP_rho=spikeLFP_rho;
    rw.firingFreq=singleUnitFiringFrequency_retinalWave;
    rw.baselineFiringFreq=singleUnitFiringFrequency_baseline;
    
end
