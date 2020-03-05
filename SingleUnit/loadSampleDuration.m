function tab=loadSampleDuration(tab,spikeFolder)

    
    tab.SampleLength(1,1)=0;
    tab.startSample(1,1)=0;
    tab.endSample(1,1)=0;
    for i=1:max(unique(tab.RecID))
    % for i=1:8
        if ~isempty(tab(tab.RecID==i,:))
            load(fullfile(spikeFolder,tab(tab.RecID==i,:).MouseID{1},'SampleDuration.mat'))
            tab(tab.RecID==i,:).SampleLength=samplesToSave;
            for j=1:height(tab(tab.RecID==i,:))
                if j==1
                    tab(tab.RecID==i,:).startSample(j)=1;
                    tab(tab.RecID==i,:).endSample(j)=tab(tab.RecID==i,:).SampleLength(j);
                else
                    tab(tab.RecID==i,:).startSample(j)=tab(tab.RecID==i,:).endSample(j-1)+1;
                    tab(tab.RecID==i,:).endSample(j)=tab(tab.RecID==i,:).SampleLength(j)+tab(tab.RecID==i,:).startSample(j);
                end
            end
        end
    end
    tab = removevars(tab,{'SampleLength','endSample'});

end