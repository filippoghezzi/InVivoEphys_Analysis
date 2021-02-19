function removeArtefact(dir,s)
    
    fid=fopen(fullfile(dir,'RawData.dat'),'r');
    frewind(fid);
%     LFPwindowDuration=ops.LFPwindow(2)+ops.LFPwindow(1);
%     LFPdata=zeros(ops.Nchan,LFPwindowDuration*ops.fs,size(stimuli,1),'int16');
    j=0;
    for i=1:numel(s.st)
        % Align file position to stimulus; remember binary file is in the
        % format Ch1,T1,Ch2,T1,...ChN,T1,Ch1,T2,...ChN,TN. First index is
        % 0, align onto (T-1)*2*Nchan with T stimulus onset sample. 
        
        T=s.st(i)-(0.0015*30000);
        offset = (2*30)*(T-1);
        
        fseek(fid, offset, 'bof');
        tmpData = fread(fid, [32 0.004*30000], '*int16');
        channel=s.cch(s.sclu(i)==s.cids);
        spikeTrace=tmpData(channel,:);
        
        if s.sclu(i)==s.suid(1)
            j=j+1;
            template(j,:)=spikeTrace;
        end
        
        
        LFPdata(:,:,i)=tmpData;
    end
    fclose(fid);
    
    %Output
%     LFPdata=double(LFPdata);


end