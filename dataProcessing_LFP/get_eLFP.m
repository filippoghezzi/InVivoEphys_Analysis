function data=get_eLFP(filename,stimuli,fs,chanN,window,varargin)
% function get_eLFP(ops,filename,stimuli)

% Obtain average LFP aligned to stimuli; exploit KS2 RawData.dat
% Inputs:
%   ops -> struct, from KS2;
%   filename -> full directory for filename;
%   stimuli -> column vector containing stimuli onset (in sample).
    
    %Set filter parameters
    if isempty(varargin)
        filter=false;
    elseif size(varargin,1)==1
        filter=true;
        if strcmp(varargin{1},'LFP')
            [b,a]=butter(3,150/(fs/2),'low');
        elseif strcmp(varargin{1},'MUA')
            [b,a]=butter(3,[400 4000]/(fs/2),'bandpass');
        end
    end



    fid=fopen(filename,'r');
    frewind(fid);
    windowDuration=window(2)+window(1);
    data=zeros(chanN,round(windowDuration*fs),size(stimuli,1));

    for trial=1:size(stimuli,1)
        % Align file position to stimulus; remember binary file is in the
        % format Ch1,T1,Ch2,T1,...ChN,T1,Ch1,T2,...ChN,TN. First index is
        % 0, align onto (T-1)*2*Nchan with T stimulus onset sample. 
        
        T=stimuli(trial)-(window(1)*fs);
        offset = (2*chanN)*(T-1);
        
        fseek(fid, offset, 'bof');
        tmpData = fread(fid, [chanN round(windowDuration*fs)], '*int16');
        tmpData=double(tmpData);
        % Filter temp data
        if filter
            tmpData=tmpData';
            filtData=filtfilt(b,a,tmpData);            
            tmpData=filtData';
        end
        
        data(:,:,trial)=tmpData;
    
    end
    fclose(fid);
    
end