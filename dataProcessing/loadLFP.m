function LFP=loadLFP(binaryFilename,stimuli,fs,chanN,window,varargin)
% function loadLFP(filename,stimuli,fs,chanN,window)

% Obtain LFP traces aligned to stimuli according to window; exploit KS2
% RawData.dat. To load a continuous LFP trace (e.g. baseline), set
% stimuli=startBaseline (in samples) and window=[0,endBaseline]).
% Inputs:
%   binaryFilename -> full directory for binary filename;
%   stimuli -> column vector containing stimuli onset (in sample).
%   fs -> sampling frequency;
%   chanN -> total number of channels;
%   window -> [pre post] values (in seconds) related to stimuli values to collect LFP
%   varargin -> empty for LFP without filtering, 'LFP' for [0 150] Hz bandpass signal; 'MUA' for [400 4000] Hz bandpass signal.
    
    fprintf(strcat('Loading LFP ...','\n'))

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

    fid=fopen(binaryFilename,'r');
    frewind(fid);
    windowDuration=window(2)+window(1);
    LFP=zeros(chanN,round(windowDuration*fs),size(stimuli,1));

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
        
        LFP(:,:,trial)=tmpData;
    
    end
    fclose(fid);
    
end