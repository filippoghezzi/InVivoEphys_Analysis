function LFP=loadLFP_baseline(binaryFilename,fs,chanN,startBaseline,endBaseline,varargin)
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
%   varargin -> empty for LFP without filtering, 'LFP' for [0 100] Hz bandpass signal; 'MUA' for [400 4000] Hz bandpass signal.
    
    fprintf(strcat('Loading LFP ...','\n'))

    %Set filter parameters
    if isempty(varargin)
        doFiltering=false;
    elseif size(varargin,1)==1
        doFiltering=true;
        if strcmp(varargin{1},'LFP')
            [b,a]=butter(3,100/(fs/2),'low');
        elseif strcmp(varargin{1},'MUA')
            [b,a]=butter(3,[400 4000]/(fs/2),'bandpass');
        end
    end

    fid=fopen(binaryFilename,'r');
    frewind(fid);
    baselineDuration=endBaseline-startBaseline; %samples
    LFP=zeros(chanN,round(baselineDuration));

    batchSize=10*fs;
    Nbatches=baselineDuration/batchSize;
    
    for batch=1:Nbatches-1
        % Align file position to stimulus; remember binary file is in the
        % format Ch1,T1,Ch2,T1,...ChN,T1,Ch1,T2,...ChN,TN. First index is
        % 0, align onto (T-1)*2*Nchan with T stimulus onset sample. 
        
        T=((batch-1)*batchSize+1)+startBaseline;
        offset = (2*chanN)*(T-1);
        
        fseek(fid, offset, 'bof');
        tmpData = fread(fid, [chanN batchSize], '*int16');
        tmpData=double(tmpData);
        
        % Filter temp data
        if doFiltering
            filtData=tmpData';            
            filtData=filtfilt(b,a,filtData);            
            tmpData=filtData';
        end
        LFP(:,(batch-1)*batchSize+1:batch*batchSize)=tmpData;
    
    end
    fclose(fid);
    
end