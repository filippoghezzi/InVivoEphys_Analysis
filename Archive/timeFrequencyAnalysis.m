function [CFS,f]=timeFrequencyAnalysis(LFP,fs)
%function timeFrequencyAnalysis(evokedLFP,t,fs)
%
% Build time-frequency matrix (CFS) using continuous wavelet decompositon. Anaysis is
% performed only on the first channel of the LFP in input independently for 
% every trial and then averaged. 
% Input:
%       evokedLFP -> [Nchan Timestamps Trial];
%       fs -> sampling frequency.
% Output:
%       CSF -> time-frequency matrix [frequency timestamps];
%       f -> frequency vector generated from cwt.
        
    LFP=LFP';

    for trial=1:size(LFP,1)
        [tmp_cfs,f] = cwt(LFP(trial,:),fs);
        tmp_cfs=abs(tmp_cfs);

        if trial==1
            CFS=tmp_cfs;
        else
            CFS=CFS+tmp_cfs;
        end
    end

    CFS=CFS./trial;
end
