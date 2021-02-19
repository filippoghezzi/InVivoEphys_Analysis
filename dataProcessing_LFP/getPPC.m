function [PPC,vectorAngle,vectorLength,pValue]=getPPC(spiketimes,filtLFP,varargin)
% Measure values of spike-LFP phase coupling as classical vector length and
% paiwise phase consistency. 
% Requires LFP signal given as input: this must be bandpass filtered in the 
% frequency where coupling needs to be measured.
% Vector length and angle are calculated using functions from P. Berens, CircStat: A Matlab Toolbox for Circular Statistics, Journal of Statistical Software, Volume 31, Issue 10, 2009.
% Pairwise phase consistency is measured according to Vinck et al., The pairwise phase consistency: A bias-free measure of rhythmic neuronal synchronization, NeuroImage (2010)
% Input: spiketimes -> spike times of a single unit in samples. 
%       filtLFP -> filtered single channel LFP signal
%
% Outputs: PPC-> numeric, paiwise phase consistency
%          vecotrAngle -> numeric, Mean direction of a sample of circular data
%          vectorLength -> numeric, Resultant vector length
%          pValue -> pValue of Rayleigh's test for nonuniformity
 
    addpath('C:\Users\Butt Lab\Documents\MATLAB\Add-Ons\CircStat2012a')
    
    %% Input 
    if nargin==3
        plotting=varargin{1};
    else
        plotting=0;
    end
    
    %% Pre-process filtered LFP
    hilbertLFP=hilbert(filtLFP);
    phaseLFP=angle(hilbertLFP);
    
    %% Get PPC
    theta=phaseLFP(spiketimes);
    vectorLength=circ_r(theta');
    vectorAngle=circ_mean(theta');
    pValue=circ_rtest(theta');
    PPC=ppc(theta');
    
    %% Plotting
    if plotting
        circ_plot(theta','hist',[],20,true,true,'linewidth',2,'color','r');
    end
end    

function PPC=ppc(theta)
    N=numel(theta);
    PPC=(2/(N*(N-1)))*sum(pdist(theta,@(x,y) cos(x)*cos(y)+sin(x)*sin(y)));
end