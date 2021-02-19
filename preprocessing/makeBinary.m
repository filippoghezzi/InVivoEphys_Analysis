function ops = makeBinary(ops,dirIN,dirOUT)
% function ops = makeBinary(ops,dirIN,dirOUT)
% Wrapper for OpenEphys2Binary.  RawData.dat will be saved after reordering channels by
% depth, according to ElectrodeMap (if reorderChannels == 1). 
% Input:   ops -> struct, needs to contain fbinary and ElectrodeMap fields.
%          dirIN ->  cell array (num files, 1), containing directories of
%                    all files in the recording. 
%          dirOUT -> output directory where RawData.dat is saved.
% Output:  ops ->  containing nSamplesBlocks for individual recordings.

    %% Set-up ops
    ops.Nchan    = numel(ops.ElectrodeMap); % total number of channels in your recording
    ops.dataRoot=dirIN;
    ops.binaryRoot=dirOUT;
    ops.fbinary             = 'RawData'; 
    ops.reorderChannels = 1;
    
    %% Convert OpenEphys data to binary data
    ops = openEphys2Binary(ops);
    
end
