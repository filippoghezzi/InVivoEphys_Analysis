function ops = openEphys2Binary(ops)
% function ops = openEphys2Binary(ops)
%
% Modified from KS2 convertOpenEphysToRawBInary to merge different
% recordings belonging to the same experiment and reorder according to
% electrode channel map.
 
fname       = fullfile(ops.binaryRoot, sprintf('%s.dat', ops.fbinary)); 

for folder=1:size(ops.dataRoot,1)
    if folder==1
        fidout=fopen(fname,'w');
    else
        fidout=fopen(fname,'A');
    end
    %
    clear fs
    for j = 1:ops.Nchan
        %%%% This condition account for problem in pausing then restarting the recording.%%%%
        if exist(fullfile(ops.dataRoot{folder},'100_CH17_2.continuous'),'file')   
           fs{j} = dir(fullfile(ops.dataRoot{folder}, sprintf('100_CH%d_2.continuous', j+16) )); 
        elseif exist(fullfile(ops.dataRoot{folder},'100_CH17.continuous'),'file') 
           fs{j} = dir(fullfile(ops.dataRoot{folder}, sprintf('100_CH%d.continuous', j+16) )); % Added +16 because using ch 17-48  
        %%%% This condition account for recordings done on Liad's rig.%%%%
        elseif exist(fullfile(ops.dataRoot{folder},'107_CH17.continuous'),'file') 
            fs{j} = dir(fullfile(ops.dataRoot{folder}, sprintf('107_CH%d.continuous', j+16) )); % Added +16 because using ch 17-48  
        end
    end
    nblocks = cellfun(@(x) numel(x), fs);
    if numel(unique(nblocks))>1
       error('different number of blocks for different channels!') 
    end
    %
    nBlocks     = unique(nblocks);
    nSamples    = 1024;  % fixed to 1024 for now!

    fid = cell(ops.Nchan, 1);

    for k = 1:nBlocks
        for j = 1:ops.Nchan
            fid{j}             = fopen(fullfile(ops.dataRoot{folder}, fs{j}(k).name));
            % discard header information
            fseek(fid{j}, 1024, 0);
        end
        %
        nsamps = 0;
        flag = 1;
        while 1
            samples = zeros(nSamples * 1000, ops.Nchan, 'int16');
            for j = 1:ops.Nchan
                collectSamps    = zeros(nSamples * 1000, 1, 'int16');

                rawData         = fread(fid{j}, 1000 * (nSamples + 6), '1030*int16', 10, 'b');

                nbatches        = ceil(numel(rawData)/(nSamples+6));
                for s = 1:nbatches
                    if s==nbatches
                        continue
                    end
                    rawSamps = rawData((s-1) * (nSamples + 6) +6+ [1:nSamples]);
                    collectSamps((s-1)*nSamples + [1:nSamples]) = rawSamps;
                end
                samples(:,j)         = collectSamps;
            end

            if nbatches<1000
                flag = 0;
            end
            if flag==0
                samples = samples(1:s*nSamples, :);
            end

            samples         = samples';

            if ops.reorderChannels
                samples=samples(ops.ElectrodeMap-16,:);
            end
            
            fwrite(fidout, samples, 'int16');

            nsamps = nsamps + size(samples,2);

            if flag==0
                break;
            end
        end
        ops.nSamplesBlocks(folder,k) = nsamps;

        for j = 1:ops.Nchan
           fclose(fid{j}); 
        end

    end

    fclose(fidout);
end