function getFilteredBinary(ops)
% function getFilteredBinary(ops)
%
% Modified from KS2.preprocessDataSub(ops) to only get BP filtered binary
% file.

    ops.fBPbinary_temp=fullfile(ops.binaryRoot,'temp_BP.dat');

    %% Set up 
    NT       = ops.NT ;
    NchanTOT = ops.NchanTOT;
    Nbatch   = ops.Nbatch;
    NTbuff   = NT + 4*ops.ntbuff;
    chanMap  = ops.chanMap;

    fid         = fopen(ops.fbinary, 'r');
    fidW        = fopen(ops.fBPbinary_temp,   'w');
    DATA        = [];

    if isfield(ops,'fslow')&&ops.fslow<ops.fs/2
        [b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass');
    else
        [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high');
    end

    for ibatch = 1:Nbatch
        %% Offset
        offset = max(0, ops.twind + 2*NchanTOT*((NT - ops.ntbuff) * (ibatch-1) - 2*ops.ntbuff));
        if offset==0
            ioffset = 0;
        else
            ioffset = ops.ntbuff;
        end
        fseek(fid, offset, 'bof');

        %% Load data
        buff = fread(fid, [NchanTOT NTbuff], '*int16');
        if isempty(buff)
            break;
        end
        nsampcurr = size(buff,2);
        if nsampcurr<NTbuff
            buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr);
        end

        if ops.GPU
            dataRAW = gpuArray(buff);
        else
            dataRAW = buff;
        end
        dataRAW = dataRAW';
        dataRAW = single(dataRAW);
%         dataRAW = dataRAW(:, chanMap);

        %% Subtract the mean from each channel
        dataRAW = dataRAW - mean(dataRAW, 1);    

        %% Butterworth Filter
        datr = filter(b1, a1, dataRAW);
        datr = flipud(datr);
        datr = filter(b1, a1, datr);
        datr = flipud(datr);

        %% CAR, common average referencing by median 
        if getOr(ops, 'CAR', 1)
            datr = datr - median(datr, 2);
        end

        datr = datr(ioffset + (1:NT),:);

    %     datr    = datr * Wrot;   ????????????????????????????????????????????

        %% Write data
        datcpu  = gather_try(int16(datr));
        fwrite(fidW, datcpu, 'int16');

    end
    
fclose(fidW);
fclose(fid);