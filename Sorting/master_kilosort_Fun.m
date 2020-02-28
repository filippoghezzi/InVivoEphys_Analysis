function master_kilosort_Fun(ID,dataFolders,sortingFolder) 

    addpath(genpath('C:\Users\Butt Lab\Documents\GitHub\Kilosort2')) % path to kilosort folder
    addpath('C:\Users\Butt Lab\Documents\GitHub\npy-matlab')

    %% Convert OpenEphys data to binary data
    ops.Nchan    = 32; % total number of channels in your recording
    ops.dataRoot=dataFolders;
    ops.binaryRoot=sortingFolder;
    ops.fbinary             = 'RawData'; % will be created for 'openEphys'	

    ops=convertOpenEphysToRawBInary_MergeRecs(ops);

    %% Set up sorting
    pathToYourConfigFile = 'C:\Users\Butt Lab\Documents\SpikeSorting'; % take from Github folder and put it somewhere else (together with the master_file)
    configFilename = strcat('StandardConfig_',ID);
    run(fullfile(pathToYourConfigFile,configFilename))
    rootH = sortingFolder;
    ops.fproc       = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
    % ops.chanMap = fullfile(pathToYourConfigFile, 'neuropixPhase3A_kilosortChanMap.mat');

    ops.trange = [0 Inf]; % time range to sort
    ops.NchanTOT    = 32; % total number of channels in your recording

    % the binary file is in this folder
    rootZ = sortingFolder;

    %% Sorting algorithm
    fprintf('Looking for data inside %s \n', rootZ)

    % is there a channel map file in this folder?
    fs = dir(fullfile(rootZ, 'chan*.mat'));
    if ~isempty(fs)
        ops.chanMap = fullfile(rootZ, fs(1).name);
    end

    % find the binary file
    fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
    ops.fbinary = fullfile(rootZ, fs(1).name);

    % preprocess data to create temp_wh.dat
    rez = preprocessDataSub(ops);

    % time-reordering as a function of drift
    rez = clusterSingleBatches(rez);
    save(fullfile(rootZ, 'rez.mat'), 'rez', '-v7.3');

    % main tracking and template matching algorithm
    rez = learnAndSolve8b(rez);

    % final merges
    rez = find_merges(rez, 1);

    % final splits by SVD
    rez = splitAllClusters(rez, 1);

    % final splits by amplitudes
    rez = splitAllClusters(rez, 0);

    % decide on cutoff
    rez = set_cutoff(rez);

    fprintf('found %d good units \n', sum(rez.good>0))

    % write to Phy
    fprintf('Saving results to Phy  \n')
    rezToPhy(rez, rootZ);

    %% Save Matlab file and eliminate temporary files
    % discard features in final rez file (too slow to save)
    rez.cProj = [];
    rez.cProjPC = [];

    % save final results as rez2
    fprintf('Saving final results in rez2  \n')
    fname = fullfile(rootZ, 'rez2.mat');
    save(fname, 'rez', '-v7.3');
    
    delete(ops.fproc); %Check HERE to eliminate temp file.

end
