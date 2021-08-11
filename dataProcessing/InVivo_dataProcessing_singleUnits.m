function s=InVivo_dataProcessing_singleUnits(s,ops,stim)

    if isempty(s.suid)
        return
    end
    
    fprintf(strcat('Processing single units...','\n'))

    %% Obtain simple information about single units
    s.suage=ones(numel(s.suid),1)*ops.age;
    s=InVivo_dataProcessing_singleUnits_findLayerDepth(s,ops);        
    s=InVivo_dataProcessing_singleUnits_getTemplateFeatures(s,ops,1,ops.dirOUT);

    %% Run PSTH analysis
    s = InVivo_dataProcessing_singleUnits_getPSTHbyCondition (s,stim,ops,ops.verbose);

    
end
