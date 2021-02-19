%%%%% Analysis to add: 
%                      - single units:  responsive, percentage of responsive neurons across development
%                      - check if spike template amplitude could be use to
%                      take away laser evoked noise from total spikes



clear 
close all
clc

tab=readtable('V1_In Vivo_SST;Ai32.csv');

spikeFolder='C:\Users\Butt Lab\Documents\SpikeSorting';
load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\A1x32_Map.mat');
tab=tab(tab.Use~=0,:);

%%%temporary
% tab=tab(tab.Temp==1,:);


%%%
%% Load sample duration
tab.SampleLength(1,1)=0;
tab.startSample(1,1)=0;
tab.endSample(1,1)=0;
for i=1:max(unique(tab.RecID))
% for i=1:8
    if ~isempty(tab(tab.RecID==i,:))
        load(fullfile(spikeFolder,tab(tab.RecID==i,:).MouseID{1},'SampleDuration.mat'))
        tab(tab.RecID==i,:).SampleLength=samplesToSave;
        for j=1:height(tab(tab.RecID==i,:))
            if j==1
                tab(tab.RecID==i,:).startSample(j)=1;
                tab(tab.RecID==i,:).endSample(j)=tab(tab.RecID==i,:).SampleLength(j);
            else
                tab(tab.RecID==i,:).startSample(j)=tab(tab.RecID==i,:).endSample(j-1)+1;
                tab(tab.RecID==i,:).endSample(j)=tab(tab.RecID==i,:).SampleLength(j)+tab(tab.RecID==i,:).startSample(j);
            end
        end
    end
end
tab = removevars(tab,{'SampleLength','endSample'});

AllSingleUnits = table;
% for i=8:max(unique(tab.RecID))
for i=6
    if ~isempty(tab(tab.RecID==i,:))
        [SingleUnits ,CSDinfo ]= Analysis_SingleAnimal(tab(tab.RecID==i,:),ElectrodeMap,spikeFolder);
    end
    save(fullfile(tab(tab.RecID==i,:).Folder{1},'CSDinfo.mat'),'CSDinfo')
end

% save('SU2.mat','AllSingleUnits')

%% Main data analysis
function [SingleUnits, CSDinfo] = Analysis_SingleAnimal(tab,ElectrodeMap,spikeFolder)
    savingFolder=tab.Folder{1};
    
    optotaggingSession=1;
    visualSession=1;
    baselineSession=1;
    tab=sortrows(tab,'Use');
    [spikes,templates,suid]=LoadSpikes(fullfile(spikeFolder,tab.MouseID{1}),ElectrodeMap);

    for i=1:height(tab)
        dataFolder=fullfile(tab.Folder{i},tab.Experiment{i});
        dataFolderContent=dir(dataFolder);

        %% Load LFP data
        dataPresence=[];
        for j=1:length(dataFolderContent)
            dataPresence=[dataPresence;strcmp('LFP_Data.mat',dataFolderContent(j).name)];
        end
        if ~any(dataPresence)
            [LFP]=loadLFPData(dataFolder,ElectrodeMap);
            save(fullfile(dataFolder,'LFP_Data.mat'),'LFP')
        else
            load(fullfile(dataFolder,'LFP_Data.mat'))
        end

        %% Load spike data
        eventIdx=loadEventsForSpikes(dataFolder);

        %% Choose analysis
        switch tab.Protocol{i}
            case 'VisualFlash'
                suModV(size(suid,1),visualSession)=0;
                CSDinfo=eLFP_Analysis(LFP.data,LFP.dataTime,LFP.eventArray,0,[],'LED', spikes(ismember(spikes(:,2),[1,2]),:), eventIdx,tab.startSample(i));
                figname=strcat(tab.MouseID{i},'- P',int2str(tab.Age(i)),' - Visual Flash LFP - Laser OFF - 1');
                sgtitle(figname)
                [suModV, PSTH] = Spike_Analysis(spikes,templates,suid,eventIdx,CSDinfo.L4,tab.startSample(i),0,savingFolder,figname,suModV,visualSession);
                figname=strcat(tab.MouseID{i},'- P',int2str(tab.Age(i)),' - Visual Flash Units - Laser OFF - 1');
                sgtitle(figname)
                export_fig(fullfile(savingFolder,figname),'-tiff','-transparent')
                close
                visualSession=visualSession+1;
            case 'VisualFlash_Opto'
                suModV(size(suid,1),visualSession)=0;
                suModI(size(suid,1),baselineSession)=0;
                CSDinfo=eLFP_Analysis(LFP.data,LFP.dataTime,LFP.eventArray,0,[],'LED', spikes(ismember(spikes(:,2),[1,2]),:),eventIdx,tab.startSample(i));
                figname=strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - Visual Flash LFP - Laser OFF');
                sgtitle(figname)
                [suModV, PSTH] = Spike_Analysis(spikes,templates,suid,eventIdx,CSDinfo.L4,tab.startSample(i),0,savingFolder,figname,suModV,visualSession);
                figname=strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - Visual Flash Units - Laser OFF');
                sgtitle(figname)
                export_fig(fullfile(savingFolder,figname),'-tiff','-transparent')
                close
                eLFP_Analysis(LFP.data,LFP.dataTime,LFP.eventArray,1,CSDinfo,'LED', spikes(ismember(spikes(:,2),[1,2]),:),eventIdx,tab.startSample(i));
                figname=strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - Visual Flash LFP - Laser ON');
                sgtitle(figname)
                [suModI, ~] = Spike_Analysis(spikes,templates,suid,eventIdx,CSDinfo.L4,tab.startSample(i),1,savingFolder,figname,suModI,baselineSession);
                figname=strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - Visual Flash Units - Laser ON');
                sgtitle(figname)
                export_fig(fullfile(savingFolder,figname),'-tiff','-transparent')
                close
            case 'Baseline'
%                 Baseline_Analysis(LFP.data,LFP.dataTime,LFP.eventArray)
            case 'Baseline_Opto'
                suModI(size(suid,1),baselineSession)=0;
                eLFP_Analysis(LFP.data,LFP.dataTime,LFP.eventArray,1,CSDinfo,'Laser', spikes(ismember(spikes(:,2),[1,2]),:),eventIdx,tab.startSample(i));
                figname=strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - LFP Laser only');
                sgtitle(figname)
                [suModI, ~] = Spike_Analysis(spikes,templates,suid,eventIdx,CSDinfo.L4,tab.startSample(i),1,savingFolder,figname,suModI,baselineSession);
                figname=strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - Units Laser only');
                sgtitle(figname)
                export_fig(fullfile(savingFolder,figname),'-tiff','-transparent')
                close
                baselineSession =baselineSession +1;
            case 'Optotagging'
                suModO(size(suid,1),optotaggingSession)=0;
                suModO=Optotagging_Analysis(spikes,templates,suid,eventIdx,tab.startSample(i),optotaggingSession,suModO,CSDinfo.L4);
                figname=strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - Optotagging',int2str(optotaggingSession));
                sgtitle(figname)
                export_fig(fullfile(savingFolder,figname),'-tiff','-transparent')
                close
                optotaggingSession=optotaggingSession+1;
            case 'ComplexSession'
                
        end

    end
    
    suModO = sum(suModO,2);
    suModV = sum(suModV,2);
    suModI = sum(suModI,2);
    
    SingleUnits = table;
    SingleUnits.ID = suid;
    SingleUnits.Optotagging = suModO;
    SingleUnits.VisualResponsive = suModV;
    SingleUnits.Inhibited = suModI;
    SingleUnits.STH = PSTH;
    SingleUnits.Mouse(:) = tab.MouseID(1);    
    
end

%% Analysis functions
function [suModV, PSTH] = Spike_Analysis(spikes,templates,suid,eventIdx,L4,start,opto,savingFolder,figname,suModV,session)
%% Plot multi unit activity histogram
    mua=spikes(ismember(spikes(:,2),[1,2]),:); % Select both single units and multi units
    mua=mua(ismember(mua(:,3),L4),:); % Only in Layer 4
    
    %% Find events with laser ON or OFF
    if ~isempty(eventIdx{1})
        if ~isempty(eventIdx{3}) %Case both laser and LED
            for i=1:size(eventIdx{3},1)
                LaserON(i,1)=any(eventIdx{3}(i,1)-eventIdx{1}(:,1)==1000);
            end
            if opto
                stimulus=eventIdx{3}(LaserON,:)+start;
            else
                stimulus=eventIdx{3}(~LaserON,:)+start;
            end
        else %Case laser only
            stimulus=eventIdx{1}+start; %3 for LED
        end
    else %Case LED only
        stimulus=eventIdx{3}+start; %3 for LED
    end
    
    muaTimes=[];
    for i=1:size(stimulus,1)
        tmpSpikes(:,4)=mua(:,4)-stimulus(i,1);
        muaTimes=[muaTimes;tmpSpikes(tmpSpikes(:,4)>-20000&tmpSpikes(:,4)<5*20000,4)];
    end
    muaTimes=double(muaTimes)./20; % Convert to double and divide by 20 to get time in ms
 
    subplot(4,1,3)
    histogram(muaTimes,600,'EdgeColor','none','FaceColor','k')
    ax=gca;
    ax.XLim=[-1000, 2000];
    ax.XTick=[];
    ax.YLabel.String='Spike #';
    ax.YLabel.FontSize=15;
    ax.XAxis.Color='none'; 
    box off
    export_fig(fullfile(savingFolder,figname),'-tiff','-transparent')
    close

%% Plot single units
    figure('units','normalized','outerposition',[0 0 1 1]);
    for i=1:size(suid,1)
        su=spikes(ismember(spikes(:,1),suid(i)),:);
        if ismember(su(1,3),L4)
            Layer='L4';
        elseif su(1,3)<min(L4)
            Layer='Ingrafranular';         
        elseif su(1,3)>max(L4)
            Layer='Supragranular';         
        end
        subplot(7,8,i)
        tempPSTH = [];
        for j=1:size(stimulus,1)
            tmpSu=su(:,4)-stimulus(j);
            raster(i,j)={tmpSu(tmpSu>-20000&tmpSu<5*20000)};
            raster{i,j}=double(raster{i,j})./20;            
            scatter(raster{i,j},ones(size(raster{i,j},1),1)*j,'.k')
            hold on 
            tempPSTH = [tempPSTH;raster{i,j}]; 
        end
        PSTH(i,1)={tempPSTH};
        if ~isempty(eventIdx{3}) 
            patch([0 100 100 0],[0 0 50 50],'y','FaceAlpha',.3,'EdgeColor','none')
        end
        if opto
            patch([-50 150 150 -50],[0 0 50 50],'c','FaceAlpha',.1,'EdgeColor','none')
        end
        hold off
        title(strcat('N',int2str(suid(i)),' - ',Layer),'FontSize',12)
        ax=gca;
        ax.XLim=[-1000, 2000];
        ax.YLim=[0, 50];
        ax.YTick=[];
        ax.YAxis.Color='none'; 
        box off
        
        %Find responsive units
        if ~opto
            pre=@(x) nnz(x>-100&x<0);
            post=@(x) nnz(x>1&x<100);
        
            isResponsive(:,1)=cellfun(pre,raster(i,:))';
            isResponsive(:,2)=cellfun(post,raster(i,:))';
        
            [~,p]=ttest(isResponsive(:,1),isResponsive(:,2));
            if p<=0.05 && sum(isResponsive(:,2))>size(isResponsive,1)*.3       
                if sum(isResponsive(:,2))>sum(isResponsive(:,1))
                    suModV(i,session)=1;
                else
%                     suModV(i,session)=2;
                end    
            end
            title(strcat('N',int2str(suid(i,1)),' -',Layer,' - Group:',int2str(suModV(i,session))),'FontSize',12)
        
        elseif opto
            pre=@(x) nnz(x>-250&x<-50);
            post=@(x) nnz(x>-50&x<150);
        
            isInhibited(:,1)=cellfun(pre,raster(i,:))';
            isInhibited(:,2)=cellfun(post,raster(i,:))';
        
            [~,p]=ttest(isInhibited(:,1),isInhibited(:,2));
            if p<=0.05 && sum(isInhibited(:,1))>size(isInhibited,1)*.3       
                if sum(isInhibited(:,1))>sum(isInhibited(:,2))
                    suModV(i,session)=1;
                else
%                     suModV(i,session)=2;
                end    
            end
            title(strcat('N',int2str(suid(i,1)),' -',Layer,' - Group:',int2str(suModV(i,session))),'FontSize',12) 
            
        end
        
    end
end

function CSDinfo=eLFP_Analysis(data,dataTime,eventArray,opto,CSDinfo,alignTrials, mua,eventIdx,start)
    %Event for spikes, used for MUA
    %% Find events with laser ON or OFF
    if ~isempty(eventIdx{1})
        if ~isempty(eventIdx{3}) %Case both laser and LED
            for i=1:size(eventIdx{3},1)
                LaserON(i,1)=any(eventIdx{3}(i,1)-eventIdx{1}(:,1)==1000);
            end
            if opto
                stimulus=eventIdx{3}(LaserON,:)+start;
            else
                stimulus=eventIdx{3}(~LaserON,:)+start;
            end
        else %Case laser only
            stimulus=eventIdx{1}+start; %3 for LED
        end
    else %Case LED only
        stimulus=eventIdx{3}+start; %3 for LED
    end
    

    %% Current source density
    LED=3;          %3, index of eventArray line for LED visual stimulation;
    laser=1;        %1, index of eventArray line for 470nm laser;
    
    if strcmp(alignTrials,'LED')
        % Align on LED onset
        LEDOnsetIdx=find([0,diff(eventArray(LED,:))]>0);
        LEDOffsetIdx=find([0,diff(eventArray(LED,:))]<0)-1;
        
        %Select trial in which laser is OFF or ON
        LEDOnsetIdx=LEDOnsetIdx(eventArray(laser,LEDOnsetIdx)==opto);
        LEDOffsetIdx=LEDOffsetIdx(eventArray(laser,LEDOffsetIdx)==opto);
        stimuliIdx=[LEDOnsetIdx',LEDOffsetIdx'];
        
    elseif strcmp(alignTrials,'Laser')
        % Align on Laser
        LaserOnsetIdx=find([0,diff(eventArray(laser,:))]>0); 
        LaserOffsetIdx=find([0,diff(eventArray(laser,:))]<0)-1;
        stimuliIdx=[LaserOnsetIdx',LaserOffsetIdx'];
        
    end

    if opto
        [CSD,~]=getAverageCSD(data,[stimuliIdx(:,1) stimuliIdx(:,2)+5100],1000,800,25,1,1,0, mua,stimulus);
    else
        [CSD,CSDinfo]=getAverageCSD(data,[stimuliIdx(:,1) stimuliIdx(:,2)+5100],1000,800,25,1,1,1, mua,stimulus);
    end

    %% Evoked LFP
    factor = 1; % the time factor (20 or 1 per msecond depending on sampling (MAKE IT AUTOMATIC ACCORDING TO DOWNSAMPLING)

    evokedLFP = [];
    i=0;
    for channel=1:size(data,1)
        if any(channel==CSDinfo.L4)
            i=i+1;
            for trial=1:size(stimuliIdx,1)
                LFPwindow=(stimuliIdx(trial,1)-1000*factor:stimuliIdx(trial,1)+100+5000);    
                evokedLFP(i,trial,1:length(LFPwindow))=data(channel,LFPwindow);
            end
        end
    end 

    %% Time-frequency plot and all plots
    figure('units','normalized','outerposition',[0 0 1 1]);
    meanLFP=mean(squeeze(mean(evokedLFP(:,:,:),1)),1);
    meanSweep=squeeze(mean(evokedLFP,1));
    subplot(4,1,1)
    hold on
    for sweep=1:size(evokedLFP,2)
        PlotAllData(meanSweep(sweep,:),1000,'Grey')
    end
    PlotAllData(meanLFP,1000,'k');
    ax=gca;
    ax.XLim=[-1000, 2000];
    ax.YLim=[min(meanLFP,[],'all')-20 max(meanLFP,[],'all')+20];
    ax.YLabel.String='Voltage (\muV)';
    ax.YLabel.FontSize=15;
    ax.XAxis.Color='none';
    
    if strcmp(alignTrials,'LED')
        patch([0 100 100 0],[min(meanLFP)-20 min(meanLFP)-20 max(meanLFP)+20 max(meanLFP)+20],'y','FaceAlpha',.2,'EdgeColor','none')
    end
    if opto
        if strcmp(alignTrials,'Laser')
            patch([0 200 200 0],[min(meanLFP)-20 min(meanLFP)-20 max(meanLFP)+20 max(meanLFP)+20],'c','FaceAlpha',.1,'EdgeColor','none')
        else
            patch([-50 150 150 -50],[min(meanLFP)-20 min(meanLFP)-20 max(meanLFP)+20 max(meanLFP)+20],'c','FaceAlpha',.1,'EdgeColor','none')
        end
    end
    box off
    
    
    subplot(4,1,4)
    imagesc([-1000 5000],[0 750],CSD,[min(CSD(:))/5,max(CSD(:))/5])
    colormap('jet')
    ax=gca;
    ax.XLim=[-1000, 2000];
    ax.YLabel.String='Cortical depth (\mum)';
    ax.YLabel.FontSize=15;
    ax.XLabel.String='Time (ms)';
    ax.XLabel.FontSize=15;
    box off
   
    subplot(4,1,2)
    timeFrequencyAnalysis(squeeze(mean(evokedLFP,1)))
    colorbar('off')
    ax=gca;
    ax.XLim=[-1000, 2000];
    ax.XTick=[];
    ax.YLabel.String='Frequency (Hz)';
    ax.YLabel.FontSize=15;
    ax.XAxis.Color='none'; 
    box off


end

function Baseline_Analysis(data,dataTime,eventArray)

end

function suModO=Optotagging_Analysis(spikes,templates,suid,eventIdx,start,session,suModO,L4)
    %% this function create the variable suMod with columns (one per session) with optotagging modulation in which 0=no modulation; 1=excitation (e.g.tagged units); 2=inhibition.
    laser=eventIdx{1}+start;

    figure('units','normalized','outerposition',[0 0 1 1]);
    for i=1:size(suid,1)
        su=spikes(ismember(spikes(:,1),suid(i)),:);
        subplot(7,7,i)
        if ismember(su(1,3),L4)
            Layer='L4';
        elseif su(1,3)<min(L4)
            Layer='Ingrafranular';         
        elseif su(1,3)>max(L4)
            Layer='Supragranular';         
        end
        for j=1:size(laser,1)
            tmpSu=su(:,4)-laser(j);
            raster(i,j)={tmpSu(tmpSu>-10000&tmpSu<10000)};
            raster{i,j}=double(raster{i,j})./20;
            scatter(raster{i,j},ones(size(raster{i,j},1),1)*j,'.k')
            hold on      
        end
        patch([0 50 50 0],[0 0 100 100],'c','FaceAlpha',.3,'EdgeColor','none')
        patch([-300 -250 -250 -300],[0 0 100 100],'r','FaceAlpha',.3,'EdgeColor','none')
        patch([250 300 300 250],[0 0 100 100],'r','FaceAlpha',.3,'EdgeColor','none')

        hold off
        ax=gca;
        ax.XLim=[-500, 500];
        ax.YLim=[0, 50];
        ax.YTick=[];
        ax.YAxis.Color='none'; 
        box off
        
        %% Find optotagged neurons
        pre=@(x) nnz(x>-300&x<-250);
        opto=@(x) nnz(x>1&x<50);
        post=@(x) nnz(x>250&x<300);
        
        isTagged(:,1)=cellfun(pre,raster(i,:))';
        isTagged(:,2)=cellfun(opto,raster(i,:))';
        isTagged(:,3)=cellfun(post,raster(i,:))';

        [p,~,stats]=anova1(isTagged,[],'off');
        if p<=0.05      
           c=multcompare(stats,'Display','off');
           if c(1,6)<0.05 && c(3,6)<0.05
               if sum(isTagged(:,2))>sum(isTagged(:,1)) && sum(isTagged(:,2))>sum(isTagged(:,3)) && sum(isTagged(:,2))>size(isTagged,1)*.2    %Higher than pre and post and spiking 30% of the trials
                   suModO(i,session)=1;
               end
%            elseif c(1,6)<0.05 && sum(isTagged(:,1))>size(isTagged,1)*.3
%                suModO(i,session)=2;
           end
        end
        title(strcat('N',int2str(suid(i,1)),' -',Layer,' - Group:',int2str(suModO(i,session))),'FontSize',12)

        
    end
    
end