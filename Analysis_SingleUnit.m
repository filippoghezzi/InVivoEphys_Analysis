clear 
close all
clc

tab=readtable('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\V1_InVivo_SST;Ai32.csv'); 

spikeFolder='C:\Users\Butt Lab\Documents\SpikeSorting';
load('C:\Users\Butt Lab\Documents\GitHub\InVivoEphys_Analysis\ElectrodeMaps\A1x32_Map.mat');
tab=tab(tab.Use~=0,:);

tab=loadSampleDuration(tab,spikeFolder);



for i=max(unique(tab.RecID))
    if ~isempty(tab(tab.RecID==i,:))
        [SingleUnits ,CSDinfo ]= Analysis_SingleAnimal(tab(tab.RecID==i,:),ElectrodeMap,spikeFolder);
    end
end


%% Main data analysis
function [SingleUnits, CSDinfo] = Analysis_SingleAnimal(tab,ElectrodeMap,spikeFolder)
    savingFolder=tab.Folder{1};
    
    optotaggingSession=1;
    visualSession=1;
    baselineSession=1;
    tab=sortrows(tab,'Use');
    [spikes,templates,suid]=LoadSpikes(fullfile(spikeFolder,tab.MouseID{1}),ElectrodeMap);
    load(fullfile(savingFolder,'CSDinfo.mat'))

    for i=1:height(tab)
        %% Load spike data
        dataFolder=fullfile(tab.Folder{i},tab.Experiment{i});
        eventIdx=loadEventsForSpikes(dataFolder);
        
        
        %% Choose analysis
        switch tab.Protocol{i}
            case 'VisualFlash'
                suModV(size(suid,1),visualSession)=0;
                [suModV, PSTH] = Spike_Analysis(spikes,templates,suid,eventIdx,CSDinfo.L4,tab.startSample(i),0,suModV,visualSession,30000);
                figname=strcat(tab.MouseID{i},'- P',int2str(tab.Age(i)),' - Visual Flash Units - Laser OFF - 1');
%                 sgtitle(figname)
%                 export_fig(fullfile(savingFolder,figname),'-tiff','-transparent')
                close
                visualSession=visualSession+1;
            case 'VisualFlash_Opto'
                suModV(size(suid,1),visualSession)=0;
                suModI(size(suid,1),baselineSession)=0;
                [suModV, PSTH] = Spike_Analysis(spikes,templates,suid,eventIdx,CSDinfo.L4,tab.startSample(i),0,suModV,visualSession,30000);
                figname=strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - Visual Flash Units - Laser OFF');
%                 sgtitle(figname)
%                 export_fig(fullfile(savingFolder,figname),'-tiff','-transparent')
                close
                [suModI, ~] = Spike_Analysis(spikes,templates,suid,eventIdx,CSDinfo.L4,tab.startSample(i),1,suModI,baselineSession,30000);
                figname=strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - Visual Flash Units - Laser ON');
%                 sgtitle(figname)
%                 export_fig(fullfile(savingFolder,figname),'-tiff','-transparent')
                close
            case 'Baseline'
%                 Baseline_Analysis(LFP.data,LFP.dataTime,LFP.eventArray)
            case 'Baseline_Opto'
                suModI(size(suid,1),baselineSession)=0;
                [suModI, ~] = Spike_Analysis(spikes,templates,suid,eventIdx,CSDinfo.L4,tab.startSample(i),1,suModI,baselineSession,30000);
                figname=strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - Units Laser only');
%                 sgtitle(figname)
%                 export_fig(fullfile(savingFolder,figname),'-tiff','-transparent')
                close
                baselineSession =baselineSession +1;
            case 'Optotagging'
                suModO(size(suid,1),optotaggingSession)=0;
                suModO=Optotagging_Analysis(spikes,templates,suid,eventIdx,tab.startSample(i),optotaggingSession,suModO,CSDinfo.L4);
                figname=strcat(tab.MouseID{i},' - P',int2str(tab.Age(i)),' - Optotagging',int2str(optotaggingSession));
%                 sgtitle(figname)
%                 export_fig(fullfile(savingFolder,figname),'-tiff','-transparent')
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
function [suModV, PSTH] = Spike_Analysis(spikes,templates,suid,eventIdx,L4,start,opto,suModV,session,sr)  

    eventIdx=eventIdx+start;
    
    [PSTH,PSTHbins,raster]=makePSTHandRaster(spikes,suid,eventIdx(:,2),sr,[1,5],0.03);
    plotRaster(raster,L4)
    plotPSTH(PSTH,PSTHbins,L4)
    findVisualResponsiveUnits(raster)
        
        %Find responsive units
        if ~opto

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