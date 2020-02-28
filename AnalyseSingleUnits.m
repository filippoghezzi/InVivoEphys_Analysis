clear 
close all
clc

load ('SU.mat')

load('InVivo_V1_SST;Ai32.mat')

for i = 1:height(AllSingleUnits)
    AllSingleUnits.Age(i)=tab(strcmp(tab.MouseID,AllSingleUnits.Mouse(i)),:).Age(1);
    
end


DevWindow=[5,8;9,13;14,18];
figure
for i = 1:3
    
    Dev = AllSingleUnits((AllSingleUnits.Age>=DevWindow(i,1)) & (AllSingleUnits.Age<=DevWindow(i,2)),:);
    Pyr_No=0;
    Pyr_Yes=0;
    SST_No=0;
    SST_Yes=0;
    for j = 1:height(Dev)
        if Dev.Optotagging(j)==0
            if Dev.VisualResponsive(j)==0 
                Pyr_No = Pyr_No+1;
            else
                Pyr_Yes = Pyr_Yes +1;
            end
        else 
            if Dev.VisualResponsive(j)==0
                SST_No = SST_No+1;
            else 
                SST_Yes = SST_Yes +1;
            end
        end
        
    end
    
    subplot(2,3,i)
    pie([Pyr_No,Pyr_Yes],{'Non Responsive','Responsive'})
    subplot(2,3,i+3)
    pie([SST_No,SST_Yes],{'Non Responsive','Responsive'})
    
    PercentageData(i,1:4)=[Pyr_No,Pyr_Yes,SST_No,SST_Yes];
end

PercentageData



%%% TO SEE SST inhibited units compared to other units


SSTinhibition=AllSingleUnits((AllSingleUnits.Age>=DevWindow(3,1)) & (AllSingleUnits.Age<=DevWindow(3,2)),:);
SSTinhibition=SSTinhibition(SSTinhibition.VisualResponsive~=0,:);
SSTinhibition=SSTinhibition(SSTinhibition.Optotagging==0,:);
j=1;
k=1;

for i = 1: height(SSTinhibition)
        SSTinhibition.PSTH(i)={histcounts(SSTinhibition.STH{i},200).*1000./25./50};
        if SSTinhibition.Inhibited(i)==0
            NonInhibited(j,:)=SSTinhibition.PSTH{i}./mean(SSTinhibition.PSTH{i}(1:20));
            j=j+1;
        else 
            Inhibited(k,:)=SSTinhibition.PSTH{i}./mean(SSTinhibition.PSTH{i}(1:20));
            k=k+1;
        end
end



figure
hold on
plot(mean(Inhibited,1),'k')
plot(mean(NonInhibited,1),'r')
