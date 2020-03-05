function doClustering(PSTH)



[coeff,score,latent,tsquared,explained]=pca(zscore(allPSTH));

for i=1:size(explained,1)
    if i==1
        VarExp(i)=explained(i);
    else
        VarExp(i)=VarExp(i-1)+explained(i);
    end
end
figure
plot(VarExp,'.-')


for j=1:100
for k=1:10
    figure
    
    idx=kmeans(score(:,1:40),k); %Kmeans only the first components that explain enough variance
%     T = clusterdata(score(:,1:10),cutoff);


    subplot(1,2,1)
    hold on
    scatter(score(idx==1,1),score(idx==1,2),'r')
    scatter(score(idx==2,1),score(idx==2,2),'b')
    scatter(score(idx==3,1),score(idx==3,2),'g')
    scatter(score(idx==4,1),score(idx==4,2),'y')
    scatter(score(idx==5,1),score(idx==5,2),'k')


%     subplot(1,2,2)
    hold on
    scatter(score(cellIdx==1,1),score(cellIdx==1,2),'k+')
    scatter(score(cellIdx==2,1),score(cellIdx==2,2),'k*')
    scatter(score(cellIdx==3,1),score(cellIdx==3,2),'kx')
    scatter(score(cellIdx==4,1),score(cellIdx==4,2),'k^')
    legend('C1','C2','C3','C4','L2/3','L4','L5','L6')
%     subplot(1,2,2)
    title(strcat('k=',int2str(k)))
    s = silhouette(score(:,1:10),idx);
%     silVal(j,k-1)=mean(s);
end
% end
figure
% boxplot(silVal)