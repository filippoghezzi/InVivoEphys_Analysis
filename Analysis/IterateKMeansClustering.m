function Klabel = IterateKMeansClustering(data, Kvalues, varargin)

%% Input parser
    p=inputParser;
    addRequired(p,'data',@(x) isnumeric(x));
    addRequired(p,'Kvalues',@(x) isnumeric(x));
    
    addOptional(p,'TrueLabels',[]);
    addOptional(p,'NumIterations',100,@(x) isnumeric(x));
    addOptional(p,'Distance','cosine',@(x) ischar(x));
    addOptional(p,'Verbose',false ,@(x) islogical(x));


    parse(p,data,Kvalues,varargin{:});
    
    nIterations=p.Results.NumIterations;
    trueLabels=p.Results.TrueLabels;
    distance=p.Results.Distance;
    verbose=p.Results.Verbose;
    
%% Iterative K-means clustering
    for k=Kvalues
        for i=1:nIterations
            [label_tmp,~] = kmeans(data, k, 'Distance',distance);
            s_tmp = silhouette(data, label_tmp, distance);
            silVals(i,k)=mean(s_tmp);
        end
    end
    
    [~,K_optimal]=max(median(silVals));
    fprintf('Found optimal K value = %d', K_optimal)
    median(silVals)
    std(silVals)
    
    if verbose
        ax=subplot(1,2,1);
        b=boxchart(silVals);
        b.BoxFaceColor='k';
        b.LineWidth=0.75;
        b.MarkerStyle='.';
        b.MarkerSize=10;
        b.MarkerColor='k';
        b.BoxWidth=0.75;
        ax.XLabel.String='K';
        ax.YLabel.String='Mean silhouette value';
        ax.Box='off';
        ax.LineWidth=0.75;
        ax.FontName='Arial';
        ax.FontSize=20;
        ax.XTickLabelRotation=0;
        ax.Title.String='K-means clustering';
    end

    
    %% Final K-means clustering with optimal K
    [Klabel,~] = kmeans(data,K_optimal,'Distance',distance,'Replicates',nIterations);
    
    if verbose
        ax=subplot(1,2,2);
        silhouette(data,Klabel,distance);
        ax.YLabel.String='Cluster ID';
        ax.Box='off';
        ax.LineWidth=0.75;
        ax.FontName='Arial';
        ax.FontSize=20;
        ax.Title.String='Silhouette plot - Optimal K';
    end
end