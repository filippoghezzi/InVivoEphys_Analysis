function response=getUnitResponsiveness(raster,condition)
%   function response=findVisualResponsiveUnits(raster)
%   
%   To determine if a cell is responsive, t-test or anova or kruskalwallis test
%   to compare the number of spikes within a given time window between conditions 
%   (baseline, ON, OFF and late response). Additional condition on reliability of firing, e.g. number
%   of trials with spikes needs to be at least cutoff of total trial number. 
%
%   Output: response 
%            value:  1 -> increse firing
%                    2 -> decrease firing
%            column: 1 -> ON
%                    2 -> OFF
%                    3 -> late



%% Set anonymous functions

    
    if strcmp(condition,'Visual') || strcmp(condition,'VisualOpto') || strcmp(condition,'Visual_K')
        baselineWindow=[-100,0];
        targetWindow=[(0:100:1900)'+1,(100:100:2000)'];
%         baseline=@(x) nnz(x>-100&x<0);
%         on=@(x) nnz(x>0&x<100);
%         off=@(x) nnz(x>200&x<300);
%         late=@(x) nnz(x>800&x<900);
        cutoff=0;
        
    elseif strcmp(condition,'WhiskerStim') || strcmp(condition,'WhiskerStim_K') 
        baselineWindow=[-100,-50];
        targetWindow=[0,50];
%         baseline=@(x) nnz(x>-100&x<0);
%         on=@(x) nnz(x>0&x<100);
%         off=@(x) nnz(x>200&x<300);
%         late=@(x) nnz(x>800&x<900);
        cutoff=0;
        
    elseif strcmp(condition,'Optotagging') 
        baselineWindow=[-50,0];
        targetWindow=[0,50];
%         baseline=@(x) nnz(x>-50&x<0);
%         on=@(x) nnz(x>0&x<50);
        cutoff=0.2;

    elseif strcmp(condition,'LaserOnly')
        baselineWindow=[-250,-50];
        targetWindow=[-50,150;150,350];
%         baseline=@(x) nnz(x>-250&x<-50);
%         on=@(x) nnz(x>-50&x<150);
%         off=@(x) nnz(x>150&x<350);
        cutoff=0.3;

    end
        
    for suIdx=1:size(raster,1)
        %% Build spike number x trial matrix
        baseline=@(x) nnz(x>baselineWindow(1)&x<baselineWindow(2));
        sumSpikes(:,1)=cellfun(baseline,raster(suIdx,:))';
        
        for responseWindow=1:size(targetWindow,1)
            target=@(x) nnz(x>targetWindow(responseWindow,1)&x<targetWindow(responseWindow,2));
            sumSpikes(:,1+responseWindow)=cellfun(target,raster(suIdx,:))';
        end
%         if exist ('off','var')
%             sumSpikes(:,3)=cellfun(off,raster(suIdx,:))';
%         end
%         if exist ('late','var')
%             sumSpikes(:,4)=cellfun(late,raster(suIdx,:))';
%         end
        
        %% Do statistics
        if size(sumSpikes,2) == 2
            response(suIdx,1)=doTtest(sumSpikes(:,1),sumSpikes(:,2),0.05,cutoff);
        elseif size(sumSpikes,2) > 2
            response(suIdx,:)=doAnova(sumSpikes,0.05,cutoff);
        end
        
    end
    
end


function outcome=doTtest(baseline,response,alpha,cutoff)
    
    if checkNormality([baseline,response])
        [p,~] = ranksum(baseline,response);
    else
        [~,p]=ttest(baseline,response);
    end   
    
    outcome=0;
    if p<=alpha
        if sum(response)>sum(baseline) && sum(response) > size(response,1)*cutoff  
            outcome=1;
        elseif sum(response)<sum(baseline) 
            outcome=2;
        end    
    end
    
end


function outcome=doAnova(data,alpha,cutoff)
% Only vs. Baseline mulitple comparison. Baseline defined as first column.

    if checkNormality(data)
        [p,~,stats]=kruskalwallis(data,[],'off');
    else
        [p,~,stats]=anova1(data,[],'off');
    end

    outcome=zeros(1,size(data,2)-1);
    if p <= alpha
        c = multcompare(stats,'Display','off');

        for i=1:size(data,2)-1
            if c(i,6) <= alpha            
                if sum(data(:,i+1)) > sum(data(:,1)) && sum(data(:,i+1)) > size(data,1)*cutoff 
                    outcome(1,i)=1;
                elseif sum(data(:,i+1)) < sum(data(:,1))
                    outcome(1,i)=2;
                end 
            end  
        end            
    end
        
end


function h=checkNormality(data)
% h, false, normal distribution
% h, true, non-normal distrubition
    for group=1:size(data,2)
        h(group)=lillietest(data(:,group));
    end
    if sum(h)==0; h=false; else; h=true; end
end