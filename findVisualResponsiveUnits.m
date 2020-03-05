function response=findVisualResponsiveUnits(raster)


    for suIdx=1:size(raster,1)

        baseline=@(x) nnz(x>-500&x<0);
        early=@(x) nnz(x>1&x<500);
        late=@(x) nnz(x>501&x<1000);

        sumSpikes(:,1)=cellfun(baseline,raster(suIdx,:))';
        sumSpikes(:,2)=cellfun(early,raster(suIdx,:))';
        sumSpikes(:,3)=cellfun(late,raster(suIdx,:))';

        %Early responders
        response(suIdx,1)=doStatistics(sumSpikes(:,1),sumSpikes(:,2),0.05);
        %Late responders
        response(suIdx,2)=doStatistics(sumSpikes(:,1),sumSpikes(:,3),0.05);

    end
    
end


function outcome=doStatistics(baseline,response,alpha)
    [~,p]=ttest(baseline,response);
    
    outcome=0;
    if p<=alpha
        if sum(response)>sum(baseline) && sum(response) > size(response,1)*.3  
            outcome=1;
        elseif sum(response)<sum(baseline) && sum(baseline) > size(baseline,1)*.3  
            outcome=2;
        end    
    end
    
end