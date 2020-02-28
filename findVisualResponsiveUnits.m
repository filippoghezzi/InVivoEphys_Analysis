function findVisualResponsiveUnits(raster)


for suIdx=1:size(raster,1)
    
    
            baseline=@(x) nnz(x>-500&x<0);
            early=@(x) nnz(x>1&x<500);
            late=@(x) nnz(x>501&x<1000);
        
            sumSpikes(:,1)=cellfun(baseline,raster(suIdx,:))';
            sumSpikes(:,2)=cellfun(early,raster(suIdx,:))';
            sumSpikes(:,3)=cellfun(late,raster(suIdx,:))';
            
            %Early responders
            [~,p]=ttest(sumSpikes(:,1),sumSpikes(:,2));
            if p<=0.05 && sum(sumSpikes(:,2))>size(sumSpikes,1)*.3       
                if sum(sumSpikes(:,2))>sum(sumSpikes(:,1))
%                     suModV(suIdx,session)=1;
                else
%                     suModV(i,session)=2;
                end    
            end
end
end