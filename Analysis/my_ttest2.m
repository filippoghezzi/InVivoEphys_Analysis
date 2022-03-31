function [h,p,stats] = my_ttest2 (x,y)
    
   
    if swtest(x) && swtest(y)
        [p,h,stats]=ranksum(x,y);
        fprintf('Wilcoxon rank sum test...\n')
    else 
        [h,p,~,stats]=ttest2(x,y);
        fprintf('Two sample t-test...\n')
    end

    
end