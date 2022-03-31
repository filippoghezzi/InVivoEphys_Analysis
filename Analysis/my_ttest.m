function [h,p,stats] = my_ttest (x,varargin)
    
    if nargin==1
        if swtest(x)
            [p,h,stats]=signrank(x);
            fprintf('Wilcoxon signed rank test...\n')
        else 
            [h,p,~,stats]=ttest(x);
            fprintf('One sample t-test...\n')
        end
    end
    
    if nargin==2
        y=varargin{1,1};
        if swtest(x) && swtest(y)
            [p,h,stats]=signrank(x,y);
            fprintf('Wilcoxon signed rank test...\n')
        else 
            [h,p,~,stats]=ttest(x,y);
            fprintf('Paired t-test...\n')
        end
    end
    
    
end