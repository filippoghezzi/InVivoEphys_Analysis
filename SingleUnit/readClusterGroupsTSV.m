function [ID,group]=readClusterGroupsTSV(filename)

    s=tdfread(filename);
    ID=(s.cluster_id)';
    
    for i=1:size(s.group,1)
        switch s.group(i,:)
            
            case 'good '
                group(1,i)=2;
                
            case 'noise'
                group(1,i)=0;
                
            case 'mua  '
                group(1,i)=1;
                
            otherwise
                group(1,i)=3;
        end
    end
    
    
end
            
                
    