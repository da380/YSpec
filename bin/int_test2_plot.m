function int_test2_plot(ilog)

     
    d = dlmread('int_test2.out1');
    d2 = dlmread('int_test2.out2');
    d3 = dlmread('int_test2.out3');
    
    f = d(:,1);
   
    
    if(ilog == 1) 
         
       
        u1 = log(sqrt(d(:,2).^2+d(:,3).^2));   
        v1 = log(sqrt(d(:,4).^2+d(:,5).^2));    
    
        u2 = log(sqrt(d(:,6).^2+d(:,7).^2));    
        v2 = log(sqrt(d(:,8).^2+d(:,9).^2));                
    
        u3 = log(sqrt(d(:,10).^2+d(:,11).^2));   
        v3 = log(sqrt(d(:,12).^2+d(:,13).^2));            
    
        u4 = log(sqrt(d(:,14).^2+d(:,15).^2));    
        v4 = log(sqrt(d(:,16).^2+d(:,17).^2));
        
        u12 = log(sqrt(d2(:,2).^2+d2(:,3).^2));   
        v12 = log(sqrt(d2(:,4).^2+d2(:,5).^2));    
    
        u22 = log(sqrt(d2(:,6).^2+d2(:,7).^2));    
        v22 = log(sqrt(d2(:,8).^2+d2(:,9).^2));                
   
        u32 = log(sqrt(d2(:,10).^2+d2(:,11).^2));   
        v32 = log(sqrt(d2(:,12).^2+d2(:,13).^2));            
    
        u42 = log(sqrt(d2(:,14).^2+d2(:,15).^2));    
        v42 = log(sqrt(d2(:,16).^2+d2(:,17).^2));
        
        u13 = log(sqrt(d3(:,2).^2+d3(:,3).^2)); 
        v13 = log(sqrt(d3(:,4).^2+d3(:,5).^2));    
    
        u23 = log(sqrt(d3(:,6).^2+d3(:,7).^2));    
        v23 = log(sqrt(d3(:,8).^2+d3(:,9).^2));                
   
        u33 = log(sqrt(d3(:,10).^2+d3(:,11).^2));   
        v33 = log(sqrt(d3(:,12).^2+d3(:,13).^2));            
    
        u43 = log(sqrt(d3(:,14).^2+d3(:,15).^2));    
        v43 = log(sqrt(d3(:,16).^2+d3(:,17).^2));
      
    end
    
       if(ilog == 0) 
   
        u1 = (sqrt(d(:,2).^2+d(:,3).^2));   
        v1 = (sqrt(d(:,4).^2+d(:,5).^2));    
    
        u2 = (sqrt(d(:,6).^2+d(:,7).^2));    
        v2 = (sqrt(d(:,8).^2+d(:,9).^2));                
    
        u3 = (sqrt(d(:,10).^2+d(:,11).^2));   
        v3 = (sqrt(d(:,12).^2+d(:,13).^2));            
    
        u4 = (sqrt(d(:,14).^2+d(:,15).^2));    
        v4 = (sqrt(d(:,16).^2+d(:,17).^2));
     
        u12 = (sqrt(d2(:,2).^2+d2(:,3).^2));   
        v12 = (sqrt(d2(:,4).^2+d2(:,5).^2));    
    
        u22 = (sqrt(d2(:,6).^2+d2(:,7).^2));    
        v22 = (sqrt(d2(:,8).^2+d2(:,9).^2));                
    
        u32 = (sqrt(d2(:,10).^2+d2(:,11).^2));   
        v32 = (sqrt(d2(:,12).^2+d2(:,13).^2));            
    
        u42 = (sqrt(d2(:,14).^2+d2(:,15).^2));    
        v42 = (sqrt(d2(:,16).^2+d2(:,17).^2));
        
      
        u13 = (sqrt(d3(:,2).^2+d3(:,3).^2)); 
        v13 = (sqrt(d3(:,4).^2+d3(:,5).^2));    
    
        u23 = (sqrt(d3(:,6).^2+d3(:,7).^2));    
        v23 = (sqrt(d3(:,8).^2+d3(:,9).^2));                
   
        u33 = (sqrt(d3(:,10).^2+d3(:,11).^2));   
        v33 = (sqrt(d3(:,12).^2+d3(:,13).^2));            
    
        u43 = (sqrt(d3(:,14).^2+d3(:,15).^2));    
        v43 = (sqrt(d3(:,16).^2+d3(:,17).^2));
        
       end
    

       
    close
    
    subplot(4,2,1)
    hold on 
    plot(f,u1)
    plot(f,u12,'r--')
    plot(f,u13,'g-.')
    xlim([min(f),max(f)])
    subplot(4,2,2)
    hold on
    plot(f,v1)
    plot(f,v12,'r--')
    plot(f,v13,'g-.')
    xlim([min(f),max(f)])
    
    subplot(4,2,3)
    hold on
    plot(f,u2)
    plot(f,u22,'r--')
    plot(f,u23,'g-.')
    xlim([min(f),max(f)])
    subplot(4,2,4)
    hold on
    plot(f,v2)
    plot(f,v22,'r--')
    plot(f,v23,'g-.')
    xlim([min(f),max(f)])
    
    
    subplot(4,2,5)
    hold on
    plot(f,u3)
    plot(f,u32,'r--')
    plot(f,u33,'g-.')
    xlim([min(f),max(f)])
    subplot(4,2,6)
    hold on
    plot(f,v3)
    plot(f,v32,'r--')
    plot(f,v33,'g-.')
    xlim([min(f),max(f)])
   
    subplot(4,2,7)
    hold on
    plot(f,u4)
    plot(f,u42,'r--')
    plot(f,u43,'g-.')
    xlim([min(f),max(f)])
    subplot(4,2,8)
    hold on
    plot(f,v4)
    plot(f,v42,'r--')
    plot(f,v43,'g-.')
    xlim([min(f),max(f)])
   
    
    
    


end