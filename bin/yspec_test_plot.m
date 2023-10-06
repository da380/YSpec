function yspec_test_plot(iplot,ilog)



if(iplot == 1)
    
   d = dlmread('yspec_test.rad');
   
   f = d(:,1);
   if(ilog == 1) 
      u1 = log(sqrt(d(:,2).^2+d(:,3).^2));
      u2 = log(sqrt(d(:,4).^2+d(:,5).^2));
   else
      u1 = (sqrt(d(:,2).^2+d(:,3).^2));
      u2 = (sqrt(d(:,4).^2+d(:,5).^2));
   end
   
   close
   subplot(2,1,1)
   plot(f,u1)
   xlim([min(f),max(f)])
   subplot(2,1,2)
   plot(f,u2)
   xlim([min(f),max(f)])
   
end


if(iplot == 2)
    
   d = dlmread('yspec_test.tor');
   
   f = d(:,1);
   
   if(ilog == 1)
   w1 = log(sqrt(d(:,2).^2+d(:,3).^2));
   w2 = log(sqrt(d(:,4).^2+d(:,5).^2));
   end
   
   
   if(ilog == 0)
   w1 = (sqrt(d(:,2).^2+d(:,3).^2));
   w2 = (sqrt(d(:,4).^2+d(:,5).^2));
   end
   
   close
   subplot(2,1,1)
   plot(f,w1)
   xlim([min(f),max(f)])
   subplot(2,1,2)
   plot(f,w2)
   xlim([min(f),max(f)])
    
end


if(iplot == 3)
    
     
    d = dlmread('yspec_test.sph');
    
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
        
    end
    
    close
    
    subplot(4,2,1)
    plot(f,u1)
    xlim([min(f),max(f)])
    subplot(4,2,2)
    plot(f,v1)
    xlim([min(f),max(f)])
    
    subplot(4,2,3)
    plot(f,u2)
    xlim([min(f),max(f)])
    subplot(4,2,4)
    plot(f,v2)
    xlim([min(f),max(f)])
    
    
    subplot(4,2,5)
    plot(f,u3)
    xlim([min(f),max(f)])
    subplot(4,2,6)
    plot(f,v3)
    xlim([min(f),max(f)])
    
    subplot(4,2,7)
    plot(f,u4)
    xlim([min(f),max(f)])
    subplot(4,2,8)
    plot(f,v4)
    xlim([min(f),max(f)])
    
    
    
    
    
    
end


end