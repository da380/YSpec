clear all

irec = 40;




f11 = 1;
f12 = 2;
f21 = 160.0;
f22 = 166.0;

t11 = 1000/f11;
t12 = 1000/f12;
t21 = 1000/f21;
t22 = 1000/f22;


tf = 100;

close

%file = 'flores.syn';

file = 'yspec.out';



for icomp = 1:3

    

    irecp = 3*(irec-1)+icomp;
    
   [dat,t] = ah_get_data(strcat(file,'.ahx'),irecp,tf,t21,t22,t11,t12);
   
   
   %[dat2,t2] = ah_get_data(strcat(file,'.ahx.syn'),irecp,tf,t21,t22,t11,t12);

   
   t = t/60;
   %t2 = t2/60;
   fac = 1;
   
   subplot(3,1,icomp)
   hold on
   plot(t,dat)
   %plot(t2,dat2,'r--')
   %plot(t2,fac*(dat2-dat)-1.2*max(abs(dat)),'g')
   xlim([min(t),max(t)])

end

