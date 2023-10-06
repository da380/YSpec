clear all






f11 = 0.2;
f12 = 1;
f21 = 160.0;
f22 = 166.0;

t11 = 1000/f11;
t12 = 1000/f12;
t21 = 1000/f21;
t22 = 1000/f22;


tf = 1000;

close

%file = 'flores.syn';

file = 'yspec.out';

irec1 = 5;
irec2 = 85;
irecs = 5;

icomp = 1;

lat11 = 2;
dlat = 2;
lat1 = lat11+(irec1-1)*dlat;
lat2 = lat11+(irec2-1)*dlat;



for irec = irec1:irecs:irec2
    

   irecp = 3*(irec-1)+icomp;
    
   [dat,t] = ah_get_data(strcat(file,'.ahx'),irecp,tf,t21,t22,t11,t12);
   
   
   [dat2,t2] = ah_get_data(strcat(file,'.ahx.syn'),irecp,tf,t21,t22,t11,t12);

   fac =  0.1;  
   if(irec == irec1)
       mm = fac*max(dat)/dlat;
   end

   dat = dat/mm;
   dat2 = dat2/mm;
   lat = lat1+(irec-irec1)*dlat;
   
   t = t/60;
   t2 = t2/60;

   hold on
   plot(t,dat+lat)
   plot(t2,dat2+lat,'r--')
   xlim([min(t),max(t)])

    
    
    
end




