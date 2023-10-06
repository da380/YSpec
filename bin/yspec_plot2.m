function yspec_plot2(t1,t2,i1,i2,comp,frac)

t1 = t1/60;
t2 = t2/60;

thick = 0.1;

pref = 'seis.out.';



% open the first file 
str=int2str(i1);
file=strcat(pref,str);  
d1=dlmread(file);

nt = length(d1(:,1));
nr = i2-i1+1;

d(1:nt,1:5,1:nr) = 0;

for i=i1:i2
    str=int2str(i);
    file=strcat(pref,str);
    clear d1
    d1=dlmread(file);
    for j = 1:nt
        for k = 1:5
           d(j,k,i-i1+1) = d1(j,k);
        end
    end
end



dt = d(2,1,1)-d(1,1,1);

t = 0;
n = 0;
while (t < t2*3600)
    n = n+1;
    t = t+dt;
end
n = min([n,length(d(:,1,1))]);




frac2 = 0.0;
t11 = max(frac2*t2,t1+0.1);
t22 = (1-frac2)*t2;

% apply a taper in the time domain
for i = 1:n
    t = (i-1)*dt/3600;
    if(t < t1)
        filt = 0;
    elseif(t >= t1 & t < t11)
        filt = pi*(t-t1)/(t11-t1);
        filt = 0.5*(1-cos(filt));
    elseif(t >= t1 & t <= t22)
        filt = 1;
    elseif(t > t22 & t < t2)
        filt = pi*(t2-t)/(t2-t22);
        filt = 0.5*(1-cos(filt));
    elseif(t > t2)
        filt = 0;
    end
    d(i,2:5,:) = filt*d(i,2:5,:);
end

m = max(max(max(abs(d(1:n,comp+1,:)))));

close; hold on
for k = 1:nr
    plot(d(:,1,1)/60,d(:,comp+1,k)+frac*(k-1)*m,...
        'linewidth',thick,'color','k')
end
xlim([t1*60,t2*60])
ylim([-m,m*nr*frac])
end
