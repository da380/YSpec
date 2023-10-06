function yspec_plot(t1,t2,i1)

t1 = t1/60;
t2 = t2/60;

str=int2str(i1);
file=strcat('seis.out.',str);  
d=dlmread(file);



close;

dt = d(2,1)-d(1,1);

t = 0;
n = 0;
while (t < t2*3600)
    n = n+1;
    t = t+dt;
end
n = min([n,length(d(:,1))]);

frac = 0.1;
t11 = max(frac*t2,t1+0.1);
t22 = (1-frac)*t2;

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
    d(i,2:5) = filt*d(i,2:5);
end


m = 1.1*max(max(abs(d(1:n,2:4))));
mg = 1.1*(max(abs(d(1:n,5))));

subplot(4,1,1)
hold on
plot(d(:,1),d(:,2),'k')
%xlim([t1*60,t2*60])
%ylim([-m,m])

subplot(4,1,2)
hold on
plot(d(:,1),d(:,3),'k')
%xlim([t1*60,t2*60])
%ylim([-m,m])

subplot(4,1,3)
hold on
plot(d(:,1),d(:,4),'k')
%xlim([t1*60,t2*60])
%ylim([-m,m])


subplot(4,1,4)
hold on
plot(d(:,1),d(:,5),'k')
%xlim([t1*60,t2*60])
ylim([-mg,mg])

end
