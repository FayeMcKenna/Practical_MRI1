function y = profplot(mxy,mz,posZ)
mxy_mag = abs(mxy);
mxy_ph = angle(mxy);
mz_mag = mz;   

%plot one magnitude of mxy 
subplot(3,1,1) % second subplot in column showing rf amp and rf phase
x=posZ; % second plot x is time
y=mxy_mag;% second plot y is phase of rfPulse
plot(x,y,'color','r') % plot red
title('mxy magnitude') % title
ylabel('Tesla') % y measure label
xlabel('PosZ'); % x axis label
legend('mx mag');% legend

subplot(3,1,2) % third subplot in column showing gradient amp
x=posZ; % third plot x is time
y=mxy_ph;% third plot y is rfamp
plot(x,y,'color','g') % plot green
title('mxy phase') % title
ylabel('VF') % y measure label
xlabel('PosZ');% x measure label
legend('mxy phase');% legend

subplot(3,1,3) % fourth subplot in column showing gradient amp
x=posZ; % fourth plot x is time
y=mz_mag;% fourth plot y is amp of gamp
plot(x,y,'color','b') % plot blue
title('mz magnitude') % title
ylabel('Tesla') % y measure label
legend('mz mag');% legend
xlabel('PosZ');% x measure labe