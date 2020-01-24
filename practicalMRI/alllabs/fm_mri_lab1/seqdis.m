function y = displaysequence(time,rfPulse,gradZ)
rfAmp = abs(rfPulse)
rfPhase = angle(rfPulse)
y = figure

%plot one magnitude of mxy 
subplot(3,1,1) % second subplot in column showing rf amp and rf phase
plot(time,rfAmp,'color','r') % plot red
title('rf amplitude') % title
ylabel('') % y measure label
xlabel(''); % x axis label
legend('');% legend

subplot(3,1,2) % second subplot in column showing rf amp and rf phase
plot(time,rfPhase,'color','g') % plot red
title('rf magnitude') % title
ylabel('angle') % y measure label
xlabel('time'); % x axis label


subplot(3,1,3) % second subplot in column showing rf amp and rf phase
plot(time,gradZ,'color','g') % plot red
title('rf magnitude') % title
ylabel('tesla/meter') % y measure label
xlabel('time'); % x axis label

