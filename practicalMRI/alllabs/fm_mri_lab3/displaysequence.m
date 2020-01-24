function y = displaysequence(time,rfPulse,gradZ)
rfAmp = abs(rfPulse);
time = time;
gradZ = gradZ;
rfPhase = angle(rfPulse);
y = figure

% make other subplots
subplot(3,1,1) % second subplot in column showing rf amp and rf phase
x=time; % second plot x is time
y=rfPhase;% second plot y is phase of rfPulse
plot(x,y,'color','r') % plot red
title('RF phase') % title
ylabel('VF') % y measure label
xlabel('time'); % x axis label
legend('RF phase');% legend

subplot(3,1,2) % third subplot in column showing gradient amp
x=time; % third plot x is time
y=rfAmp;% third plot y is rfamp
plot(x,y,'color','g') % plot green
title('RF amplitude') % title
ylabel('%B1+ in Tesla') % y measure label
xlabel('time');% x measure label
legend('RF amplitude');% legend

subplot(3,1,3) % fourth subplot in column showing gradient amp
x=time; % fourth plot x is time
y=gradZ;% fourth plot y is amp of gamp
plot(x,y,'color','b') % plot blue
title('gradient amplitude') % title
ylabel('Tesla per meter') % y measure label
legend('gradAmp');% legend
xlabel('time');% x measure labe
end 
        