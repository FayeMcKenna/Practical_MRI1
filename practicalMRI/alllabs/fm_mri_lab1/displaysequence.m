function y = displaysequence(time,rfPulse,gradZ)
rfAmp = abs(rfPulse)
time = time
gradZ = abs(gradZ)
rfPhase = angle(rfPulse)
y = figure
subplot(4,1,1)%(1st plot is combination)
x1 = time; % first x axis is time
y1 = rfPhase; % first y axis is the phase of rfPulse
line(x1,y1,'color','r') % put in red
ax1 = gca; %current axes
ax1.XColor = 'r'; %red
ax1.YColor = 'r'; %red
ax1_pos = ax1.Position; % position of first axes in combo plot
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none'); % position 2nd plot in combo plot and in black
x2 = time;% second x axis is time
y2 = abs(rfPulse); % second y axis is amp of rfPulse
line(x2,y2,'Parent',ax2,'Color','g') % in green


ax3 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none'); % position 2nd plot in combo plot and in black
x3 = time;% second x axis is time
y3 = gradZ; % second y axis is amp of rfPulse
line(x3,y3,'Parent',ax3,'Color','b') % in blue
title('RF amplitude, RF phase and gradient amplitude') %title
ylabel('Tesla per meter ') % y label
%legend('rfamp','rfphase','gradamp'); % legend
xlabel('time'); % x label

% make other subplots
subplot(4,1,2) % second subplot in column showing rf amp and rf phase
x=time; % second plot x is time
y=rfPhase;% second plot y is phase of rfPulse
plot(x,y,'color','r') % plot red
title('RF phase') % title
ylabel('VF') % y measure label
xlabel('time'); % x axis label
legend('RF phase');% legend

subplot(4,1,3) % third subplot in column showing gradient amp
x=time; % third plot x is time
y=rfAmp;% third plot y is rfamp
plot(x,y,'color','g') % plot green
title('RF amplitude') % title
ylabel('%B1+ in Tesla') % y measure label
xlabel('time');% x measure label
legend('RF amplitude');% legend

subplot(4,1,4) % fourth subplot in column showing gradient amp
x=time; % fourth plot x is time
y=gradZ;% fourth plot y is amp of gamp
plot(x,y,'color','b') % plot blue
title('gradient amplitude') % title
ylabel('Tesla per meter') % y measure label
legend('gradAmp');% legend
xlabel('time');% x measure labe
end 
        