% f = 15000; % frequency
% t = -3:.01:3;
% a = sin(2*pi*t);
% b = sinc(t);
% % c = ramp(n);
% plot(t,b);
% plot(t,a);

for i=1:15000 %i starts at 1 go's to 15000
    t(i) = (-3:.01:3)*i;
    rfPulse(i) = sinc(t(i)); %B1+ in Tesla
    gradAmp(i) = 5*10^-2;                       %Tesla per meter
end
plot(t,rfPulse);



t = (-3:.00040001:3);
rfPulse = sinc(t); %B1+ in Tesla

plot(t,rfPulse);
