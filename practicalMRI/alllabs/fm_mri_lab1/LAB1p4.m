%% Exercise 1.4
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

%Allocate the memory needed
rfPulse     = zeros(1,10000); %variable to hold a RF waveform with 10000 samples
gradAmp     = zeros(1,10000); %variable to hold a gradient waveform with 10000 samples
time        = zeros(1,10000); %variable to hold 10000 time points

posZ       = zeros(1,400); %variable to hold 400 positions allong the z direction
mFinal     = zeros(1,400); %variable to hold the final magnetization calculated for each position

dB0AtPosJ  = 0 ; % klunge

%Generate the RF waveform (also generates the time line for plotting)
for i=1:10000 %i starts at 1 go's to 10000
    rfPulse(i) = exp(-((i-5000)/2000)^2)*10^-5;%B1+ in Tesla
    gradAmp(i) = 5*10^-2; %Tesla per meter
    time(i)    = i*10^-7;                       %Time in seconds
end


%Generate a list of sampling points allong the z direction
for i=1:400 %i starts at 1 go's to 200
    posZ(i)  = (i-200)*10^-4; %Distance from iso center in meters
end  


%Question A: (plot the wave form)

%plot in one figure RF amplitude, RF phase and gradient amplitude by time,
figure
subplot(3,1,1) % RF amplitude subplot 
x=time; % second plot x is time
y=abs(rfPulse);% second plot y is mag of rfPulse
plot(x,y,'color','r') % plot red
title('RF Amplitude') % title
ylabel('Tesla') % y measure label
xlabel('time'); % x axis label
legend('rfamp');% legend

subplot(3,1,2) % RF phase subplot 
x=time; % third plot x is time
y=angle(rfPulse);% third plot y is phase
plot(x,y,'color','g') % plot green
title('RF Phase') % title
ylabel('VF') % y measure label
xlabel('time');% x measure label
legend('RFphase');% legend

subplot(3,1,3) % gradient amplitude subplot 
x=time; % fourth plot x is time
y=gradAmp;% fourth plot y is amp 
plot(x,y,'color','b') % plot blue
title('Gradient Amplitude') % title
ylabel('Tesla') % y measure label
legend('GradAmp');% legend
xlabel('time');% x measure label

%Question B: (compute the magnetization / time)

for j=1:400 %loop over different positions allong the z direction
    m =  smalltipangle(dB0AtPosJ,rfPulse(1), 0); % first time point assuming the transverse magnetization 0
    for i=2:10000 %(i starts at 2) for every point in the magnetization process
        dB0AtPosJ = posZ(j)*gradAmp(i); % dB0 of j is its pos * its gradient
         m =  smalltipangle(dB0AtPosJ,rfPulse(i), m);% total magentization is updated for one time in rfPulse
    end
    mFinal(j) = m; %magnetization calculated for each j point
end    

 
%Question C: Plot in a single figure (two subplots underneath one another) the magnitude and phase of
%the mr signal as a function of the position along the z direction
figure
subplot(2,1,1) % RF amplitude subplot 
plot(posZ,abs(mFinal),'color','r') % plot red
title('Magnetitude of MR signal') % title
ylabel('Tesla') % y measure label
xlabel('posZ'); % x axis label
legend('Magnetitude of MR signal');% legend

subplot(2,1,2) % RF phase subplot 
plot(posZ,angle(mFinal),'color','g') % plot green
title('Phase of MR signal') % title
ylabel('VF') % y measure label
xlabel('posZ');% x measure label
legend('Phase of MR signal');% legend

% %Question D:

% how much signal would this excitation produce?
 abs(sum(mFinal)) 
% This excitation would produce a signal of .0163.  The plotted slice 
% profile shows the magnetization amp to be a max of 0.14,with a range
% below that so it seems about right.
