%% Exercise 1.5
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

%Allocate the memory needed
rfPulse     = zeros(1,15000); %variable to hold a RF waveform with 10000 samples
gradAmp     = zeros(1,15000); %variable to hold a gradient waveform with 10000 samples
time        = zeros(1,15000); %variable to hold 10000 time points


posZ       = zeros(1,400); %variable to hold 400 positions allong the z direction
mFinal     = zeros(1,400); %variable to hold the final magnetization calculated for each position

dB0AtPosJ  = 0 ; % klunge

%Generate the RF waveform
for i=1:15000 %i starts at 1 go's to 10000
    rfPulse(i) = exp(-((i-5000)/2000)^2)*10^-5; %B1+ in Tesla
    gradAmp(i) = 5*10^-2;                       %Tesla per meter
end

%Generates the time line for plotting
for i=1:15000 %i starts at 1 go's to 15000
    time(i)    = i*10^-7;                       %Time in seconds
end

%Generate a list of sampling points allong the z direction
for i=1:400 %i starts at 1 go's to 200
    posZ(i)  = (i-200)*10^-4; %Distance from iso center in meters
end  

% Question A: 
% add a slice re-focusing gradient to the wave form 

for i=10001:15000 %i starts at 1 go's to 10000
    %rfPulse(i) = rfPulse(i)*0; %B1+ in Tesla
    gradAmp(i) = - 5*10^-2; % applying opposite
end

% Question B: 
% plot in one figure RF amplitude, RF phase and gradient amplitude by time,
% and subplots in column
figure
% make other subplots
subplot(3,1,1) % second subplot in column showing rf amp and rf phase
plot(time,angle(rfPulse),'color','r') % plot red
title('RF phase') % title
ylabel('VF') % y measure label
xlabel('time'); % x axis label
legend('RF phase');% legend

subplot(3,1,2) % third subplot in column showing gradient amp
plot(time,abs(rfPulse),'color','g') % plot green
title('RF amplitude') % title
ylabel('%B1+ in Tesla') % y measure label
xlabel('time');% x measure label
legend('RF amplitude');% legend

subplot(3,1,3) % fourth subplot in column showing gradient amp
plot(time,gradAmp,'color','b') % plot blue
title('gradient amplitude') % title
ylabel('Tesla per meter') % y measure label
legend('gradAmp');% legend
xlabel('time');% x measure label
 
 %Question C: % simulate the slice profile

for j=1:400 %loop over different positions allong the z direction
    m =  smalltipangle(dB0AtPosJ,rfPulse(1), 0); % first time point assuming the transverse magnetization 0
    for i=2:15000 %(i starts at 2) for every point in the magnetization process
        dB0AtPosJ = posZ(j)*gradAmp(i);% dB0 of j is its pos * its gradient
         m =  smalltipangle(dB0AtPosJ,rfPulse(i), m);% total magentization is updated for one time in rfPulse
    end
    mFinal(j) = m; %magnetization calculated for each j point
end    

% % % Question D: Plot in a single figure (two sub plots underneath one another) the magnitude and phase of
% the mr signal as a function of the position along the z direction

figure
subplot(2,1,1)
plot(posZ,abs(mFinal));
title('Magnitude') % title
ylabel('Tesla') % y measure label
legend('magnititude');% legend
xlabel('position');% x measure label

subplot(2,1,2);
plot(posZ,angle(mFinal));
title('Phase') % title
ylabel('VF') % y measure label
legend('phase');% legend
xlabel('position');% x measure label


% Question E: How much signal does the excitation produce now? 
%Why is there more signal now than in Exercise 1.4: There is more signal
%because the the refocusing gradient has realligned the phases so there is
%more signal detected.   12.5981

abs(sum(mFinal))

% Question F : What is the slice thickness produced by this excitation? Is this a nice slice profile?

%slice thickness 
k = max(mFinal)*.5; % the maximum of the final magnetization /2 (for area)

thick =find(abs(k)> 0); % find where within that area is greater than 0
thick = sum(thick); % sum together

%The slice thickness produced by this excitation is 1.

% slice profile is a pretty nice slice excitation profile, ideally it would be rectangular but that is impossible for the scanner
%this one has slight slopes.
figure
plot(posZ,abs(mFinal));







