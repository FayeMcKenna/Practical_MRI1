%% Exercise 5.1
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same


dt    = 10^-5; 
gamma = 42.577*10^6;

te = [6 8 12 16 24 32 48 64 96 128 192 256]*10^-3; % te in ms

% set T1 and T2
T1 = 850;
T2 = 60;

%Question A

%Allocate the memory needed
nTimeSteps  = 70000;
rfPulseE    = zeros(1,200); %exciation
rfPulseR    = zeros(1,200); %refcosuing
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform with 10000 samples
gradAmp     = zeros(3,nTimeSteps); %variable to hold a gradient waveform with 10000 samples
time        = zeros(1,nTimeSteps); %variable to hold 10000 time points

tiSteps    = size(te,2);
nSpins     = 1000;
spins      = zeros(1,nSpins); %variable to hold nPosSteps positions allong the z direction
mFinalVect = zeros(tiSteps,2); %variable to hold the final magnetization calculated for each position
tmp        = zeros(nSpins,2);

%Generates the time line for plotting
for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = i*dt;                       %Time in seconds
end

%Generate off-resonances 
for i=1:nSpins
    spins(i) = 2.5*(rand-0.5)*10^-7;     
end

%Question A:

%Generate the excitation pulse 1ms 3TBW Sinc pulse.
for i=1:100 %i 
    rfPulse1(i) = (sin(pi*i/100)^2)*sinc(pi*(i-35)/60)*10^-5; %B1+ in Tesla          
end

%Generate the RF excitation waveform to be 90 degrees
 A = sum(rfPulse1);
rfPulse_FA = gamma*A*dt; % equation where flip angle=gamma*sum of rfPulse*dt

degree = rfPulse_FA*360; % degree = 25.5194 %rfPulse2 = rfPulse*pi/rfPulse_FA;

rfPulse1=90/degree*rfPulse1; %add on remaining degree that pulse needs to generate for 180 degrees

% use Hann window funtion to apodize excite and invert rfPulses

l = 100 ;%lenght of the pulse

h = hann(l); % use hann function

hfunc = transpose(h); % transpose to match rfPulse array

rfPulse1 = times(hfunc,rfPulse1); % create the hann function for inversion rfPulse

%Generate the refocusing waveform 1ms 3TBW Sinc pulse
for i=1:100 %i 
    rfPulse2(i) = (sin(pi*i/100)^2)*sinc(pi*(i-35)/60)*10^-5; %B1+ in Tesla          
end

%Generate the RF Refocusing waveform to be 180 degrees
A = sum(rfPulse2);
rfPulse_FA = gamma*A*dt; % equation where flip angle=gamma*sum of rfPulse*dt

degree = rfPulse_FA*360; % degree = 25.5194 %rfPulse2 = rfPulse*pi/rfPulse_FA;

rfPulse2=180/degree*rfPulse2; %add on remaining degree that pulse needs to generate for 180 degrees

% use Hann window funtion to apodize excite and invert rfPulses

l = 100 ;%lenght of the pulse

h = hann(l); % use hann function

hfunc = transpose(h); % transpose to match rfPulse array

rfPulse2 = times(hfunc,rfPulse2); % create the hann function for inversion rfPulse

%after pulses 1 ms gap
for i=1:100 %i 
    rfPulse3(i) = 0; %B1+ in Tesla          
end

%add pulse parts together
 rfPulse = [rfPulse1, rfPulse3,rfPulse1, rfPulse3];

%% Question B:

for measIndex = 1:size(te,2)
T1 = 850;
T2 = 60;
    for spinIndex=1:nSpins
     gradAtPosJ = spins(spinIndex);
        [mT,mZ] =  bloch(dt, gradAtPosJ,rfPulse(1), 850*1.0e-3, 60*1.0e-3, 0, 1);
        for i=2:tiSteps %i starts at 2
                [mT,mZ] =  bloch(dt, gradAtPosJ,rfPulse(i), 850*1.0e-3, 60*1.0e-3, mT, mZ);
         end

    tmp(spinIndex,:) = [mT,mZ];
    end
    mFinalVect(measIndex,:) = [mT, mZ];


end

mz = mFinalVect(:,2); % mz array
mxy = mFinalVect(:,1); % mx array


%% Plot the measured signal (amplitude and phase) as a function of TE.

figure
profplot(mxy,mz,te); 

mxy=abs(mxy)
t2=zeros(1,12);
tet= transpose(te);

%t2mxy1= tet/mxy1;

for j=mxy
    t2mxy= exp(-tet/mxy); 
end

t2mxy = t2mxy(:,3);

% Plot the measured signal (amplitude and phase) as a function of T2.
figure
profplot(t2mxy,mz,tet); 
title('T2 ') % title
ylabel(' T2') % y measure label
xlabel('TE');% x measure label
legend('t2');% legend

%Fit the appropriate signal equation to the measured data points to T2

%fit data
f = fit(t2mxy, tet, 'poly2' );
%Plot your fit and the data.
figure
plot( f, t2mxy, tet )
title('T2 ') % title
ylabel(' T2') % y measure label
xlabel('TE');% x measure label
legend('t2');% legend



