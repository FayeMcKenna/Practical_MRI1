% Use an 1ms apodized sinc, time bandwidth factor 3, 
%to excite a 5mm thick transverse slice (z-gradient) using a flip-angle of 90 .
%Use 0.5ms for the slice refocusing gradient

%% Exercise 3.3
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same


dt    = 10^-6; 
gamma = 42.577*10^6;

ti = [25 50 100 200 400 800 1600 3200 6400]*10^-3; % ti in ms

%Allocate the memory needed
nTimeSteps  = 30000;
rfPulseI    = zeros(1,400); %inversion
rfPulseE    = zeros(1,200); %exciation
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform with 10000 samples
gradAmp     = zeros(1,nTimeSteps); %variable to hold a gradient waveform with 10000 samples
time        = zeros(1,nTimeSteps); %variable to hold 10000 time points

tiSteps  = size(ti,2);
posZ       = zeros(1,tiSteps); %variable to hold nPosSteps positions allong the z direction
mFinalVect = zeros(tiSteps,2); %variable to hold the final magnetization calculated for each position
mFinalVect2 = zeros(tiSteps,2); % check final vector after rpulse excitation

%Generates the time line for plotting
for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = 133.3*i*dt;                       %Time in seconds
end

%Generate a list of sampling points allong the z direction
for i=1:tiSteps %i starts at 1 go's to 200
    posZ(i)  = (i-100)*10^-4; %Distance from iso center in meters
end

% set T1 and T2
T1 = 850;
T2 = 60;


%Generate the RF Inversion waveform
for i=1:30000 %i starts at 1 go's to 10000
    rfPulse(i) = (sin(pi*i/30000)^2)*sinc(pi*(i-15000)/15000)*10^-5; %B1+ in Tesla  
    gradAmp(i) = 10*10^-3;                       %Tesla per meter
end


figure
plot(time,rfPulse);
%Generate the RF inverse waveform 

% use Hann window funtion 

l = 30000 ;%lenght of the pulse

h = hann(l); % use hann function

hfunc = transpose(h); % transpose to match rfPulse array

rfPulseI = times(hfunc,rfPulse); % create the new hann function rfPulse

figure
plot(time,rfPulseI);

%make pulse 180 inversion
A = sum(rfPulseI);
rfPulse_FA = gamma*A*dt; % equation where flip angle=gamma*sum of rfPulse*dt

degree = rfPulse_FA*360; % degree = 25.5194 %rfPulse2 = rfPulse*pi/rfPulse_FA;

rfPulseI=180/degree*rfPulseI; %add on remaining degree that pulse needs to generate for 180 degrees

save('rfPulseI');

figure
plot(time,rfPulseI);

% time bandwidth factor number of zeros in waveform

fftrfpI=fft(rfPulseI);

figure
plot(fftrfpI, rfPulseI);