%% Exercise 4.1
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same


dt    = 10^-5; 
gamma = 42.577*10^6;

%load the voxel model
load('PD.mat');
load('T1.mat');
load('T2.mat');

%Allocate the memory needed
nTimeSteps  = 200;
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform
gradAmp     = zeros(3,nTimeSteps); %variable to hold a gradient waveform
time        = zeros(1,nTimeSteps); %variable to hold the time points


xSteps  = size(T1,1);
ySteps  = size(T1,2);
zSteps  = 1;

dX = 4.0e-3;
dY = 4.0e-3;
dZ = 1.0e-4;

% 3D positions in space
pos = zeros(3,xSteps,ySteps,zSteps);
for k=1:xSteps
    for j=1:ySteps
        for i=1:zSteps
            pos(1,k,j,i) = (k-xSteps/2)*dX;
            pos(2,k,j,i) = (j-ySteps/2)*dY;
            pos(3,k,j,i) = (i-zSteps/2)*dZ;
        end
    end
end


%Generates the time line for plotting
for i=1:nTimeSteps 
    time(i)    = i*10^-7;                       %Time in seconds
end

%Generate the excitation waveform
for i=1:100
    rfPulse(i) = (sin(pi*i/100)^2)*sinc(pi*(i-50)/50)*10^-5;
    gradAmp(3,i) = 5*10^-2; %Z gradients in Tesla per meter
end

% evaluate rfPulse
figure
plot(time,rfPulse);

% flip angle 
%Generate the excitation waveform % flip angle 90
A = sum(rfPulse);
rfPulse_FA = gamma*A*dt; % equation where flip angle=gamma*sum of rfPulse*dt

degree = rfPulse_FA*360; % degree = 25.5194 %rfPulse2 = rfPulse*pi/rfPulse_FA;

rfPulseE=90/degree*rfPulse; %add on remaining degree that pulse needs to generate for 180 degrees

% evaluate rfPulseE
figure
plot(time,rfPulseE);

%hann function to rfPulse
l = 200 ;%length of the pulse

h = hann(l); % use hann function

hfunc = transpose(h); % transpose to match rfPulse array

rfPulseE = times(hfunc,rfPulseE); % create the hann function for inversion rfPulse

% evaluate rfPulseE
figure
plot(time,rfPulseE);

%to excite a 5mm thick transverse slice (z-gradient)
gradAmp= 5*10^-2;
%dF = 1/sum(time); % DF is dF = 1/50*10^-5 %1/T0 is 1.5

%look at magnetization profile of sinc pulse
magfft=fft(rfPulse);

figure
plot(time,magfft);

rb= % look at transform mag 0-rb hz shown
dF = rb/2; %

thickS = 2*pi*dF/gamma*gradAmp;

% time bandwidth factor 3, 
dt = 
timeBWf=dT*dF;

% Use 0.5ms for the slice refocusing gradient
%readout pre-phasor (oriented along x)
%Generate the readout gradients
for i=151:(150+size(T1,1))
     %gradAmp(1,i) = ... %X gradients in Tesla per meter
end


%Generate the pre-phasor and slice refocusing gradients
for i=101:150
    %gradAmp(1,i) =  ... %X gradients in Tesla per meter
    %gradAmp(3,i) =  ... %Z gradients in Tesla per meter         
end


%% Question B:


%% Question C:


%% Question D:

% Allocate memory for measured signal
kSpace  = zeros(size(T1,1),size(T1,2));

tic
for k=1:xSteps  
    for j=1:ySteps
        for i=1:zSteps
            
            %dB0 = ... dot product between G(t).r 
            [mT,mZ] =  bloch(dt, dB0,rfPulse(1), T1(k,j),T2(k,j), 0, 1);   % start from fully relaxed spin state
    
            for t=2:nTimeSteps %t starts at 2

                %dB0 = ... dot product between G(t).r 
                %[mT,mZ] =  bloch....
                
                if(adc(t)>0)
                    %Sum the signal over all spins and store in k-space
                    %Don't forget to wigh the signal with the PD.
                    %kSpace(24,round(adc(t))) = kSpace... 
                end
                
            end  %end of time loop          
            
            
        end %end of i loop
    end %end of j loop
end %end of k loop
toc

%% Question D: