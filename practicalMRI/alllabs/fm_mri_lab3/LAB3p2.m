%% Exercise 3.1
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same


dt    = 10^-6; 
gamma = 42.577*10^6;


%Allocate the memory needed
nTimeSteps  = 30000;
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform with 10000 samples
gradAmp     = zeros(1,nTimeSteps); %variable to hold a gradient waveform with 10000 samples
time        = zeros(1,nTimeSteps); %variable to hold 10000 time points

nPosSteps  = 200;
posZ       = zeros(1,nPosSteps); %variable to hold nPosSteps positions allong the z direction
mFinalVect = zeros(nPosSteps,2); %variable to hold the final magnetization calculated for each position


%Generates the time line for plotting
for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = i*10^-7;                       %Time in seconds
end

%Generate a list of sampling points allong the z direction
for i=1:nPosSteps %i starts at 1 go's to 200
    posZ(i)  = (i-100)*10^-4; %Distance from iso center in meters
end  

%Generate the RF Inversion waveform
for i=1:20000 %i starts at 1 go's to 10000
    rfPulse(i) = (sin(pi*i/20000)^2)*sinc(pi*(i-10000)/19000)*10^-5; %B1+ in Tesla  
    gradAmp(i) = 10*10^-3;                       %Tesla per meter
end

for i=20001:30000 %i starts at 1 go's to 10000
     gradAmp(i) = -gradAmp(1000);  
end

%Question 0:
%plot sequence the sequce

 y = displaysequence(time,rfPulse,gradAmp);


%Question B:

%function [ mtOut, mzOut  ] = bloch(dt,dB0,B1,T1,T2,mt,mz) %1sec
tic
for j=1:nPosSteps
T1 = 2;
T2 = 2;
gradAtPosJ =  posZ(j)*gradAmp(1); % dB0 of j is its pos 1 * its gradient 1 
    [mT,mZ] = bloch(dt,gradAtPosJ,rfPulse(1),T1,T2,0,1);% 
    for i=2:nTimeSteps %(i starts at 2) for every point in the magnetization process
        gradAtPosJ = posZ(j)*gradAmp(i); % dB0 of j is its pos * its gradient
        [mT,mZ] =  bloch(dt,gradAtPosJ,rfPulse(i),T1,T2,mT,mZ);% total magentization is updated for one time in rfPulse
    end
    mFinalVect(j,:) = [mT, mZ] ; %
end
toc


%Question B

%plotslice profile of rfPulse
mxy = mFinalVect(:,1); % mx array
mz = mFinalVect(:,2); % mz array
figure
profplot(mxy,mz,posZ); 

% compare profile to Lab3p1 ref


%Question C

%Elapsed time is 3.222726 seconds compared to 2.573273 seconds ( compared
%to 0.062418 seconds fastsim).  The slice profile is still somewhat accurate, but the mxy magnitude shows a split profile.


