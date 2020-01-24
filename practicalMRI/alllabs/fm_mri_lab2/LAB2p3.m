%% Exercise 2.3
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

%load sequence
load('sinc_excitation.mat');
nTimeSteps = size(time,2);

%Allocate the memory needed
nPosSteps  = 200;
posZ       = zeros(1,nPosSteps); %variable to hold the positions allong the z direction
mFinalVect = zeros(nPosSteps,2); %variable to hold the final magnetization calculated for each position


%Generate a list of sampling points allong the z direction
for i=1:nPosSteps 
    posZ(i)  = (i-100)*10^-4; %Distance from iso center in meters
end  

%Question A:
 y = displaysequence(time,rfPulse,gradAmp);

%Question B:
tic
for j=1:nPosSteps
dB0AtPosJ =  posZ(j)*gradAmp(1); % dB0 of j is its pos 1 * its gradient 1 

    [mT,mZ] = fastsim(dB0AtPosJ,rfPulse(1), 0,1);% 
    for i=2:size(rfPulse,2) %(i starts at 2) for every point in the magnetization process
        dB0AtPosJ = posZ(j)*gradAmp(i); % dB0 of j is its pos * its gradient
        [mT,mZ] =  fastsim(dB0AtPosJ,rfPulse(i), mT,mZ);% total magentization is updated for one time in rfPulse
    end
    mFinalVect(j,:) = [mT, mZ] ; %
end
toc

%%Question C:

mxy = mFinalVect(:,1); % mx array
mz = mFinalVect(:,2); % mz array
figure
profplot(mxy,mz,posZ);

%check results with pre-simulated data 
load('Lab2p1_ref.mat');

mxy = abs(Lab2p1_ref); % mx array
mz = angle(Lab2p1_ref); % mz array
figure
profplot(mxy,mz,posZ);
%plot one magnitude of mxy 
figure
subplot(3,1,1) % second subplot in column showing rf amp and rf phase
x=posZ; % second plot x is time
y=abs(Lab2p1_ref);% second plot y is phase of rfPulse
plot(x,y,'color','r') % plot red
title('magnitude') % title
ylabel('Tesla') % y measure label
xlabel('PosZ'); % x axis label
legend('RF mag');% legend

subplot(3,1,2) % third subplot in column showing gradient amp
x=posZ; % third plot x is time
y=angle(Lab2p1_ref);% third plot y is rfamp
plot(x,y,'color','g') % plot green
title('phase') % title
ylabel('VF') % y measure label
xlabel('PosZ');% x measure label
legend('RF phase');% legend

% %Question D: The code is faster Elapsed time is 0.062418 seconds. The
% slice profile is the same at the presimulated data.

