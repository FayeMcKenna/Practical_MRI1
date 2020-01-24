%% Exercise 2.1
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

%load sequence
load('sinc_excitation.mat');
nTimeSteps = size(time,2);

%Allocate the memory needed
nPosSteps  = 200;
posZ       = zeros(1,nPosSteps); %variable to hold the positions allong the z direction
mFinal     = zeros(nPosSteps,3); %variable to hold the final magnetization calculated for each position
gradZ      = zeros(1,nPosSteps); % variable to hold gradient 
%Generate a list of sampling points allong the z direction
for i=1:nPosSteps
    posZ(i)  = (i-100)*10^-4; %Distance from iso center in meters
    gradZ(i) =  1*10^-2;
end  

%Question B

y = displaysequence(time,rfPulse,gradAmp);

%Question D:
tic
m = [0,0,1]'; % start large tip angle starting from [0,0,1]':
for j=1:nPosSteps
dB0AtPosJ =  -.0099*.01; % dB0 of j is its pos 1 * its gradient 1 
m = [0,0,1]'; % start large tip angle starting from [0,0,1]'    
    m = largetipangle(dB0AtPosJ,rfPulse(1), m);% 
    for i=2:15000 %(i starts at 2) for every point in the magnetization process
        dB0AtPosJ = posZ(j)*gradAmp(i); % dB0 of j is its pos * its gradient
         m =  largetipangle(dB0AtPosJ,rfPulse(i), m);% total magentization is updated for one time in rfPulse
    end
    mFinal(j,:) = m; % final magnetization
end    
toc
mx = mFinal(:,1); % mx array
my = mFinal(:,2); % my array
mz = mFinal(:,3); % mz array
mxy = mx+1i*my;%m(1)+1i*m(2) mxy array

%Question E: plot the magnitude and phase of the mr signal as a function
%along the z position

%plot one magnitude of mxy 
figure
subplot(3,1,1) % second subplot in column showing rf amp and rf phase
plot(posZ,abs(mxy),'color','r') % plot red
title('magnitude') % title
ylabel('Tesla') % y measure label
xlabel('PosZ'); % x axis label
legend('RF mag');% legend

subplot(3,1,2) % third subplot in column showing gradient amp
plot(posZ,angle(mxy),'color','g') % plot green
title('phase') % title
ylabel('VF') % y measure label
xlabel('PosZ');% x measure label
legend('RF phase');% legend

subplot(3,1,3) % fourth subplot in column showing gradient amp
plot(posZ,abs(mz),'color','b') % plot blue
title('z magnitude') % title
ylabel('Tesla') % y measure label
legend('zmag');% legend
xlabel('PosZ');% x measure label

% check results with pre-simulated data #yes it matches

load('Lab2p1_ref.mat');
figure
%plot one magnitude of mxy 
subplot(2,1,1) % second subplot in column showing rf amp and rf phase
y=abs(Lab2p1_ref);% second plot y is phase of rfPulse
plot(posZ,y,'color','r') % plot red
title('magnitude') % title
ylabel('Tesla') % y measure label
xlabel('PosZ'); % x axis label
legend('RF mag');% legend

subplot(2,1,2) % third subplot in column showing gradient amp
y=angle(Lab2p1_ref);% third plot y is rfamp
plot(posZ,y,'color','g') % plot green
title('phase') % title
ylabel('VF') % y measure label
xlabel('PosZ');% x measure label
legend('RF phase');% legend


%Question F: make rf pulse x20 of original 

rfPulse20 = rfPulse*20;
tic
m = [0,0,1]'; % start large tip angle starting from [0,0,1]':
for j=1:nPosSteps
dB0AtPosJ =  -.0099*.01; % dB0 of j is its pos 1 * its gradient 1 
m = [0,0,1]'; % start large tip angle starting from [0,0,1]'    
    m = largetipangle(dB0AtPosJ,rfPulse20(1), m);% 
    for i=2:15000 %(i starts at 2) for every point in the magnetization process
        dB0AtPosJ = posZ(j)*gradAmp(i); % dB0 of j is its pos * its gradient
         m =  largetipangle(dB0AtPosJ,rfPulse20(i), m);% total magentization is updated for one time in rfPulse
    end
    mFinal(j,:) = m; % final magnetization
end    
toc

mx2 = mFinal(:,1); % mx2 array
my2 = mFinal(:,2); % my2 array
mz2 = mFinal(:,3); % mz2 array
mxy2 = mx2+1i*my2;%m(1)+1i*m(2) mxy2 array

% save as lab2p1f20.mat
save('lab2p1f20.mat','mxy2','mz');

%normalize new results

scaledmxy2 = (mxy2-min(mxy2(:))) ./ (max(mxy2(:)-min(mxy2(:))));
min(scaledmxy2(:)) % the min is 0
max(scaledmxy2(:)) % the max 1

%make one figure comparing two rfPulses
figure
% subplot(4,1,1)%(1st plot is combination)
x1 = posZ; % first x axis is time
y1 = abs(mxy); % first y axis is the phase 
line(x1,y1,'color','r') % put in red
ax1 = gca; %current axes
ax1.XColor = 'r'; %red
ax1.YColor = 'r'; %red
ax1_pos = ax1.Position; % position of first axes in combo plot
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none'); % position 2nd plot in combo plot and in black
x2 = posZ;% second x axis is time
y2 = abs(mxy2); % second y axis is amp of rfPulse
line(x2,y2,'Parent',ax2,'Color','g') % in green
title('RF Pulse and 20x RF Pulse') %title
ylabel('Tesla per meter ') % y label
legend('rf20'); % legend
xlabel('posZ'); % x label

% why are these different? Would the small-tip-angle approximation be able to predict this result?
% These are so different because the amplitude of the RF is proportional to
% the flip angle.  The small tip angle would not be able to predict these
% results because it requires information from all three directions.

%Question G: use tic toc function to measure how long the nested loops take
% Elapsed time is 6.931698 seconds.

