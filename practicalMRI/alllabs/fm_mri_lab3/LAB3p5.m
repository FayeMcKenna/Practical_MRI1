%% Exercise 3.5
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same


dt    = 10^-5; 
gamma = 42.577*10^6;

te = [2 3 4 6 8 12 16 24 32 48 64 96 128 192 256]*10^-3; % te in ms



%Allocate the memory needed
nTimeSteps  = 10000;
rfPulseE    = zeros(1,200); %exciation
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform with 10000 samples
gradAmp     = zeros(1,nTimeSteps); %variable to hold a gradient waveform with 10000 samples
time        = zeros(1,nTimeSteps); %variable to hold 10000 time points
time1       = zeros(1,9800);

tiSteps  = size(te,2);
posZ       = zeros(1,tiSteps); %variable to hold nPosSteps positions allong the z direction
mFinalVect = zeros(tiSteps,2); %variable to hold the final magnetization calculated for each position

% time line for plotting
for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = i*dt;                       %Time in seconds
end

% create new posZ with randomness many spins
for i=1:5000 %i starts at 1 go's to 200 *333.3
    posZ(i)  = i*.00000025*rand(); %Distance from iso center in meters
end

for i=5001:10000 %i starts at 1 go's to 200 *333.3
    posZ(i)  = i*-.00000025*rand(); %Distance from iso center in meters
end


% set T1 and T2
T1 = 850;
T2 = 60;


%% Excitation pulse

%  90 , 2ms TBW factor 3 excitation Sinc pulse. 
%Generate the excitation pulse 2ms 3TBW Sinc pulse.
for i=1:200 %i starts at 12000 when other pulse ends
    rfPulseE(i) = (sin(pi*i/200)^2)*sinc(pi*(i-125)/75)*10^-5; %B1+ in Tesla  
                          
end

for i=1:200 %i starts at 1 go's to 15000
    timeE(i)    = i*dt;                       %Time in seconds
end
plot(timeE, rfPulseE) 


%Generate the RF excitation waveform to be 180 degrees
 A = sum(rfPulseE);
rfPulse_FA = gamma*A*dt; % equation where flip angle=gamma*sum of rfPulse*dt

degree = rfPulse_FA*360; % degree = 25.5194 %rfPulse2 = rfPulse*pi/rfPulse_FA;

rfPulsedE=90/degree*rfPulseE; %add on remaining degree that pulse needs to generate for 180 degrees


%% RFpulses and Gradients 

% add together rfpulsedI time1 rdPulsedE time2 time3 to make one rfpulse series;
 
rfPulse = [rfPulsedE,time1]; %whole rfpulse series

% add all gradAmps together during rfPulses and readout

%excitation SS gradAmp
for i=1:200 %i starts at 1 go's to 15000
    gradAmp(i) = 10*10^-3;                       %Time in seconds
end

%readout gradAmp 699401:700000
for i=9400:10000 %i starts at 1 go's to 15000
    gradAmp(i) = 10*10^-3;                       %Time in seconds
end


 y = displaysequence(time,rfPulse,gradAmp);

%% magentization
% function needs input:[ mtOut, mzOut  ] = bloch(dt,dB0,B1,T1,T2,mt,mz) %1sec
for j=1:tiSteps
T1 = 850;
T2 = 60;
gradAtPosJ =  posZ(j)*gradAmp(1); % dB0 of j is its pos 1 * its gradient 1 
    [mT,mZ] = bloch(dt,gradAtPosJ,rfPulse(1),T1,T2,0,1);% 
    for i=2:nTimeSteps %(i starts at 2) for every point in the magnetization process
        gradAtPosJ = posZ(j)*gradAmp(i); % dB0 of j is its pos * its gradient
        [mT,mZ] =  bloch(dt,gradAtPosJ,rfPulse(i),T1,T2,mT,mZ);% total magentization is updated for one time in rfPulse
    end
    mFinalVect(j,:) = [mT, mZ] ; %
end

%% look at results
mz = mFinalVect(:,2); % mz array
mxy = mFinalVect(:,1); % mx array
% figure
% profplot(mxy,mz,posZ); 

%find ?e t/T2 from the equation mxy = xy0(1-exp(-TR/T1)) exp(-TE/T2) 
mxy = mFinalVect(:,1); % mx array
mxy=abs(mxy)
t2=zeros(1,15);
tet2star= transpose(te);

%t2mxy1= tet/mxy1;

for j=mxy
    t2starmxy= exp(-tet2star/mxy); % need to edit this to be one column 
end

t2starmxy= t2starmxy(:,3);

save('t2starmxy.mat')
save('tet2star.mat')

% Plot the measured signal (amplitude and phase) as a function of Te.
figure
profplot(t2starmxy,mz,tet2star); 

%Fit the appropriate signal equation to the measured data points to T2star

% transpose posZ for fitting
% posZ_v =posZ';
%fit data
fitt2star = fit( t2starmxy, tet2star, 'poly2' );
%Plot your fit and the data.
figure
plot( fitt2star, t2starmxy, tet2star )

% it fits very well

%compare both t2 and t2star
load('t2mxy.mat')
load('tet2.mat')

% plot all together and compare vs TE

x=tet2star;
y1=t2mxy;
y2=t2starmxy;

figure
plot(x,y1,x,y2);
title('Comparison of TE T2 T2star') % title
ylabel(' measure') % y measure label
xlabel('TE');% x measure label
legend('t2','t2star');% legend

% load T2 results, how much shorter is T2*? fit same equation
f = fit( t2starmxy, tet2star, 'poly2' );

%Plot your fit and the data.

tet2= transpose(te);
plot( fitt2star, t2mxy, tet2 )
title('T2star fit on T2 data') % title
ylabel(' measure') % y measure label
xlabel('TE');% x measure label
legend('t2','t2star');% legend

%about 9 useconds

