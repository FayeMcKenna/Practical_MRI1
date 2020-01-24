%% Exercise 3.4
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


%Generates the time line for plotting
for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = i*dt;                       %Time in seconds
end

%Generate a list of sampling points allong the z direction
for i=1:tiSteps %i starts at 1 go's to 200
    posZ(i)  = (i-100)*10^-4; %Distance from iso center in meters
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

% for i=1:200 %i starts at 1 go's to 15000
%     timeE(i)    = i*dt;                       %Time in seconds
% end
% plot(timeE, rfPulseE) 


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

%find ?e t/T2 from the equation mxy = xy0(1-exp(-TR/T1)) exp(-TE/T2) 
mxy = mFinalVect(:,1); % mx array
mxy1=abs(mxy)
t2=zeros(1,15);
tet2= transpose(te);

%t2mxy1= tet/mxy1;

for j=mxy
    t2mxy= exp(-tet2/mxy1);
end

save('t2mxy.mat');

t2mxy = t2mxy(:,15);

% Plot the measured signal (amplitude and phase) as a function of T2.

figure
profplot(t2mxy,mz,tet2); 
title('T2') % title
ylabel(' T2') % y measure label
xlabel('TE');% x measure label
legend('t2');% legend

%Fit the appropriate signal equation to the measured data points to T2

% transpose posZ for fitting
posZ_v =posZ';
%fit data
fitt2 = fit(t2mxy, tet2, 'poly2'); %'rat23' 
%Plot your fit and the data.
figure
plot( fitt2, t2mxy, tet2 )
title('T2 FIT') % title
ylabel(' T2') % y measure label
xlabel('TE');% x measure label
legend('t2');% legend

save('fitt2.mat')
save('tet2.mat')


save('tet2')
save('t2mxy')
save('fitt2')

% it fits well

%% run with other te and noise

te = [12, 24, 48]*10^-3;;

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


%Generates the time line for plotting
for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = i*dt;                       %Time in seconds
end

%Generate a list of sampling points allong the z direction
for i=1:tiSteps %i starts at 1 go's to 200
    posZ(i)  = (i-100)*10^-4; %Distance from iso center in meters
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

%% magentization
% function needs input:[ mtOut, mzOut  ] = bloch(dt,dB0,B1,T1,T2,mt,mz) %1sec
for j=1:tiSteps
T1 = 850;
T2 = 60;
gradAtPosJ =  posZ(j)*gradAmp(1); % dB0 of j is its pos 1 * its gradient 1 
    [mT,mZ] = bloch(dt,gradAtPosJ,rfPulse(1),T1,T2,0,1);% 
    for i=2:nTimeSteps %(i starts at 2) for every point in the magnetization process
        gradAtPosJ = posZ(j)*gradAmp(i)*(0.02*rand - 0.01); % dB0 of j is its pos * its gradient
        [mT,mZ] =  bloch(dt,gradAtPosJ,rfPulse(i),T1,T2,mT,mZ);% total magentization is updated for one time in rfPulse
    end
    mFinalVect(j,:) = [mT, mZ] ; %
end

%% look at results
mz = mFinalVect(:,2); % mz array

%find ?e t/T2 from the equation mxy = xy0(1-exp(-TR/T1)) exp(-TE/T2) 
mxy = mFinalVect(:,1); % mx array
mxy=abs(mxy)
t2=zeros(1,15);
tet= transpose(te);

%t2mxy1= tet/mxy1;

for j=mxy
    t2mxy= exp(-tet/mxy); 
end

t2mxy = t2mxy(:,3);

% Plot the measured signal (amplitude and phase) as a function of T2.

figure
profplot(t2mxy,mz,tet); 
title('T2 with noise and small TE') % title
ylabel(' T2') % y measure label
xlabel('TE');% x measure label
legend('t2');% legend

%Fit the appropriate signal equation to the measured data points to T2

% transpose posZ for fitting
posZ_v =posZ';
%fit data
f = fit(t2mxy, tet, 'poly2' );
%Plot your fit and the data.
figure
plot( f, t2mxy, tet )
title('T2 with noise and small TE FIT') % title
ylabel(' T2') % y measure label
xlabel('TE');% x measure label
legend('t2');% legend

% T2 fits worse with noise but line is simple


