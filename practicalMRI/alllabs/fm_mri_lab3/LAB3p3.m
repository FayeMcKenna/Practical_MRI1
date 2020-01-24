
%% Exercise 3.3
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

dt    = 10^-5; 
gamma = 42.577*10^6;

ti = [25 50 100 200 400 800 1600 3200 6400]*10^-3; % ti in ms

%Allocate the memory needed
nTimeSteps  = 700000;
rfPulseI    = zeros(1,400); %inversion
rfPulseE    = zeros(1,200); %excitation
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform with 10000 samples
gradAmp     = zeros(1,nTimeSteps); %variable to hold a gradient waveform with 10000 samples
time        = zeros(1,nTimeSteps); %variable to hold 10000 time points
time1       = zeros(1,465870); % time between rf pulses
time2       = zeros(1,232930); %time between excitation pulse and readout
time3       = zeros(1,600); %time of readout
tiSteps  = size(ti,2);
posZ       = zeros(1,tiSteps); %variable to hold nPosSteps positions allong the z direction
mFinalVect = zeros(tiSteps,2); %variable to hold the final magnetization calculated for each position
%mFinalVect2 = zeros(tiSteps,2); % check final vector after rpulse excitation

%Generates the time line for plotting
for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = 133*i*dt;                       %Time in seconds
end

%Generate a list of sampling points allong the z direction
for i=1:tiSteps %i starts at 1 go's to 200
    posZ(i)  = (i-100)*10^-4; %Distance from iso center in meters
end

% set T1 and T2
T1 = 850;
T2 = 60;

%% Inversion pulse

%generate inversion pulse 4ms TBW factor 5 apodized Sinc pulse
for i=1:400 %i starts at 1 go's to 10000
    rfPulseI(i) = (sin(pi*i/300)^2)*sinc(pi*(i-100)/100)*10^-5; %B1+ in Tesla  50
end

%timeline to check rfpulseE
for i=1:400 %i starts at 1 go's to 15000
    timeI(i)    = i*dt;                       %Time in seconds
end

%Generate the RF inverse waveform to be 180 degrees
 A = sum(rfPulseI);
rfPulse_FA = gamma*A*dt; % equation where flip angle=gamma*sum of rfPulse*dt

degree = rfPulse_FA*360; % degree = 25.5194 %rfPulse2 = rfPulse*pi/rfPulse_FA;

rfPulsedI=180/degree*rfPulseI; %add on remaining degree that pulse needs to generate for 180 degrees

% use Hann window funtion to apodize excite and invert rfPulses

l = 400 ;%lenght of the pulse

h = hann(l); % use hann function

hfunc = transpose(h); % transpose to match rfPulse array

rfPulsehI = times(hfunc,rfPulsedI); % create the hann function for inversion rfPulse

%check inversion pulse

save('rfPulsehI')

%% Excitation pulse

%Generate the excitation pulse 2ms 3TBW Sinc pulse.
for i=1:200 %i starts at 12000 when other pulse ends
    rfPulseE(i) = (sin(pi*i/200)^2)*sinc(pi*(i-125)/100)*10^-5; %B1+ in Tesla  
                          %Tesla per meter
end

%timeline to check rfpulseE
for i=1:200 %i starts at 1 go's to 15000
    timeE(i)    = i*dt;                       %Time in seconds
end

%Generate the RF excitation waveform to be 180 degrees
 A = sum(rfPulseE);
rfPulse_FA = gamma*A*dt; % equation where flip angle=gamma*sum of rfPulse*dt

degree = rfPulse_FA*360; % degree = 25.5194 %rfPulse2 = rfPulse*pi/rfPulse_FA;

rfPulsedE=45/degree*rfPulseE; %add on remaining degree that pulse needs to generate for 180 degrees

%check excitation pulse

save('rfPulsedE')

%% add pulses and gradients together

% add together rfpulsedI time1 rdPulsedE time2 time3 to make one rfpulse series;
 
rfPulse = [rfPulsehI,time1,rfPulsedE,time2,time3]; %whole rfpulse series

% add all gradAmps together during rfPulses and readout
%inversion SS gradAmp
for i=1:400 %i starts at 1 go's to 15000
    gradAmp(i) = 10*10^-3;                       %Time in seconds
end

%excitation SS gradAmp
for i=466271:466470 %i starts at 1 go's to 15000
    gradAmp(i) = 10*10^-3;                       %Time in seconds
end

%readout gradAmp 699401:700000
for i=699401:700000 %i starts at 1 go's to 15000
    gradAmp(i) = 10*10^-3;                       %Time in seconds
end


%% Get magnetization

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

% Plot the measured signal (amplitude and phase) as a function of TI.
mxy = mFinalVect(:,1); % mx array
mz = mFinalVect(:,2); % mz array
figure
profplot(mxy,mz,posZ); 

%Fit the appropriate signal equation to the measured data points to estimate the T1
%(Matlab has build in functions for fitting that you can use). 
%How well does this correspond to the ground truth?

%T1 is longitudinal Mz magnetization 
absmz= transpose(mz);
absmz= abs(absmz)';
% need to transpose both vectors to fit
posZ_v =posZ';
%fit data
f = fit( absmz, posZ_v, 'rat23' );
%Plot your fit and the data.
figure
plot( f, absmz, posZ_v )

% it corresponds well

%%  What if you only take TI = {200, 800, 3200} ms and add some noise to each data point (0.02 ? rand   0.01)? 
%Repeat this several times, and record the values. How accurate are your T1 estimations now?

ti = [200, 800, 3200]*10^-3; % ti in ms

%% Get magnetization


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

% Plot the measured signal (amplitude and phase) as a function of TI.
mxy = mFinalVect(:,1); % mx array
mz = mFinalVect(:,2); % mz array
figure
profplot(mxy,mz,posZ); 

%Fit the appropriate signal equation to the measured data points to estimate the T1
%(Matlab has build in functions for fitting that you can use). 
%How well does this correspond to the ground truth?

%T1 is longitudinal Mz magnetization 
absmz= abs(mz);
% need to transpose both vectors to fit
absmz= transpose(absmz)';
posZ_v =posZ';
%fit data
f = fit( absmz, posZ_v, 'rat23' );
%Plot your fit and the data.
figure
plot( f, absmz, posZ_v )

%T1 estimates are worse with noise