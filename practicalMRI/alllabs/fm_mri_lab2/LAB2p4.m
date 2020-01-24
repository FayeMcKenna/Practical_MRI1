%% Exercise 2.4
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

dt    = 10^-7; 
gamma = 42.577*10^6;

%Allocate the memory needed
nTimeSteps  = 30000;
rfPulse     = zeros(1,nTimeSteps); %variable to hold a RF waveform with 30000 samples
gradAmp     = zeros(1,nTimeSteps); %variable to hold a gradient waveform with 30000 samples
time        = zeros(1,nTimeSteps); %variable to hold 30000 time points

nPosSteps  = 200;
posZ       = zeros(1,nPosSteps); %variable to hold 400 positions allong the z direction
mFinalVect = zeros(nPosSteps,2); %variable to hold the final magnetization calculated for each position
mFinal = zeros(1,nPosSteps);

%Generate a list of sampling points allong the z direction
for i=1:nPosSteps %i starts at 1 go's to 200
    posZ(i)  = (i-100)*10^-4; %Distance from iso center in meters
end  

%1 make a 2ms inversion pulse with an approximate slice thickness of 1cm.
for i=1:20000 %i starts at 1 go's to 10000
    rfPulse(i) = ((sin(pi*i/20000)^2)*sinc(pi*(i-10000)/10000)*10^-5); %B1+ in Tesla ((sin(pi*i/20000)^2)*sinc(pi*(i-10000)/5000)*10^-5)
end

%Generates the time line for plotting
for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = i*10^-4;                       %Time in seconds *.66
end

%create inversion pulse
A = sum(rfPulse);
rfPulse_FA = gamma*A*dt; % equation where flip angle=gamma*sum of rfPulse*dt

degree = rfPulse_FA*360; % degree = 25.5194 %rfPulse2 = rfPulse*pi/rfPulse_FA;

rfPulse=180/degree*rfPulse; %add on remaining degree that pulse needs to generate for 180 degrees

% use Hann window funtion 

l = 30000 ;%lenght of the pulse

h = hann(l); % use hann function

hfunc = transpose(h); % transpose to match rfPulse array

rfPulse = times(hfunc,rfPulse); % create the new hann function rfPulse

%check inversion rf Pulse
figure
plot(time,rfPulse);

% Create slice select gradient 

% slice select gradient to excite 1 cm thick
%part1
%N = time * BW is the TBW product equation
N = 4;
timePulse = 10^-4;
bW = N/timePulse;
%part2
gamma = 42.577*10^6;
dZ = .01;
gZa= bW/(gamma *dZ);

%slice selection gradient 
for i=1:20000 %i starts at 1 go's to 10000 
    gradAmp(i) = gZa;%Tesla per meter
end

% add a slice re-focusing gradient to the wave form 1 cm thick

for i=20001:30000 %i starts at 1 go's to 10000
    gradAmp(i) = -gZa; % applying opposite
end

% Profile of sequence
y = displaysequence(time,rfPulse,gradAmp);

% run with simulation
for j=1:nPosSteps
dB0AtPosJ =  posZ(j)*gradAmp(1); % dB0 of j is its pos 1 * its gradient 1 

    [mT,mZ] = fastsim(dB0AtPosJ,rfPulse(1), 0,1);% 
    for i=2:size(rfPulse,2) %(i starts at 2) for every point in the magnetization process
        dB0AtPosJ = posZ(j)*gradAmp(i); % dB0 of j is its pos * its gradient
        [mT,mZ] =  fastsim(dB0AtPosJ,rfPulse(i), mT,mZ);% total magentization is updated for one time in rfPulse
    end
    mFinalVect(j,:) = [mT, mZ] ; %
end

%plotslice profile of rfPulse with fastsim
mxy = mFinalVect(:,1); % mx array
mz = mFinalVect(:,2); % mz array
figure
profplot(mxy,mz,posZ); 

%Question A: Explain in a comment why you think this is (close to) the best possible solution.
% Even though it is impossible to get all of the areas in the slice to be
% stimulated perfectly by the inversion pulse, taking the sum of the range of FA of the slice this 
%is the best estimate. 

%Question B:What is happening outside the slice?
%Do you think this will be a problem? Explain why this may/may-not be a problem.

% There is some magnetization outside of the slice this could be a problem
% if the change in spin/phase affects spins/phase inside the selected slice on the
% edge.

%Question C: Do you think it is possible to completely invert all the spins in the slice? Explain why you think this is the case.

% No, I do not think it is possible to completely invert all the spins in
% the slice because there is naturally a range of frequencies.

%Question D: Simulate the same pulse using the small-tip-angle approximation (previous chapter).

% simulate rfPulse2 with smalltipangle function
dB0AtPosJ  = 0 ; % klunge

for j=1:200 %loop over different positions allong the z direction
    m =  smalltipangle(dB0AtPosJ,rfPulse(1), 0); % first time point assuming the transverse magnetization 0
    for i=2:30000 %(i starts at 2) for every point in the magnetization process
        dB0AtPosJ = posZ(j)*gradAmp(i); % dB0 of j is its pos * its gradient
         m =  smalltipangle(dB0AtPosJ,rfPulse(i), m);% total magentization is updated for one time in rfPulse
    end
    mFinal(j) = m; %magnetization calculated for each j point
end    

%plot rfPulse2 with small tipangle mFinal slice profile

figure
subplot(2,1,1)
plot( posZ, angle(mFinal))
ylabel('VF') % y measure label
title('phase');% legend
xlabel('position');% x measure label

subplot(2,1,2)
plot( posZ, abs(mFinal))
ylabel('Tesla') 
title('magnitude')
xlabel('position');% x measure label


%What are the differences in the slice profile compared to the large-tip-angle simulation?
%Is this even physically possible? Explain your observations in a comment.

% The mz is ignored in the small tip angle approximation so you do not see
% the inversion of the pulse, but only the magnetization as you would see
% in the mx direction.  Although there must be an inversion, it is not detected. 