%% Exercise 2.4
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same


dt    = 10^-7; 
gamma = 42.577*10^6;


%Allocate the memory needed
nTimeSteps   = 30000;
rfPulse      = zeros(1,nTimeSteps); %variable to hold a RF waveform with 30000 samples
rfPulse45    = zeros(1,nTimeSteps); %variable to hold a RF waveform with 30000 samples
rfPulse180   = zeros(1,nTimeSteps); %variable to hold a RF waveform with 30000 samples
rfPulseA     = zeros(1,60000);% add all together 
rfPulseN     = zeros(1,20000);% nonselective pulse
rfPulseB     = zeros(1,20000);% time between pulses
gradAmp      = zeros(1,nTimeSteps); %variable to hold a gradient waveform with 30000 samples
gradAmp2     = zeros(1,nTimeSteps); %variable to hold a gradient waveform with 30000 samples
gradAmpB     =zeros(1,20000);% time between gradients
gradAmpA     = zeros(1,60000); % add all together 
gradAmp2x    = zeros(1,30000); % 2x inversion gradient
time         = zeros(1,nTimeSteps); %variable to hold 30000 time points
timeA        = zeros(1,80000);% add all together 
timeN        = zeros(1,20000) % time of non-selective rfPulse

nPosSteps  = 200;
posZ       = zeros(1,nPosSteps); %variable to hold 400 positions allong the z direction
mFinalVect = zeros(nPosSteps,2); %variable to hold the final magnetization calculated for each position
%mFinalVect45 = zeros(nPosSteps,2); %variable to hold the final magnetization calculated for each position

%Generates the time line for plotting
for i=1:nTimeSteps %i starts at 1 go's to 15000
    time(i)    = i*10^-4;                       %Time in seconds
end

%Generate a list of sampling points allong the z direction
for i=1:nPosSteps %i starts at 1 go's to 200
    posZ(i)  = (i-100)*10^-4; %Distance from iso center in meters
end  

%% Question A: 

%1)Create a 180 inversion (2ms pulse duration, time bandwidth factor 3, apodized Sinc, 1.0cm slice thickness)
for i=1:20000 %i starts at 1 go's to 10000
    rfPulse(i) = ((sin(pi*i/20000)^2)*sinc(pi*(i-12000)/9000)*10^-5); %B1+ in Tesla ((sin(pi*i/20000)^2)*sinc(pi*(i-10000)/5000)*10^-5)
end

%create the 180 inversion 
A = sum(rfPulse);
rfPulse_FA = gamma*A*dt; % equation where flip angle=gamma*sum of rfPulse*dt

degree = rfPulse_FA*360; % degree = 25.5194 %rfPulse2 = rfPulse*pi/rfPulse_FA;

rfPulse=180/degree*rfPulse; %add on remaining degree that pulse needs to generate for 180 degrees

% use Hann window function 

l = 30000 ;%lenght of the pulse

h = hann(l); % use hann function

hfunc = transpose(h); % transpose to match rfPulse array

rfPulse180 = times(hfunc,rfPulse); % create the new hann function rfPulse

% slice select gradient to excite 1 cm thick for rfPulse180 %N = time * BW is the TBW product equation
%part1
N = 3;
timePulse = 10^-4;
bW = N/timePulse;
%part2
gamma = 42.577*10^6;
dZ = .01;
gZa= bW/(gamma *dZ);

%slice selection gradient for rfPulse180
for i=1:20000 %i starts at 1 go's to 10000 
    gradAmp(i) = gZa;%Tesla per meter
end

% add a slice re-focusing gradient for rfPulse180

for i=20001:30000 %i starts at 1 go's to 10000
    gradAmp(i) = -gZa; % applying opposite
end

%%B) Create a 45 excitation pulse (2ms pulse duration, time bandwidth factor, 2, 
%apodized Sinc, 0.5cm slice thickness), with a 2ms delay between them.

for i=1:20000 %i starts at 1 go's to 10000
    rfPulse(i) = ((sin(pi*i/20000)^2)*sinc(pi*(i-10000)/15000)*10^-5); %B1+ in Tesla ((sin(pi*i/20000)^2)*sinc(pi*(i-10000)/5000)*10^-5)
end

%create the 45 inversion 
A = sum(rfPulse);
rfPulse_FA = gamma*A*dt; % equation where flip angle=gamma*sum of rfPulse*dt

degree = rfPulse_FA*360; % degree = 25.5194 %rfPulse2 = rfPulse*pi/rfPulse_FA;

rfPulse=45/degree*rfPulse; %add on remaining degree that pulse needs to generate for 180 degrees

% use Hann window function 

l = 30000 ;%lenght of the pulse

h = hann(l); % use hann function

hfunc = transpose(h); % transpose to match rfPulse array

rfPulse45 = times(hfunc,rfPulse); % create the new hann function rfPulse

% slice select gradient to excite .5 cm thick for rfPulse45 %N = time * BW is the TBW product equation
%part1
N = 2;
timePulse = 10^-4;
bW = N/timePulse;
%part2
gamma = 42.577*10^6;
dZ = .005;
gZa2= bW/(gamma *dZ);

%slice selection gradient for rfPulse180
for i=1:20000 %i starts at 1 go's to 10000 
    gradAmp2(i) = gZa2;%Tesla per meter
end

% add a slice re-focusing gradient for rfPulse180
for i=20001:30000 %i starts at 1 go's to 10000
    gradAmp2(i) = -gZa2; % applying opposite
end

% check profile
y = displaysequence(time,rfPulse180,gradAmp);

z = displaysequence(time,rfPulse45,gradAmp2);

%% Add sequences together 

%Generates the time for all pulses
for i=1:80000 %i starts at 1 go's to 15000
    timeA(i)    = i*10^-4;                       %Time in seconds
end

for i=1:20000 %i starts at 1 go's to 10000
    rfPulseB(i) = 0; %B1+ in Tesla ((sin(pi*i/20000)^2)*sinc(pi*(i-10000)/5000)*10^-5)
end

for i=1:20000 %i starts at 1 go's to 10000
    gradAmpB(i) = 0; %B1+ in Tesla ((sin(pi*i/20000)^2)*sinc(pi*(i-10000)/5000)*10^-5)
end

% create whole rfPulse sequence and whole gradAmp sequence
rfPulseA = [rfPulse180,rfPulseB,rfPulse45];
gradAmpA = [ gradAmp, gradAmpB,gradAmp2];

%display the whole sequence
x = displaysequence(timeA,rfPulseA,gradAmpA);

%% run whole pulse and grad amp sequence

for j=1:nPosSteps
dB0AtPosJ =  posZ(j)*gradAmpA(1); % dB0 of j is its pos 1 * its gradient 1 
    [mT,mZ] = fastsim(dB0AtPosJ,rfPulseA(1), 0,1);% 
    for i=2:size(rfPulse180,2) %(i starts at 2) for every point in the magnetization process
        dB0AtPosJ = posZ(j)*gradAmpA(i); % dB0 of j is its pos * its gradient
        [mT,mZ] =  fastsim(dB0AtPosJ,rfPulseA(i), mT,mZ);% total magentization is updated for one time in rfPulse
    end
    mFinalVect(j,:) = [mT, mZ] ; %
end

mxy = mFinalVect(:,1); % mx array
mz = mFinalVect(:,2); % mz array

%check  pulse slice profile
figure
sliceProf=profplot(mxy,mz,posZ); 
save('sliceProf')

%% 3)Design the gradients such that the transverse magnetization outside of the excitation slice is de-phased as 
%much as possible (rapid phase oscillations).

% create spoiler gradient!
for i=1:10000 %i starts at 1 go's to 10000
    gradAmpB1(i) = gZa; %B1+ in Tesla ((sin(pi*i/20000)^2)*sinc(pi*(i-10000)/5000)*10^-5)
end

for i=1:10000 %i starts at 1 go's to 10000
    gradAmpB2(i) = -gZa; %B1+ in Tesla ((sin(pi*i/20000)^2)*sinc(pi*(i-10000)/5000)*10^-5)
end

% add spoiler gradient to other slice select gradients
gradAmpC = [ gradAmp, gradAmpB1,gradAmpB2,gradAmp2];

%display the whole sequence to check new gradients
x = displaysequence(timeA,rfPulseA,gradAmpC);

%4)Simulate the slice profile from -1 to 1cm in steps of 0.01cm.

for j=1:nPosSteps
dB0AtPosJ =  posZ(j)*gradAmpC(1); % dB0 of j is its pos 1 * its gradient 1 
    [mT,mZ] = fastsim(dB0AtPosJ,rfPulseA(1), 0,1);% 
    for i=2:size(rfPulse180,2) %(i starts at 2) for every point in the magnetization process
        dB0AtPosJ = posZ(j)*gradAmpC(i); % dB0 of j is its pos * its gradient
        [mT,mZ] =  fastsim(dB0AtPosJ,rfPulseA(i), mT,mZ);% total magentization is updated for one time in rfPulse
    end
    mFinalVect2(j,:) = [mT, mZ] ; %
end

mxy2 = mFinalVect2(:,1); % mx array
mz2 = mFinalVect2(:,2); % mz array

%check  pulse slice profile
figure
sliceProfspoil=profplot(mxy2,mz2,posZ); 

save('sliceProfspoil');

%5) Save the results (slice profile) to file.


%Question B: Were you able to perfectly invert all the spins inside the imaging slice? Explain why this may or may-not be a problem.

% I was still unable to perfectly invert all the spins inside the imaging
% slice, this isn't necessarily a problem as long as you take into account
% that all the spins cannot be uniform in a slice in your model.

%Question C:

%Do the results improve if the inversion slice is twice as thick (2cm)?
N = 3;
timePulse = 10^-4;
bW = N/timePulse;
%part2
gamma = 42.577*10^6;
dZ = .02; % twice as thick
gZa= bW/(gamma *dZ);

%slice selection gradient for rfPulse180
for i=1:20000 %i starts at 1 go's to 10000 
    gradAmp2x(i) = gZa;%Tesla per meter
end

% add a slice re-focusing gradient for rfPulse180
for i=20001:30000 %i starts at 1 go's to 10000
    gradAmp2x(i) = -gZa; % applying opposite
end

% create whole rfPulse sequence and whole gradAmp sequence
rfPulseA = [rfPulse180,rfPulseB,rfPulse45];
gradAmpA = [ gradAmp2x, gradAmpB,gradAmp2];

%display the whole sequence
x = displaysequence(timeA,rfPulseA,gradAmpA);

for j=1:nPosSteps
dB0AtPosJ =  posZ(j)*gradAmpA(1); % dB0 of j is its pos 1 * its gradient 1 
    [mT,mZ] = fastsim(dB0AtPosJ,rfPulseA(1), 0,1);% 
    for i=2:size(rfPulse180,2) %(i starts at 2) for every point in the magnetization process
        dB0AtPosJ = posZ(j)*gradAmpA(i); % dB0 of j is its pos * its gradient
        [mT,mZ] =  fastsim(dB0AtPosJ,rfPulseA(i), mT,mZ);% total magentization is updated for one time in rfPulse
    end
    mFinalVect3(j,:) = [mT, mZ] ; %
end

mxy3 = mFinalVect3(:,1); % mx array
mz3 = mFinalVect3(:,2); % mz array

%check  pulse slice profile
figure
sliceProf2x=profplot(mxy3,mz3,posZ); 
save('sliceProf2x');

%What if we use a non-selective pulse? 

% add non-selective pulse

for i=1:20000 %i starts at 1 go's to 10000
    rfPulseN(i) = ((sin(pi*i/10000)^2)*sinc(pi*(i-5000)/5000)*10^-5); %B1+ in Tesla ((sin(pi*i/20000)^2)*sinc(pi*(i-10000)/5000)*10^-5)
end

%Generates the time for all pulses
for i=1:20000 %i starts at 1 go's to 15000
    timeN(i)    = i*10^-4;                       %Time in seconds
end

l = 20000 ;%lenght of the pulse

h = hann(l); % use hann function

hfunc = transpose(h); % transpose to match rfPulse array

rfPulseN = times(hfunc,rfPulseN); % create the new hann function rfPulse

%check nonselective rfPulse
% figure
% plot(timeN,rfPulseN);

% add together new sequence of rfPulse and gradAMp
rfPulseA = [rfPulse180,rfPulseN,rfPulse45];
gradAmpA = [ gradAmp, gradAmpB,gradAmp2];

p = displaysequence(timeA,rfPulseA,gradAmpA);

% find slice profile
for j=1:nPosSteps
dB0AtPosJ =  posZ(j)*gradAmpA(1); % dB0 of j is its pos 1 * its gradient 1 
    [mT,mZ] = fastsim(dB0AtPosJ,rfPulseA(1), 0,1);% 
    for i=2:size(rfPulse180,2) %(i starts at 2) for every point in the magnetization process
        dB0AtPosJ = posZ(j)*gradAmpA(i); % dB0 of j is its pos * its gradient
        [mT,mZ] =  fastsim(dB0AtPosJ,rfPulseA(i), mT,mZ);% total magentization is updated for one time in rfPulse
    end
    mFinalVect4(j,:) = [mT, mZ] ; %
end

mxy4 = mFinalVect4(:,1); % mx array
mz4 = mFinalVect4(:,2); % mz array

%check  pulse slice profile
figure
sliceProfNS=profplot(mxy4,mz4,posZ); 
save('sliceProfNS');

%Plot the results (abs(mt), arg(mt), and mz) produce by all three RF pulse configurations in one single figure. 
%Make sure your figure shows the correct dimensions on each axis and contains a legend.

figure
subplot(9,1,1) % second subplot in column showing rf amp and rf phase
plot(posZ,abs(mxy2),'color','r') % plot red
title('Slice Selection with Gradient Spoiling MXY magnitude') % title
ylabel('Tesla') % y measure label
xlabel('PosZ'); % x axis label
legend('mxy mag');% legend

subplot(9,1,2) % third subplot in column showing gradient amp
plot(posZ,abs(mxy2),'color','g') % plot green
title('Slice Selection with Gradient Spoiling MXY phase') % title
ylabel('VF') % y measure label
xlabel('PosZ');% x measure label
legend('mxy phase');% legend

subplot(9,1,3) % fourth subplot in column showing gradient amp
plot(posZ,mz2,'color','b') % plot blue
title('Slice Selection with Gradient Spoiling mz magnitude') % title
ylabel('Tesla') % y measure label
legend('mz mag');% legend
xlabel('PosZ');% x measure labe

subplot(9,1,4) % second subplot in column showing rf amp and rf phase
plot(posZ,abs(mxy3),'color','r') % plot red
title('Slice Selection with 2x Inversion slice thickness MXY magnitude') % title
ylabel('Tesla') % y measure label
xlabel('PosZ'); % x axis label
legend('mxy mag');% legend

subplot(9,1,5) % third subplot in column showing gradient amp
plot(posZ,abs(mxy3),'color','g') % plot green
title('Slice Selection with 2x Inversion slice thickness MXY phase') % title
ylabel('VF') % y measure label
xlabel('PosZ');% x measure label
legend('mxy phase');% legend

subplot(9,1,6) % fourth subplot in column showing gradient amp
plot(posZ,mz3,'color','b') % plot blue
title('Slice Selection with 2x Inversion slice thickness mz magnitude') % title
ylabel('Tesla') % y measure label
legend('mz mag');% legend
xlabel('PosZ');% x measure labe

subplot(9,1,7) % second subplot in column showing rf amp and rf phase
plot(posZ,abs(mxy4),'color','r') % plot red
title('Slice Selection with Non-selective RF Pulse MXY Amagnitude') % title
ylabel('Tesla') % y measure label
xlabel('PosZ'); % x axis label
legend('mxy mag');% legend

subplot(9,1,8) % third subplot in column showing gradient amp
plot(posZ,abs(mxy4),'color','g') % plot green
title('Slice Selection with Non-selective RF Pulse MXY phase') % title
ylabel('VF') % y measure label
xlabel('PosZ');% x measure label
legend('mxy phase');% legend

subplot(9,1,9) % fourth subplot in column showing gradient amp
plot(posZ,mz4,'color','b') % plot blue
title('Slice Selection with Non-selective RF Pulse mz magnitude') % title
ylabel('Tesla') % y measure label
legend('mz mag');% legend
xlabel('PosZ');% x measure labe

%Question D:  When using the 2cm thick inversion, calculate the total signal from spins outside of the imaging slice 
%(sum of all mt for positions outside the imaging slice) and compare this to the signal from the imaging slice. 
%What fraction of the total signal originates from the imaging slice?

% .9769  comes from the slice

% find what PosZ is in slice 
 plot(posZ, abs(mxy3)) %.002 to -.002 is 80-120

sliceMxy=abs(sum(mxy3(80:120)));
oslice1=abs(sum(mxy(1:79)));
oslice2=abs(sum(mxy(121:200)));
osliceMxy = oslice1 + oslice2;

tot= abs(sum(mxy3)) ;

 fracSig= sliceMxy/tot ; % .9769 

%Question E:Suppose we had much stronger gradients (100mT/m), 
%would that help us eliminate the signal form outside of the slice?

%Yes, but the gradients can only turn on and off at a maximal rate
%determined by the hardware.

%Question F: In the absence of stronger gradients one could make the gradients longer.
%How long would the minimum delay between the inversion and excitation need to be in order 
%to make sure that the contamination from out of slice signal is less than 1%?

% TR is what determins the minimum delay between rfPulses, Scan time =
% TR·Ny·NSA.  

