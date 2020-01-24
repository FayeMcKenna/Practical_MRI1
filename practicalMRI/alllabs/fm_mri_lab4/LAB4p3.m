%% Exercise 4.3
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
nTimeSteps  = 200;%*48;%9600
rfPulse      = [];
% rfPulse1     = zeros(1,nTimeSteps); %variable to hold a RF waveform
gradAmp1     = zeros(3,nTimeSteps); %variable to hold a gradient waveform
gradAmp     = []; %variable to hold a gradient waveform
 %adc         = zeros(2,nTimeSteps); %variable to hold a gradient waveform
time        = zeros(1,nTimeSteps); %variable to hold the time points


xSteps  = size(T1,1);   %Number of simulated "spins" in the x directions 
ySteps  = size(T1,2);   %Number of simulated "spins" in the y directions 
zSteps  = 1;            %Number of simulated "spins" in the z directions 

dX = 4.0e-3;            %Distance between simulated "spins" in the x directions  [meter]
dY = 4.0e-3;            %Distance between simulated "spins" in the y directions  [meter]
dZ = 1.0e-4;            %Distance between simulated "spins" in the z directions  [meter]

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

%Generates the time line for sequence plotting
for i=1:nTimeSteps 
    time(i)    = i*dt;                       %Time in seconds
end

%% first excitation pulse

%Generate the excitation pulse 2ms 3TBW Sinc pulse.
for i=1:100 %i 
    rfPulse1(i) = (sin(pi*i/100)^2)*sinc(pi*(i-35)/60)*10^-5; %B1+ in Tesla          
end

%Generate the RF inverse waveform to be 30 degrees
 A = sum(rfPulse1);
rfPulse_FA = gamma*A*dt; % equation where flip angle=gamma*sum of rfPulse*dt

degree = rfPulse_FA*360; % degree = 25.5194 %rfPulse2 = rfPulse*pi/rfPulse_FA;

rfPulse1=30/degree*rfPulse1; %add on remaining degree that pulse needs to generate for 180 degrees

% use Hann window funtion to apodize excite and invert rfPulses

l = 100 ;%lenght of the pulse

h = hann(l); % use hann function

hfunc = transpose(h); % transpose to match rfPulse array

rfPulse1 = times(hfunc,rfPulse1); % create the hann function for inversion rfPulse

%after pulses 
for i=1:100 %i 
    rfPulse2(i) = 0; %B1+ in Tesla          
end
%inbetween pulses 
for i=1:50 %i 
    rfPulse3(i) = 0; %B1+ in Tesla          
end

rfPulse = [rfPulse1, rfPulse2, rfPulse3,rfPulse1, rfPulse2];

%Generates the time line for sequence plotting
for i=1:450
    timeB(i)    = i*dt;                       %Time in seconds
end


figure; plot(timeB,rfPulse);

%% Gradients
%part1
N = 3;
timePulse = 1*10^-5;
bW = N/timePulse;
%part2
gamma = 42.577*10^6;
dZ = .005; % 5 mm thick
gZa= bW/(gamma *dZ);

%Generate gradient for slice selection 1 ms
for i=1:100 % 1 ms
    gradAmp1(3,i) = gZa; %Z gradients in Tesla per meter               
end

% find gradamp you need for readout from equation 
gradread=2*bW/gamma*48; % pi is out of gamma in.5 ms

% prephase readout same 2xarea of readout to get te at center of readout
gprephase= 2* gradread;

%Generate the pre-phasor and slice refocusing gradients  in x and z, should
% Gs refoc same area of slice select gradient but in .5 ms
gssarea= 1* gZa ;% 
gAmprefoc= gssarea/.5; % is .5 instead of 1


%Generate the readout gradients; 
for i=151:(150+size(T1,1))
     gradAmp1(1,i) = gradread; %X gradients in Tesla per meter
end

% phase encoding gradient maximum
Gp = .5*(48-1)*(1/48); % from equation 

%Generate the readout gradients and add the ADC event
for i=151:(150+size(T1,1))
     gradAmp1(1,i) = gradread; %X gradients in Tesla per meter
     %adc(i) = .... later
end

%Generate the pre-phasor, phase encoding and  slice refocusing gradients
for i=101:150
    gradAmp1(1,i) =  -gprephase; %X gradients in Tesla per meter
    gradAmp1(2,i) =  Gp; %Y gradients in Tesla per meter %% should this be positive and negative?
    gradAmp1(3,i) =  -gAmprefoc; %Z gradients in Tesla per meter         
end

for i=1:50
    gradAmpB(2,i)=0;
    gradAmpB(1,i)=0;
    gradAmpB(3,i)=0;
end

% add all gradAmps together
gradAmp=[gradAmp1, gradAmpB, gradAmp1];
%% Generate the ADC event
adc2= zeros(1,152); % non-phase encoding

for i=(1:48)
    adc1(i)=i;% phase encoding
end

adcB= []; % adding two adc events together

for i=(1:50)
    adcB(1,i)=0; % other non-encoding TR
    adcB(2,i)=0;
end

% adc matrix 
adc3=[adc2,adc1];
adc4 = [adc3;adc3];
adc = [adc4,adcB,adc4];

%%  %plot the compleate sequence

figure
subplot(3,1,1)
plot(timeB,rfPulse)
title('RF pulse') % title
ylabel('%B1+ in Tesla') % y measure label
xlabel('time'); % x axis label

subplot(3,1,2)
plot(timeB,gradAmp)
title('Gradients') % title
ylabel('Tesla per meter') % y measure label
xlabel('time'); % x axis label
legend('Gx','Gz','Gy');% legend

subplot(3,1,3)
plot(timeB,adc)
title('ADC') % title
ylabel('amplitude') % y measure label
xlabel('time'); % x axis label


%% Question B:

kSpace  = zeros(size(T1,1),size(T1,2));

tic

for k=1:xSteps
for t=101:150
gradAmp(2,i) =  t; %Y gradients in Tesla per meter
    disp(k);
    for j=1:ySteps
        for i=1:zSteps
            dB0 = pos(:,k,j)'*gradAmp(:,1); %dB0 = pos(:,k,j)'*gradAmp(:,1);
            [mT,mZ] =  bloch(dt, dB0,rfPulse(1), T1(k,j),T2(k,j), 0, 1);   % start from fully relaxed spin state
            for t=2:nTimeSteps %t starts at 2
                dB0 = pos(:,k,j)'*gradAmp(:,t); %dB0 = pos(:,k,j)'*gradAmp
                [mT,mZ] =  bloch(dt, dB0,rfPulse(t), T1(k,j),T2(k,j), mT, mZ);
                if(adc(1,t)>0)
                    kSpace(round(adc(1,t)),round(adc(2,t))) = kSpace(round(adc(1,t)),round(adc(2,t)))+mT*PD(k,j);
                end  %end of time loop
            end %end of i loop
        end %end of j loop
    end %end of k loop
end
end %end of TR loop
toc


%% Question C:

img = fftshift(ifft2(fftshift(kSpace*sqrt(size(kSpace,1)*size(kSpace,2)))));
figure(1); 
subplot(1,2,1);imagesc(abs(img)); colormap gray; axis equal tight off;title('image')
subplot(1,2,2);imagesc(log(abs(kSpace))); colormap gray; axis equal tight off;title('k-space')


% The artifacts are from transverse magnetization from previous rfPulse