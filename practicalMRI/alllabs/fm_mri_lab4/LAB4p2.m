
%% 4.2
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
adc         = zeros(1,nTimeSteps); %variable to hold a gradient waveform
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


%Generates the time line for plotting
for i=1:nTimeSteps 
    time(i)    = i*10^-7;                       %Time in seconds
end


%% excitation pulse

%Generate the excitation pulse 2ms 3TBW Sinc pulse.
for i=1:100 %i 
    rfPulse(i) = (sin(pi*i/100)^2)*sinc(pi*(i-35)/60)*10^-5; %B1+ in Tesla          
end


%check tbw
figure
plot(time,rfPulse);

%Generate the RF inverse waveform to be 90 degrees
 A = sum(rfPulse);
rfPulse_FA = gamma*A*dt; % equation where flip angle=gamma*sum of rfPulse*dt

degree = rfPulse_FA*360; % degree = 25.5194 %rfPulse2 = rfPulse*pi/rfPulse_FA;

rfPulse=90/degree*rfPulse; %add on remaining degree that pulse needs to generate for 180 degrees

% use Hann window funtion to apodize excite and invert rfPulses

l = 100 ;%lenght of the pulse

h = hann(l); % use hann function

hfunc = transpose(h); % transpose to match rfPulse array

rfPulse = times(hfunc,rfPulse); % create the hann function for inversion rfPulse

%check hann rfpulse
figure
plot(time,rfPulse);

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
    gradAmp(3,i) = gZa; %Z gradients in Tesla per meter               
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
     gradAmp(1,i) = gradread; %X gradients in Tesla per meter
end

% phase encoding gradient maximum
Gp = .5*(48-1)*(1/48); % .4896 from equation 

%Generate the readout gradients and add the ADC event
for i=151:(150+size(T1,1))
     gradAmp(1,i) = gradread; %X gradients in Tesla per meter
     %adc(i) = .... later
end

%Generate the pre-phasor, phase encoding and  slice refocusing gradients
for i=101:150
    gradAmp(1,i) =  -gprephase; %X gradients in Tesla per meter
    gradAmp(2,i) =  Gp; %Y gradients in Tesla per meter %% should this be positive and negative?
    gradAmp(3,i) =  -gAmprefoc; %Z gradients in Tesla per meter         
end



%% Question B:

%Generate the ADC event
adc2= zeros(1,152);

for i=(1:48)
    adc1(i)=i;
end

adc=[adc2,adc1];

for i = 150:200
    gradAmp(2,i) = .001;
end

%% Question C:Display the complete GRE sequence diagram (RF, Gx,Gy,Gz, ADC) in a single figure. 
%our ADC will look like a slope. Make sure to add the appropriate labels and dimensions to the figure.

figure
subplot(3,1,1)
plot(time,rfPulse)
title('RF pulse') % title
ylabel('%B1+ in Tesla') % y measure label
xlabel('time'); % x axis label

subplot(3,1,2)
plot(time,gradAmp)
title('Gradients') % title
ylabel('Tesla per meter') % y measure label
xlabel('time'); % x axis label
legend('Gx','Gz','Gy');% legend

subplot(3,1,3)
plot(time,adc)
title('ADC') % title
ylabel('amplitude') % y measure label
xlabel('time'); % x axis label


%% Question B:

kSpace  = zeros(size(T1,1),size(T1,2));

tic
for trIndex=1:size(T1,2)

disp(trIndex);    

%update the phase endcoding gradient
for t=101:150
    gradAmp(2,i) =  t; %Y gradients in Tesla per meter
end

for k=1:xSteps  
    for j=1:ySteps
        for i=1:zSteps       
            dB0 = pos(:,k,j)'*gradAmp(:,1); %dB0 = pos(:,k,j)'*gradAmp(:,1);
            [mT,mZ] =  bloch(dt, dB0,rfPulse(1), T1(k,j),T2(k,j), 0, 1);   % start from fully relaxed spin state
            for t=2:nTimeSteps %t starts at 2
                dB0 = pos(:,k,j)'*gradAmp(:,t); %dB0 = pos(:,k,j)'*gradAmp
                [mT,mZ] =  bloch(dt, dB0,rfPulse(t), T1(k,j),T2(k,j), mT, mZ); 
                if(adc(t)>0)
                    %Sum the signal over all spins and store in k-space
                    %Don't forget to wigh the signal with the PD.
                    kSpace(round(adc(t)),trIndex) = kSpace(round(adc(t)),trIndex)+mT*PD(k,j); %kSpace(24,round(adc(t))) = kSpace(24,round(adc(t)))+mT*PD(k,j);
                end    
            end  %end of time loop            
        end %end of i loop
    end %end of j loop
end %end of k loop
end %end of TR loop
toc

%% Question C:

img = fftshift(ifft2(fftshift(kSpace*sqrt(size(kSpace,1)*size(kSpace,2)))));
figure(1); 
subplot(1,2,1);imagesc(abs(img)); colormap gray; axis equal tight off;title('image')
subplot(1,2,2);imagesc(log(abs(kSpace))); colormap gray; axis equal tight off;title('k-space')


% we got a PDW image with long TR and a short TE
