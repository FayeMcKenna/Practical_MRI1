%% Exercise 1.7
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

%Allocate the memory needed
rfPulse     = zeros(1,15000); %variable to hold a RF waveform with 10000 samples
rfPulse2     = zeros(1,15000); %variable to hold a RF waveform with 10000 samples
gradAmp     = zeros(1,15000); %variable to hold a gradient waveform with 10000 samples
time        = zeros(1,15000); %variable to hold 10000 time points


posZ       = zeros(1,400); %variable to hold 400 positions allong the z direction
mFinal     = zeros(1,400); %variable to hold the final magnetization calculated for each position
mFinal2     = zeros(1,400); %variable to hold the final magnetization calculated for each position
dB0AtPosJ  = 0 ; % klunge


%Generate a list of sampling points allong the z direction
for i=1:400 %i starts at 1 go's to 200
    posZ(i)  = (i-200)*10^-4; %Distance from iso center in meters
end  

%Question A:Generate wave form 1 ms with 2 side lobes and slice selective gradient

for i=1:10000 %i starts at 1 go's to 10000
    rfPulse(i) = ((sin(pi*i/12000)^2)*sinc(pi*(i-6000)/6000)*10^-5); %B1+ in Tesla ((sin(pi*i/20000)^2)*sinc(pi*(i-10000)/5000)*10^-5)
end


%Generates the time line for plotting
for i=1:15000 %i starts at 1 go's to 15000    
    time(i)    = i*10^-4; %Time in seconds 
end

% slice select gradient to excite 1 cm thick
%part1
%N = time * BW is the TBW product equation
N = 4;
timePulse = 1*10^-4;
bW = N/timePulse;
%part2
gamma = 42.577*10^6;
dZ = .01;
gZa= bW/(gamma *dZ);

%slice selection gradient 
for i=1:10000 %i starts at 1 go's to 10000
    gradAmp(i) =  gZa;                       %Tesla per meter
end

% add a slice re-focusing gradient to the wave form 1 cm thick
for i=10001:15000 %i starts at 1 go's to 10000
    gradAmp(i) = -gZa; % applying opposite
end


% use Hann window funtion 
l = 15000 ;%lenght of the pulse
h = hann(l); % use hann function
hfunc = transpose(h); % transpose to match rfPulse array
rfPulse = times(hfunc,rfPulse); % create the new hann function rfPulse

%Question B:plot in one figure RF amplitude, RF phase and gradient amplitude by time,
figure
subplot(3,1,1) % RF amplitude subplot 
x=time; % second plot x is time
y=abs(rfPulse);% second plot y is mag of rfPulse
plot(x,y,'color','r') % plot red
title('RF Amplitude') % title
ylabel('Tesla') % y measure label
xlabel('time'); % x axis label
legend('rfamp');% legend

subplot(3,1,2) % RF phase subplot 
x=time; % third plot x is time
y=angle(rfPulse);% third plot y is phase
plot(x,y,'color','g') % plot green
title('RF Phase') % title
ylabel('VF') % y measure label
xlabel('time');% x measure label
legend('RFphase');% legend

subplot(3,1,3) % gradient amplitude subplot 
x=time; % fourth plot x is time
y=gradAmp;% fourth plot y is amp 
plot(x,y,'color','b') % plot blue
title('Gradient Amplitude') % title
ylabel('Tesla') % y measure label
legend('GradAmp');% legend
xlabel('time');% x measure label

% Get magnetization of hann rfPulse

for j=1:400 %loop over different positions allong the z direction
    m =  smalltipangle(dB0AtPosJ,rfPulse(1), 0); % first time point assuming the transverse magnetization 0
    for i=2:15000 %(i starts at 2) for every point in the magnetization process
        dB0AtPosJ = posZ(j)*gradAmp(i);% dB0 of j is its pos * its gradient
         m =  smalltipangle(dB0AtPosJ,rfPulse(i), m);% total magentization is updated for one time in rfPulse
    end
    mFinal(j) = m; %magnetization calculated for each j point
end    

%Question C: Simulate the slice profile and plot the results. In addition, 
% load the results from exercise 1.6 and plot these in the same window. 
%What has changed compared to Exercise 1.6? Explain
% What has changed is that the side lobes are smoothed outa little more because the hann
% window function works to filter the data. 

% comparison of Hann vs not Hann funciton on rfPulse
load('q6phase');
load('q6amp');

figure
subplot(4,1,1)
plot(posZ,abs(mFinal));
ylabel('Tesla') % y measure label
legend('magnetization');% legend
xlabel('position');% x measure label

subplot(4,1,2);
plot(posZ,angle(mFinal));
ylabel('VF') % y measure label
legend('phase');% legend
xlabel('position');% x measure label

subplot(4,1,3)
plot(posZ,q6amp);
ylabel('Tesla') % y measure label
legend('magnetization');% legend
xlabel('position');% x measure label

subplot(4,1,4);
plot(posZ,q6phase);
ylabel('VF') % y measure label
legend('phase');% legend
xlabel('position');% x measure label

%Question D:What is the slice thickness produced by this excitation? What
%is the RF bandwidth? The slice thickness is 1cm, the BW is 8.8883e-05
%Is this a nice slice profile? %How could we make it even better?  This is a nice profile but it could use
%additional filtering to make it a more 'ideal' pulse with sharper slopes.

%slice thickness

k = max(mFinal)*.5; % the maximum of the final magnetization /2 (for area)

thick =find(abs(k)> 0); % find where within that area is greater than 0
thick = sum(thick); % sum together to find slice thickness 


%bandwidth

bandW = 1/sum(time);

recSig = abs(sum(mFinal)); % how much signal you would receive 

%Question E: Make the necessary changes to the RF pulse to move the slice 1cm from iso center. 
% gamma = 42.577*10^6; %Hz per Tesla
% gradZ = gZa;%Tesla per meter

% find deltaf using equation to select frequencies 1 cm from isocenter

deltaF=(gamma*gZa*.01); % frequencies 1cm away from current location we want

%Generate the RF waveform for increasing N from 4 to 8
for i=1:10000 %i starts at 1 go's to 10000
    rfPulse2(i) = ((sin(pi*i/12000)^2)*sinc(pi*(i-6000)/3000)*10^-5); %B1+ in Tesla ((sin(pi*i/20000)^2)*sinc(pi*(i-10000)/5000)*10^-5)
end

%Generates the time line for plotting
for i=1:15000 %i starts at 1 go's to 15000    
    time(i)    = i*10^-4; %Time in seconds .67
end

% use Hann window funtion 
l = 15000 ;%lenght of the pulse
h = hann(l); % use hann function
hfunc = transpose(h); % transpose to match rfPulse array
rfPulse2 = times(hfunc,rfPulse2); % create the new hann function rfPulse

for j=1:400 %loop over different positions allong the z direction
    m =  smalltipangle(dB0AtPosJ,rfPulse2(1), 0); % first time point assuming the transverse magnetization 0
    for i=2:15000 %(i starts at 2) for every point in the magnetization process
        dB0AtPosJ = posZ(j)*gradAmp(i);% dB0 of j is its pos * its gradient
         m =  smalltipangle(dB0AtPosJ,rfPulse2(i), m);% total magentization is updated for one time in rfPulse
    end
    mFinal2(j) = m; %magnetization calculated for each j point
end    

recSig2 = abs(sum(mFinal2));
%In addition, why does the slice profile start to look a little asymmetric? Do you think this is real or an artifact from the simulation?

figure
subplot(2,1,1)
plot(posZ,abs(mFinal2));
ylabel('Tesla') % y measure label
legend('magnetization');% legend
xlabel('position');% x measure label

subplot(2,1,2);
plot(posZ,angle(mFinal2));
ylabel('VF') % y measure label
legend('phase');% legend
xlabel('position');% x measure label


% Question F : What is the slice thickness produced by this excitation? Is
% this a nice slice profile? 1 cm thick still an ok slice profile

k2 = max(mFinal2)*.5 % the maximum of the final magnetization /2 (for area)

thick2 =find(abs(k2)> 0) % find where within that area is greater than 0
thick2 = sum(thick2) % sum together

% slice profile is an ok slice excitation profile, ideally
% rectangular but with slight slopes. 
% this looks asymmetrical because it doesnt not start at isocenter (m0).



%Question F:Create an RF waveform that excites two slices simultaneously each ±1cm from iso center.
%Try to design your pulse such that the signal from one slice is 180  out of phase compared to the other.

