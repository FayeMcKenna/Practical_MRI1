%% Exercise 1.6
clear all; clc; close all; % clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

%Allocate the memory needed
rfPulse     = zeros(1,15000); %variable to hold a RF waveform with 10000 samples
gradAmp     = zeros(1,15000); %variable to hold a gradient waveform with 10000 samples
time        = zeros(1,15000); %variable to hold 10000 time points


posZ       = zeros(1,400); %variable to hold 400 positions allong the z direction
mFinal     = zeros(1,400); %variable to hold the final magnetization calculated for each position

dB0AtPosJ  = 0 ; % klunge


%Generate a list of sampling points allong the z direction
for i=1:400 %i starts at 1 go's to 200
    posZ(i)  = (i-200)*10^-4; %Distance from iso center in meters
end  

%Question A: use MATLAB to generate a 1ms Sinc pulse with 2 side lobes on either side 
%and a slice selection gradient & a refocusing gradient to excite 1cm thick slice

for i=1:10000 %i starts at 1 go's to 10000
    rfPulse(i) = ((sin(pi*i/12000)^2)*sinc(pi*(i-6000)/6000)*10^-5); %B1+ in Tesla ((sin(pi*i/20000)^2)*sinc(pi*(i-10000)/5000)*10^-5)
end


%Generates the time line for plotting
for i=1:15000 %i starts at 1 go's to 15000    
    time(i)    = i*10^-4; %Time in seconds .67
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


% Question B: 
% plot in one figure RF amplitude, RF phase and gradient amplitude by time,
% and subplots in column
figure
% first subplot
subplot(3,1,1) % second subplot in column showing rf amp and rf phase
plot(time,angle(rfPulse),'color','r') % plot red
title('RF phase') % title
ylabel('VF') % y measure label
xlabel('time'); % x axis label
legend('RF phase');% legend

subplot(3,1,2) % second subplot in column showing gradient amp
plot(time,abs(rfPulse),'color','g') % plot green
title('RF amplitude') % title
ylabel('%B1+ in Tesla') % y measure label
xlabel('time');% x measure label
legend('RF amplitude');% legend

subplot(3,1,3) % third subplot in column showing gradient amp
plot(time,gradAmp,'color','b') % plot blue
title('gradient amplitude') % title
ylabel('Tesla per meter') % y measure label
legend('gradAmp');% legend
xlabel('time');% x measure label
 

%find magnetization

for j=1:400 %loop over different positions allong the z direction
    m =  smalltipangle(dB0AtPosJ,rfPulse(1), 0); % first time point assuming the transverse magnetization 0
    for i=2:15000 %(i starts at 2) for every point in the magnetization process
        dB0AtPosJ = posZ(j)*gradAmp(i);% dB0 of j is its pos * its gradient
         m =  smalltipangle(dB0AtPosJ,rfPulse(i), m);% total magentization is updated for one time in rfPulse
    end
    mFinal(j) = m; %magnetization calculated for each j point
end    


%Is this a good slice profile? It is a pretty good slice profile but not a
%sharp rectangle

figure
subplot(2,1,1)
plot(posZ,abs(mFinal));
title('Magnitude') % title
ylabel('Tesla') % y measure label
legend('magnititude');% legend
xlabel('position');% x measure label

subplot(2,1,2);
plot(posZ,angle(mFinal));
title('Phase') % title
ylabel('VF') % y measure label
legend('phase');% legend
xlabel('position');% x measure label

%Question C: Save the results in a separate file.
q6amp= abs(mFinal);
q6phase = angle(mFinal);
save('q6phase');
save('q6amp');


