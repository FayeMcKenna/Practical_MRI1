% 4ms TBW factor 5 apodized Sinc pulse, 
% sinc pulse 

clear all; clc; close all; % clean up
tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same

%Allocate the memory needed
rfPulse     = zeros(1,15000); %variable to hold a RF waveform with 10000 samples
gradAmp     = zeros(1,15000); %variable to hold a gradient waveform with 10000 samples
time        = zeros(1,15000); %variable to hold 10000 time points

% Generate wave form 1 ms with 2 side lobes and slice selective gradient
% t = 0:.00040001:6; 	% need to make this 1 ms
% dT = t(2)-t(1);
% time = [];
% Generate the RF waveform
for i=1:15000 %i starts at 1 go's to 15000
    rfPulse(i) = .05*sinc(1/6)*10-5; % exp(-((i-5000)/2000)^2)*10^-5
end

%Generates the time line for plotting
for i=1:15000 %i starts at 1 go's to 15000
    time(i)    = i*10^-7;    
end    
    
%slice selection gradient 
for i=1:10000 %i starts at 1 go's to 10000
    gradAmp(i) = 5*10^-2;                       %Tesla per meter
end

% add a slice re-focusing gradient to the wave form 1 cm thick
for i=10001:15000 %i starts at 1 go's to 10000
    gradAmp(i) = - 5*10^-2; % applying opposite
end

y = displaysequence(time,rfPulse,gradAmp);

%% Exercise 1.5 sine wave

% tmp = matlab.desktop.editor.getActive;  % get location of this script
% cd(fileparts(tmp.Filename));            % set working directory to same
% 
% %Allocate the memory needed
% rfPulse     = zeros(1,15000); %variable to hold a RF waveform with 10000 samples
% gradAmp     = zeros(1,15000); %variable to hold a gradient waveform with 10000 samples
% time        = zeros(1,15000); %variable to hold 10000 time points
% 
% 
% posZ       = zeros(1,400); %variable to hold 400 positions allong the z direction
% mFinal     = zeros(1,400); %variable to hold the final magnetization calculated for each position
% 
% dB0AtPosJ  = 0 ; % klunge
% 
% %Generate the RF waveform
% for i=1:15000 %i starts at 1 go's to 10000
%     rfPulse(i) = exp(-((i-5000)/2000)^2)*10^-5; %B1+ in Tesla
%     gradAmp(i) = 5*10^-2;                       %Tesla per meter
% end
% 
% %Generates the time line for plotting
% for i=1:15000 %i starts at 1 go's to 15000
%     time(i)    = i*10^-7;                       %Time in seconds
% end
% 
% %Generate a list of sampling points allong the z direction
% for i=1:400 %i starts at 1 go's to 200
%     posZ(i)  = (i-200)*10^-4; %Distance from iso center in meters
% end  
% 
% % Question A: 
% % add a slice re-focusing gradient to the wave form 
% 
% for i=10001:15000 %i starts at 1 go's to 10000
%     rfPulse(i) = rfPulse(i)*0; %B1+ in Tesla
%     gradAmp(i) = - 5*10^-2; % applying opposite
% end
% 
% y = displaysequence(time,rfPulse,gradAmp);


% and excite the spins using a 45 , 2ms TBW factor 3 Sinc pulse. 