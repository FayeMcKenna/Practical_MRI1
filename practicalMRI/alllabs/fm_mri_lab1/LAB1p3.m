%% Exercise 1.3

clear all; clc; close all; %clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same


%Question B:

 load('rfPulse'); %load the rfPulse 

%Question C:
    m    = [0,0,0,0,0,0,0,0,0,0,0,0];       %Allocate memory for answer
    time = [0,1,2,3,4,5,6,7,8,9,10,11]*10^-7; %Time points for plotting

% code necessary to evaluate the magnetization of waveform as a function of time
m( 1) =  smalltipangle(0,rfPulse(1), 0 ); %1st time step
m (2) = smalltipangle(0,rfPulse(2), m(1)); %2nd time step  
m (3) = smalltipangle(0,rfPulse(3), m(2)); %3rd time step
m (4) = smalltipangle(0,rfPulse(4), m(3)); %4th time step
m (5) = smalltipangle(0,rfPulse(5), m(4)); %5th time step
m (6) = smalltipangle(0,rfPulse(6), m(5)); %6th time step
m (7) = smalltipangle(0,rfPulse(7), m(6)); %7th time step
m (8) = smalltipangle(0,rfPulse(8), m(7)); %8th time step
m (9) = smalltipangle(0,rfPulse(9), m(8)); %9th time step
m (10) = smalltipangle(0,rfPulse(10), m(9)); %10th time step
m (11) = smalltipangle(0,rfPulse(11), m(10)); %11th time step
m (12) = smalltipangle(0,rfPulse(12), m(11)); %12th time step

%Question D:

load('labSim.mat'); % load pre-simulated data

%Question E:
% verify all code together to plot simulated and pre-simulated results
% together (each curve has different color)

%plot of sim and pre-sim
x = time;
y1 = abs(m);
y2 = labSim;

figure
plot(x,y1,x,y2);

legend('sim','presim'); 
title('Magnetization of Sim and PreSim Data');
ylabel('magnetization');
xlabel('Time in *10^-7 seconds');
