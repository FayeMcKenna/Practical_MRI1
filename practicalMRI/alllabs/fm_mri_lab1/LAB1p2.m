%% Exercise 1.2
clear all; clc; close all; %clean up

tmp = matlab.desktop.editor.getActive;  % get location of this script
cd(fileparts(tmp.Filename));            % set working directory to same


rfPulse      = [ 0,0.5, 1, 1, 1, 1, 1, 1, 1, 1, 0.5, 0]*10^-6; %B1+  in Tesla
timeLine     = [ 0,  1, 2, 3, 4, 5, 6, 7, 8, 9,  10,11]*10^-7; %Time in seconds



%Question A:

% waveform plot
plot(timeLine,rfPulse) 
title('RF Waveform')
xlabel('Time in seconds')
ylabel('B1+  in Tesla')
legend('rfPulse'); 

%Question B:

%The durration of the RF pulse is 11* 10-7 seconds in total. Maybe
%the MR system could accurately produce this shape, it is a rectangular
%shape which is common but has a steep gradient.
% range from .1 - 3 ms. 


%Question C:

save('rfPulse') 



