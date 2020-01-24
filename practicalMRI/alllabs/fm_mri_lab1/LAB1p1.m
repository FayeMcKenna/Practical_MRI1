%% Exercise 1.1
clear all; clc; close all; %clean up



%Question A:
flipAngle = .8; %flip angle in degrees.

cos(pi*flipAngle/180) %example how to use the cosine fuction in matlab

 
%Question B:

% A 8% flip angle will still produce a Mz=.1 error rate(assuming the M0=1, but anything higher will
% reduce M0.  The small tip angle approximation can be used because small excitations 
% will make Mz about equal to M0 during an RF pulse that is only produces a small net rotation. 
