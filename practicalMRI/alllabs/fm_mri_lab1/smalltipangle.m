
function [ m ] = smalltipangle(dB0,B1,m)    
m0    = 1; % gamma * dt * m0 + B1 + m
dt    = 10^-7;       %0.1 micro second
gamma = 42.577*10^6; %MHz per Tesla
m  = 1i*dt*gamma*(B1*m0-dB0*m) + m ; % m  = 1i*dt*gamma*(B1*m0-dB0*m) + m ;   
    % m = eq. 4 form the text

end