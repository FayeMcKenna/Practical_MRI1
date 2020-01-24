function [ m ] = largetipangle(dB0,B1,m)
    dt    = 10^-7;       %0.1 micro second
    gamma = 42.577*10^6; %MHz per Tesla
    mat = [0 dB0 -B1; -dB0 0 B1; B1 B1 0];
    m = m + gamma*dt*mat*m; % update m accoridng to equation 6 % do I need to make B1y and B1x different?
 
end


