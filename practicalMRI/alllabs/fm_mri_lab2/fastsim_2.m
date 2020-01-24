function [ mtOut, mzOut ] = fastsim(dB0,rf,mT,mZ);

    dt    = 10^-7;       %0.1 micro second
    gamma = 42.577*10^6; %MHz per Tesla
    %step 1
    mtOut = (dB0(mxy)+B1(mz));
    mzOut = (Re(B1)*Im(mxy)-Im(B1)*Re(mxy));
    %step 2
    mtOut = mT - 1i*gamma*dt(mtOut);
    mzOut = mZ + gamma*dt(mzOut);
    
    %mxyOut = mxy - 1i*gamma*dt(dB0(mxy)+B1(mz)); %*1i
    %mzOut = mz + gamma*dt(Re(B1)*Im(mxy)-Im(B1)*Re(mxy));
%   m = [mxyOut,mzOut]; %update m accoring to eqation 8
end
    
    