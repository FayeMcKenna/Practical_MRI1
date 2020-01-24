function [ mxyOut, mzOut ] = complexlargetipangle(dB0,B1,mxy,mz)

    dt    = 10^-7;       %0.1 micro second
    gamma = 42.577*10^6; %MHz per Tesla
    mxyOut = mxy - 1i*gamma*dt*(dB0*mxy+B1*mz); %
    mzOut = mz + gamma*dt*(real(B1)*imag(mxy)-imag(B1)*real(mxy));
%   m = [mxyOut,mzOut]; %update m accoring to eqation 8
end
    
    
    

% function [ mxyOut, mzOut ] = complexlargetipangle(dB0,B1,mxy,mz);


%     mOut = [mxy,mz]; %update m accoring to eqation 8
%     mxy = mxy - gamma*dt(dB0(mxy)+B1(mz)); %*1i
%     mz = mz + gamma*dt((B1)*mxy)-(B1*mz);

%     mxy = mxy - 1i*gamma*dt(dB0(mxy)+B1(mz)); %*1i
%     mz = mz + gamma*dt(Re(B1)*Im(mxy)-Im(B1)*Re(mxy));
%     mOut = [mxy,mz]; %update m accoring to eqation 8
