function [ mtOut, mzOut  ] = fastsim(dB0,B1,mt,mz) %1sec

    dt    = 10^-7;       %0.1 micro second
    gamma = 2*pi*42.577*10^6; %Hz per Tesla

    %step 1
    mtOut = dB0*mt+B1*mz;
    mzOut = real(B1)*imag(mt)-imag(B1)*real(mt);
    %step 2
    mtOut = mt - 1i*gamma*dt*mtOut;
    mzOut = mz + gamma*dt*mzOut;


end

