function [ mtOut, mzOut  ] = bloch(dt,dB0,B1,T1,T2,mxy,mz) %1sec

   gammaBar = 2*pi*42.577*10^6; %Hz per Tesla
 
   % should go equation 11 and 12, then 10 and 13, then 14 and 15, then 16
   mtOut = cos(gammaBar*dt*imag(B1)*real(mxy))+ 1i *cos(gammaBar*dt*real(B1)*imag(mxy)); %equation 12
   mtOut = mtOut+1i*sin(gammaBar*dt*B1*mz); %equation 11
   mtOut = mtOut*exp(1i*gammaBar*dt*dB0); %equation 10 updated transverse magnetization
   mtOut = mtOut*exp(-dt/T2); %equation 13
   
   mzOut = cos(gammaBar*B1*dt)*mz; %equation 14
   mzOut = mzOut+ sin(gammaBar*dt*imag(B1))*real(mxy)-sin(gammaBar*dt*real(B1))*imag(mxy); %equation 15
   mzOut = mzOut*exp(-dt/T1)+(1-exp(-dt/T1)); % equation 16

   
end

%calulate mtOut and mzOut:
   