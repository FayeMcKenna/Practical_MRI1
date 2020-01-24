function [ mtOut, mzOut  ] = bloch(dt,dB0,B1,T1,T2,mt,mz) %1sec

   gamma = 42.577*10^6; %Hz per Tesla
   
   beta = 2*pi*gamma*B1*dt;
   
   mtOut = cos(imag(beta))*real(mt) + 1i*cos(real(beta))*imag(mt);         %resual mt after rf rotation
%   mtOut = mtOut + 1i*sin(2*pi*gamma*abs(B1)*dt)*exp(-1i*angle(B1))*mz;     %newly created transverse magnetization
   mtOut = mtOut + 1i*sin(2*pi*gamma*B1*dt)*mz;   
    
   mzOut =  cos(abs(beta))*mz;                                             %resual mz after rf rotation
   %mzOut =  mzOut - sin(imag(beta))*real(mt) - sin(real(beta))*imag(mt);   %newly created logitudinal magnetization
   mzOut =  mzOut + sin(imag(beta))*real(mt) - sin(real(beta))*imag(mt);

   mtOut = mtOut*exp(-2i*pi*gamma*dB0*dt);                                 %B0
   mtOut = mtOut*exp(-dt/T2);                                              %T2
   mzOut = mzOut*exp(-dt/T1) + (1-exp(-dt/T1));                            %T1 %assuming m0=1

end