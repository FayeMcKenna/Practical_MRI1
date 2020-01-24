% add a slice re-focusing gradient to the wave form 

for i=10001:15000 %i starts at 1 go's to 10000
    rfPulse(i) = rfPulse(i)*0; %B1+ in Tesla
    gradAmp(i) = - 5*10^-2; % applying opposite
end

%with a thickness of 1.

k = max(mFinal)*.5; % the maximum of the final magnetization /2 (for area)

thick =find(abs(k)> 0); % find where within that area is greater than 0
thick = sum(thick); % sum together

%bandwidth

BW = 1/sum(time);
