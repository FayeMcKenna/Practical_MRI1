for nn = 1:12
    randomfor(:,nn) = randi([max(A(:,nn)),min(A(:,nn))]); 
end

r = randi([.00000025 -.00000025],1,10000);


r=.00000025.*rand()

for i=1:5000 %i starts at 1 go's to 200 *333.3
    posZ(i)  = i*.00000025*rand(); %Distance from iso center in meters
end

for i=5001:10000 %i starts at 1 go's to 200 *333.3
    posZ(i)  = i*-.00000025*rand(); %Distance from iso center in meters
end

