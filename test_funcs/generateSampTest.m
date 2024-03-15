function Y = generateSampTest(fTable,M,d)

Y = zeros(M,d);
L = size(fTable,1);
x_points = 0:0.01:1;
for i = 1:d
    cdf = (cumsum(fTable(:,i)) - fTable(1,i))/(L-1);
    for j = 1:M
        u  = rand();
        b2 = find(cdf - u >= 0,1); b1 = b2-1;
        Y(j,i) = x_points(b1) + (u - cdf(b1))/(cdf(b2)-cdf(b1)) * (x_points(b2) - x_points(b1));
    end
end
        
        
    
    


