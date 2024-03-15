function error = mise(f,f_approx,x_j,d)
%Compute the MISE error between the density function f and the estimator
%f_approx

L = 100;
N = size(x_j,1);
p = net(sobolset(d),L);

errorSum = 0;
for i = 1:L
    points = mod(x_j + p(i,:),1);
    errorSum = errorSum + sum((f(points)-f_approx.approximate(points)).^2);
end
error = errorSum/(N*L);
    
    