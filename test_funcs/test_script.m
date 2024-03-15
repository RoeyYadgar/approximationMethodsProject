%%
omega = inline('2*pi^2*(x.^2-x+1/6)');
omega2 = @(x) ((2*pi)^alpha2)/((-1)^(alpha2/2+1)*factorial(alpha2))*bernoulli(alpha,mod(x,1));
gamma = gamma1' ; beta = ones(d, 1);
[z, e2] = fastrank1pt(N, d, omega, gamma, beta);

%%
x_j1 = mod((1:N)'.*z1/N,1);

M = 10^4;
lambdas = 10.^(0:-0.5:-7);
miseError1 = zeros(1,length(lambdas));
miseError2 = zeros(1,length(lambdas));
for i = 1:length(lambdas)
    [miseError1(i),~] = run_approximation(fx,M,x_j1,lambdas(i),d,gamma1,alpha1);
    [miseError2(i),~] = run_approximation_determnisitic(fx,M,x_j1,lambdas(i),d,gamma1,alpha1);
    %[miseError1(i),~] = run_approximation(f2,M,(0:N-1)'/N,lambdas(i),1,[1],2);
end 
plot(log10(miseError1))
hold on
plot(log10(miseError2))
clc