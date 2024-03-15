%% Run some 'Sanity Checks' and examples 
d = 1;
res = 0.01;
alpha = 2;
N = 11;
M = 1000;
x_j = (0:(N-1))'/N;
gamma = 1;
lambda = 0.00001;
%% Regularization effect on bias
%As we decrease the regularization value lambda we should see a decrease in
%the bias term of the MISE error

lambdas = [1 1e-1 1e-2 1e-4];
miseError = zeros(size(lambdas));
f = @(x) 1+ bernoulli(4,x);
figure
x = (0:0.01:1)';
plot(f(x));
hold on
for i = 1:length(lambdas)
    lambda = lambdas(i);
    [error,f_approx] = run_approximation(f,M,x_j,lambda,d,gamma,alpha);
    miseError(i) = error;
    plot(f_approx.approximate(x))
end
title(['Regularization effect on bias - example' newline 'M = ' num2str(M) ' N = ' num2str(N) ' \alpha = ' num2str(alpha)]); 
legend('Original Density Function f(x) = 1 + Bernouli(4,x)',['\lambda=' num2str(lambdas(1))],['\lambda=' num2str(lambdas(2))],['\lambda=' num2str(lambdas(3))],['\lambda=' num2str(lambdas(4))])

%% Regularization effect on variance
%As we decrease the regularization value lambda we should see an increase in
%the variance term of the MISE error

figure
lambdas = [1 1e-4];
plot(f(x));
hold on
miseError = zeros(2,3);
legendText = ["Original Density Function f(x) = 1 + Bernouli(4,x)"];
for i = 1:length(lambdas)
    lambda = lambdas(i);
    for j = 1:3
        [miseError(i,j),f_approx] = run_approximation(f,M,x_j,lambda,d,gamma,alpha);
        plot(f_approx.approximate(x),'--')
        legendText(end+1) = string(['\lambda=' num2str(lambdas(i)) ',iteration ' num2str(j) ', error=' num2str(miseError(i,j))]);
    end
end
title(['Regularization effect on variance - example' newline 'M = ' num2str(M) ' N = ' num2str(N) ' \alpha = ' num2str(alpha)]); 
legend(legendText);

%% Smoothness of density effect on error
%Run the estimation with different 'levels' of B-splines as the density
%function. the first 2 B-splines are not in the korobov space for alpha = 4
%and hence we should see a large bias error.
%overall we see a decrease of the MISE error as the function gets smoother

N = 50;
x_j = (0:(N-1))'/N;
lambda = 0;
alpha = 4;
x = (0:0.001:1)';

f = @(x) 10*(0.45<=x & x<=0.55);
fig1 = figure;fig2 = figure;
title('Smoothness Effect on MISE Error')

figure(fig1);subplot(2,2,1)
plot(x,f(x))
hold on
[miseError,f_approx] =  run_approximation_determnisitic(f,M,x_j,lambda,d,gamma,alpha);
plot(x,f_approx.approximate(x),'--')
%legend(["Non Continuous Density Function","Non Continuous Density Function Approximation" + " , Error " + num2str(miseError)]);
title("Non Continuous Density Function Approximation" + " , Error " + num2str(miseError))
figure(fig2);subplot(2,2,1)
plot(x,f(x)-f_approx.approximate(x));title("Non Continuous Density Function Approximation Difference")


figure(fig1);subplot(2,2,2)
f = @(x) 100.*(x-0.4).*(0.4 <= x & x <= 0.5) - 100.*(x-0.6).*(0.5 < x & x <= 0.6);
plot(x,f(x))
hold on
[miseError,f_approx] =  run_approximation_determnisitic(f,M,x_j,lambda,d,gamma,alpha);
plot(x,f_approx.approximate(x),'--')
%legend(["C_{0} Density Function","C_{0} Density Function Approximation" + " , Error " + num2str(miseError)]);
title("C_{0} Density Function Approximation" + " , Error " + num2str(miseError))
figure(fig2);subplot(2,2,2)
plot(x,f(x)-f_approx.approximate(x));title("C_{0} Density Function Approximation Difference")

figure(fig1);subplot(2,2,3)
f1 = @(x) (f(x-0.05).*(0.65-x)/(0.2) + f(x+0.05).*(x-0.35)/(0.2));
plot(x,f1(x))
hold on
[miseError,f_approx] =  run_approximation_determnisitic(f1,M,x_j,lambda,d,gamma,alpha);
plot(x,f_approx.approximate(x),'--')
%legend(["C_{1} Density Function","C_{1} Density Function Approximation" + " , Error " + num2str(miseError)]);
title("C_{1} Density Function Approximation" + " , Error " + num2str(miseError))
figure(fig2);subplot(2,2,3)
plot(x,f1(x)-f_approx.approximate(x));title("C_{1} Density Function Approximation Difference")

figure(fig1);subplot(2,2,4)
f2 = @(x) (f1(x-0.05).*(0.7-x)/(0.3) + f1(x+0.05).*(x-0.3)/(0.3));
plot(x,f2(x))
hold on
[miseError,f_approx] =  run_approximation_determnisitic(f2,M,x_j,lambda,d,gamma,alpha);
plot(x,f_approx.approximate(x),'--')
%legend(["C_{1} Density Function","C_{1} Density Function Approximation" + " , Error " + num2str(miseError)]);
title("C_{2} Density Function Approximation" + " , Error " + num2str(miseError))
figure(fig2);subplot(2,2,4)
plot(x,f2(x)-f_approx.approximate(x));title("C_{2} Density Function Approximation Difference")
