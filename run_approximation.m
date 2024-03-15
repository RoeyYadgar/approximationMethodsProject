function [miseError,f_approx] = run_approximation(f,M,x_j,lambda,d,gamma,alpha)
%run a 'simulation' of estimating the density function f and compute the
%MISE error, repeat the simulation until the averaged MISE error is accurate
%enough

if(class(f) == "func_interpolation") %if f is of type func_interpolation then convert it into a function handle
    fTable = f.dataTable;
    f = @(x) f.interpolate(x);
end

S = 1;

miseError = 0;
while(length(miseError) == 1 || sqrt(var(miseError)/S)*2 > 0.1*mean(miseError)) %repeat simulating until the standard deviation of the averaged mise error satisfies 2*sigma > 0.1 * E(MISE)
    %Generate samples, compute estimator, and calculate MISE error
    display("Running Iteration: #" + num2str(S));
    Y = generateSamp(f,[0,1],M,d,x_j);    
    [f_approx,~] = approximate_density(Y,x_j,gamma,alpha,lambda);   
    miseError(S) = mise(f,f_approx,x_j,d);
    S = S + 1;
end
miseError = mean(miseError);



