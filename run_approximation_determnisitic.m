function [miseError,f_approx] = run_approximation_determnisitic(f,M,x_j,lambda,d,gamma,alpha)
%run a simulation of estimating the density function f and compute the
%MISE error, with 'infinite' amount of samples (by replacing the values of 
%the vector b with the L2 inner product between the denstiy function and
%the kernel at points x_j, instead of averaging of M samples)

b = inner_prod_kernel_l2(f,d,x_j,gamma,alpha); 

if(class(f) == "func_interpolation") %if f is of type func_interpolation then convert it into a function handle
    fTable = f.dataTable;
    f = @(x) f.interpolate(x);
end


kernel_L2_prod = korKernel(x_j,x_j,gamma.^2,alpha*2); %Inner product of the kernel in l2 is sum over l of (beta_l)^2 * phy_l(x)*phy_l(x') (which means gamma for korobov kernel gamma is squared and alpha is doubled) 
A = kernel_L2_prod + lambda * korKernel(x_j,x_j,gamma,alpha);

c = A^(-1)*b;
f_approx = density_approximation(x_j,c,gamma,alpha);


miseError = mise(f,f_approx,x_j,d);
