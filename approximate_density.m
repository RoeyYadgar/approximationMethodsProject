function [f_approx,c] = approximate_density(Y,x_j,gamma,alpha,lambda)
%Compute the density estimator given the samples Y, with the corresponding values for x_j,gamma,alpha and lambda  

kernel_L2_prod = korKernel(x_j,x_j,gamma.^2,alpha*2); %Inner product of the kernel in l2 is sum over l of (beta_l)^2 * phy_l(x)*phy_l(x') (which means gamma for korobov kernel gamma is squared and alpha is doubled) 

%Compute the coefficients c by solving Ac = b where A_jk = <K(x_j,),
%K(x_k,)>_L2 + lambda * K(x_j,x_k), bj = sum_i(K(x_j,Y_i))/M
A = kernel_L2_prod + lambda * korKernel(x_j,x_j,gamma,alpha);
b = 1/length(Y) * sum(korKernel(Y,x_j,gamma,alpha))';

c = A^(-1)*b;

f_approx = density_approximation(x_j,c,gamma,alpha); 


end

      