function I = inner_prod_kernel_l2(f,d,x_j,gamma,alpha)
%Compute the L2 inner product between the function f and the korobov kernel
%K(x_j,y)
%(assuming the function is a product of functions of each dimension)

I = ones(size(x_j,1),1);
M = 10e3;
x_points = linspace(0,1,M);

x_points_bern = (0:0.001:1)';
bern = bernoulli(alpha,x_points_bern);

if(d == 1)
    f_samps = f(x_points);
    %evalute kernel value with linear interpolation of the correspoind
    %beroulli polynomial
    korKernel_i = 1 + gamma*((2*pi)^alpha)/((-1)^(0.5*alpha+1)*factorial(alpha))*interp1(x_points_bern,bern,mod(x_points - x_j,1));
    I = (korKernel_i*f_samps')/M;
else 
    %since both the target funcion f and the kernel are product of
    %functions in each dimension the L2 inner product is the product of the
    %inner product in each dimension
    for i = 1:d
        f_samps = interp1(f.x_points,f.dataTable(:,i),x_points);
        %evalute kernel value with linear interpolation of the correspoind
        %beroulli polynomial
        korKernel_i = 1 + gamma(i)*((2*pi)^alpha)/((-1)^(0.5*alpha+1)*factorial(alpha))*interp1(x_points_bern,bern,mod(x_points - x_j(:,i),1));
        I = I.*(korKernel_i*f_samps')/M;
    end
    
end
    
    
    
    
