function val = korKernel(x1,x2,gamma,alpha)
%Evaluate the korobov kernel value at all the pair of points x1 and x2 

val = 1;
%if the number of points the kernel value is calculated is not too big - use the bernoulli function at each point,
%otherwise use linear interpolation with 101 points to interpolate the
%value at each point,addionally since the korobov kernel is a product of
%kernels in dimension (assuming the coeffiencts gamma are product weights)
%the value of the kernel is the product of the value of the kernel in each
%dimension
if(length(x1)*length(x2) <= 1e3) 
    for j = 1:length(gamma)
        korKernel_j = 1 + gamma(j)*((2*pi)^alpha)/((-1)^(0.5*alpha+1)*factorial(alpha))*bernoulli(alpha,mod(x1(:,j)-x2(:,j)',1));
        val = val.*korKernel_j;
    end
else
    x_points = (0:0.01:1)';
    bern = bernoulli(alpha,x_points);
    for j = 1:length(gamma)
        korKernel_j = 1 + gamma(j)*((2*pi)^alpha)/((-1)^(0.5*alpha+1)*factorial(alpha))*interp1(x_points,bern,mod(x1(:,j)-x2(:,j)',1));
        val = val.*korKernel_j;
    end
end






