classdef density_approximation
    %Inerpolates a high-dimensional function (in the hypersphere [0,1]^d) that is a sum of the korKernel
    %(a periodic kernel in the hypersphere and also is a product of
    %function at each dimension)
    
    properties
        d
        x_j
        kernelInterpolator
        c
    end
    
    methods
              
        function f_approx = density_approximation(x_j,c,gamma,alpha) %Constructor function

            d = size(x_j,2);
            x_points = (0:0.01:1)';
            %Compute kernel table (used for linear interpolation at each
            %dimension for fast evaluation)
            kernelTable = zeros(length(x_points),d);
            for i = 1:d
                kernelTable(:,i) = korKernel(0,x_points,gamma(i),alpha);
            end

            f_approx.d = d;
            f_approx.x_j = x_j;
            f_approx.kernelInterpolator = func_interpolation(kernelTable,x_points);
            f_approx.c = c;

        end
        
        function val = approximate(f_approx,x) %Evaluate the estimator value at the points: x
            val = zeros(size(x,1),1);
            for j = 1:size(f_approx.x_j,1)%Sum over each x_j the value for K(mod(x-x_j,1))
                val = val + f_approx.c(j)*f_approx.kernelInterpolator.interpolate(mod(x - f_approx.x_j(j,:),1));
            end

        end
    end
    
    
end
