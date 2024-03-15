classdef func_interpolation
    %Inerpolates a high-dimensional function (in the hypersphere [0,1]^d) that is a product of functions
    %of each dimension (f(x1,..,xn) = Prod(g_i(x_i))
    properties
        d
        dataTable
        x_points;
    end
    
    methods
        
        function f = func_interpolation(dataTable,x_points) %Constructor function
            d = size(dataTable,2);
            f.d = d;
            f.dataTable = dataTable;
            f.x_points = x_points;
        end
    
        function val = interpolate(f,x)
            val = ones(size(x,1),1);
            for i = 1:f.d %Compute product of linear interpolation in each dimenison
                val = val.*interp1(f.x_points,f.dataTable(:,i),x(:,i));
            end

        end
    end
    
    
end
