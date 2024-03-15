function z = cbc_alg(d,gamma,alpha,n)
%Calculate rank 1 lattice vector generator with the component by component
%algorithm. according to "Fast Component-by-component Construction of
%Lattice Algorithms for Multivariate Approximation
%with POD and SPOD weights"

    z = zeros(1,d);
    z(1) = 1;
    w_ds = zeros(n,1);
    for s = 2:d
        for k = 0:(n-1)
            w_ds(k+1) = W_ds(z(1:(s-1)),k,d,s,gamma,alpha,n);
        end
        v_ds = gamma(s)*w_ds;
        
        T_ds = zeros(n,1);
        for z_s = 0:(n-1)
           T_ds(z_s+1) = 1/n*psi(z_s,(0:(n-1)),alpha,n)*v_ds + 2/n*omega(z_s,(0:(n-1)),alpha,n)*w_ds;
        end
        [~,minInd] = min(T_ds);
        z(s) = minInd-1; 
    end


end


function val = omega(z,k,alpha,n)
    val = ((2*pi)^alpha)/((-1)^(alpha/2+1)*factorial(alpha))*bernoulli(alpha,mod(k*z/n,1));
end
function val = psi(z,k,alpha,n)
    val = omega(z,k,alpha,n).^2 - 2*zeta(2*alpha);
end

function val = W_ds(z,k,d,s,gamma,alpha,n)
    prod1 = 1;
    prod2 = 1;
    for j = 1:(s-1)
        prod1 = prod1 * (1 + gamma(j)*omega(z(j),k,alpha,n))^2;
    end
    for j = (s+1):d
        prod2 = prod2 * (1 + 2*zeta(2*alpha)*(gamma(j)^2));
    end
    val = gamma(s)*prod1*prod2;

end

