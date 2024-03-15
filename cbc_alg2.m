function z = cbc_alg2(d,gamma,alpha,n)
%Calculate rank 1 lattice vector generator with the component by component
%algorithm. according to "Fast Component-by-component Construction of
%Lattice Algorithms for Multivariate Approximation
%with POD and SPOD weights", computes the same output as cbc_alg more
%effiecently

    omega_mat = zeros(n-1,n);
    for z = 1:(n-1)
        omega_mat(z,:) = omega(z,(0:(n-1)),alpha,n);
    end
    psi_mat = omega_mat.^2 - 2*zeta(2*alpha);
    
    z = zeros(1,d);
    z(1) = 1;
    P = (1 + gamma(1)*omega_mat(z(1),:)).^2;
    for s = 2:d
        w_ds = gamma(s)*P';
        v_ds = gamma(s)^2*P';
        T_ds = 1/n*(psi_mat*v_ds + 2 * omega_mat*w_ds);
        [~,minInd] = min(T_ds);
        z(s) = minInd; 
        P = P.*((1 + gamma(s)*omega_mat(z(s),:)).^2);
        
    end


end


function val = omega(z,k,alpha,n)
    val = ((2*pi)^alpha)/((-1)^(alpha/2+1)*factorial(alpha))*bernoulli(alpha,mod(k*z,n)/n);
end


