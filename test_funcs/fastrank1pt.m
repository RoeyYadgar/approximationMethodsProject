function [z, e2] = fastrank1pt(n, s_max, omega, gamma, beta)

if ~isprime(n), error('n must be prime'); end
z = zeros(s_max, 1);
e2 = zeros(s_max, 1);

m = (n-1)/2;           % assume the $\omega$ function symmetric around $1/2$
E2 = zeros(m, 1);      % the vector $\tilde{\vec{E}}^2$ in the text
cumbeta = cumprod(beta);

g = generatorp(n);     % generator $g$ for $\{1, 2, \ldots, n-1\}$
perm = zeros(m, 1);    % permutation formed by positive powers of $g$
perm(1) = 1; for j=1:m-1, perm(j+1) = mod(perm(j)*g, n); end
perm = min(n - perm, perm);    % map everything back to $[1, n/2)$
psi = omega(perm/n);   % the vector $\vec{\psi}'$
psi0 = omega(0);       % zero index: $\psi(0)$
fft_psi = fft(psi);

q = ones(m, 1);        % permuted product vector $\vec{q}'$ (without zero index)
q0 = 1;                % zero index of permuted product vector: $q(0)$

for s = 1:s_max
    % step 2a: circulant matrix-vector multiplication
    E2 = ifft(fft_psi .* fft(q));
    E2 = real(E2); % remove imaginary rounding errors
    % step 2b: choose $w_s$ and $z_s$ which give minimal value
    [min_E2, w] = min(E2); % pick index of minimal value
    if s == 1, w = 1; noise = abs(E2(1) - min_E2), end
    z(s) = perm(w);
    % extra: we want to know the exact value of the worst-case error
    e2(s) = -cumbeta(s) + ( beta(s) * (q0 + 2*sum(q)) + ...
               gamma(s) * (psi0*q0 + 2*min_E2) ) / n;
    % step 2c: update $\vec{q}$
    q = (beta(s) + gamma(s) * psi([w:-1:1 m:-1:w+1])) .* q;
    q0 = (beta(s) + gamma(s) * psi0) * q0;
    fprintf('s=%4d, z=%6d, w=%6d, e2=%.4e, e=%.4e\n', ...
             s, z(s), w, e2(s), sqrt(e2(s)));
end


function g = generatorp(p)

if ~isprime(p), error('p is not a prime'); end;

primef = unique(factor(p-1));
g = 2; i = 1;
while i <= length(primef)
    if powmod(g, (p-1)/primef(i), p) == 1
        g = g + 1; i = 0;
    end;
    i = i + 1;
end;


function y = powmod(x, a, n)

if ~usejava('jvm')
    y = 1; u = x;
    while a > 0
        if mod(a, 2) == 1
            y = mod(y * u, n); % this could overflow
        end;
        u = mod(u * u, n); % this could overflow
        a = floor(a / 2);
    end;
else
    x_ = java.math.BigInteger(num2str(x));
    a_ = java.math.BigInteger(num2str(a));
    n_ = java.math.BigInteger(num2str(n));
    y = x_.modPow(a_, n_).doubleValue();
end;