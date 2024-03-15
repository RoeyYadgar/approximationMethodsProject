dims = [6 15];
N = 11;
mkdir('figs');

for d = dims %Run Each Figure For Both d = 6 and d = 15
    %Compute gamma and x_j for both alpha = 2 and alpha = 4
    alpha1 = 2;
    gamma1 = 1./((1:d).^alpha1);
    z1 = cbc_alg2(d,gamma1,alpha1,N);
    x_j1 = mod((1:N)'.*z1/N,1);

    alpha2 = 4;
    gamma2 = 1./((1:d).^alpha2);
    z2 = cbc_alg2(d,gamma2,alpha2,N);
    x_j2 = mod((1:N)'.*z2/N,1);

    %In order to compute the value of the target function at each point, we
    %use linear interpolation for each dimension and take the product of
    %the result (since the function is a product of functions of each
    %dimension)
    x_points = (0:0.01:1)';
    fTable = bernoulli(4,x_points)./((1:d).^4)+1;
    fx = func_interpolation(fTable,x_points);
    %% Generate Fig 1 - 4 - MISE error as a function of M (with different values of lambda and alpha)
    M = 10.^(1:6);
    lambdas = [0.8 0.4 0.1 0.01 1e-3 1e-4];
    miseError1 = zeros(length(M),length(lambdas));
    miseError2 = zeros(length(M),length(lambdas));
    
    %Compute approximation error for each value of M,lambda and alpha
    for i = 1:length(M)
        for j = 1:length(lambdas)
            [err,~] = run_approximation(fx,M(i),x_j1,lambdas(j),d,gamma1,alpha1);
            miseError1(i,j) = err;
            [err,~] = run_approximation(fx,M(i),x_j2,lambdas(j),d,gamma2,alpha2);
            miseError2(i,j) = err;

        end
    end

    refrence = 10^(-0.5)./(M/M(1));
    %Generate the figure
    figure('Renderer', 'painters', 'Position', [10 10 1210 610])

    subplot(1,2,1)
    loglog(M,miseError1);
    hold on
    loglog(M,refrence,'black');
    xlabel('M'); ylabel('MISE');title("d = " + num2str(d) + ", \alpha = " + num2str(alpha1));
    ylim([min(miseError1(:)) max(miseError1(:))]);
    legends = ["\lambda = "  + string(num2str(lambdas')) ; "Reference rate M^{-1}"]; legend(legends);

    subplot(1,2,2)
    loglog(M,miseError2);
    hold on
    loglog(M,refrence,'black');
    xlabel('M'); ylabel('MISE');title("d = " + num2str(d) + ", \alpha = " + num2str(alpha2));
    ylim([min(miseError2(:)) max(miseError2(:))]);
    legends = ["\lambda = "  + string(num2str(lambdas')) ; "Reference rate M^{-1}"]; legend(legends);

    saveas(gcf,"figs/Mise - M plot , d = " + num2str(d) + ".png")

    %% Generate Fig 5 & 6 - Mise error as a function of lambda (with M = 10^4 and M = inf)

    M = 10^4;
    lambdas = 10.^(0:-0.5:-7);
    miseError1 = zeros(2,length(lambdas));
    miseError2 = zeros(2,length(lambdas));
    %Compute approximation error for each value of M,lambda and alpha
    for i = 1:length(lambdas)
        [miseError1(1,i),~] = run_approximation(fx,M,x_j1,lambdas(i),d,gamma1,alpha1);
        [miseError1(2,i),~] = run_approximation_determnisitic(fx,M,x_j1,lambdas(i),d,gamma1,alpha1);
        [miseError2(1,i),~] = run_approximation(fx,M,x_j2,lambdas(i),d,gamma2,alpha2);
        [miseError2(2,i),~] = run_approximation_determnisitic(fx,M,x_j2,lambdas(i),d,gamma2,alpha2);
    end 

    refrence = 0.1.*(lambdas.^2/lambdas(1)^2);
    %Generate the figure
    figure('Renderer', 'painters', 'Position', [10 10 1210 610])

    subplot(1,2,1)
    loglog(lambdas,miseError1); set(gca, 'XDir','reverse');
    hold on
    loglog(lambdas,refrence,'black');
    xlabel('\lambda'); ylabel('MISE');title("d = " + num2str(d) + ", \alpha = " + num2str(alpha1) + ",N = 11, M = 10^{4}");
    ylim([min(miseError1(:))*0.5 max(miseError1(:))*2]);
    legends = ["MISE M=10^{4}";"MISE M=\infty"; "Reference rate \lambda^{2}"]; legend(legends);


    subplot(1,2,2)
    loglog(lambdas,miseError2); set(gca, 'XDir','reverse');
    hold on
    loglog(lambdas,refrence,'black');
    xlabel('\lambda'); ylabel('MISE');title("d = " + num2str(d) + ", \alpha = " + num2str(alpha2) + ",N = 11, M = 10^{4}");
    ylim([min(miseError2(:))*0.5 max(miseError2(:))*2]);
    legends = ["MISE M=10^{4}";"MISE M=\infty"; "Reference rate \lambda^{2}"]; legend(legends);

    saveas(gcf,"figs/Mise - lambda plot , d = " + num2str(d) + ".png")

    %% Generate Fig 7 & 8 - Mise Error as a function of M (with optimal lambda 'rate') 
    M = 10.^(3:7);
    lambdas1 = 1000*M.^(-1/(1+1/alpha1));
    lambdas2 = 5000*M.^(-1/(1+1/alpha2));

    miseError1 = zeros(1,length(M));
    miseError2 = zeros(1,length(M));
    %Compute approximation error for each value of M,lambda and alpha
    for i = 1:length(M)
        [miseError1(i),~] = run_approximation(fx,M(i),x_j1,lambdas1(i),d,gamma1,alpha1);
        [miseError2(i),~] = run_approximation(fx,M(i),x_j2,lambdas2(i),d,gamma2,alpha2);    
    end

    refrence1 = 10^(0.5) * M.^(-1/(1+1/alpha1)) / M(1).^(-1/(1+1/alpha1));
    refrence2 = 10^(0.5) * M.^(-1/(1+1/alpha2)) / M(1).^(-1/(1+1/alpha2));
    %Generate the figure
    figure('Renderer', 'painters', 'Position', [10 10 1210 610])

    subplot(1,2,1)
    loglog(M,miseError1); 
    hold on
    loglog(M,refrence1,'black');
    xlabel('\M'); ylabel('MISE');title("d = " + num2str(d) + ", \alpha = " + num2str(alpha1) + ",N = 11");
    ylim([min(miseError1(:))*0.5 max(miseError1(:))*2]);
    legends = ["MISE" ; "Reference rate $M^{-\frac{1}{1+1/\alpha}}$"]; legend(legends,'Interpreter','latex');


    subplot(1,2,2)
    loglog(M,miseError2); 
    hold on
    loglog(M,refrence2,'black');
    xlabel('\M'); ylabel('MISE');title("d = " + num2str(d) + ", \alpha = " + num2str(alpha2) + ",N = 11");
    ylim([min(miseError2(:))*0.5 max(miseError2(:))*2]);
    legends = ["MISE" ; "Reference rate $M^{-\frac{1}{1+1/\alpha}}$"]; legend(legends,'Interpreter','latex');

    saveas(gcf,"figs/Mise - M optimal lambda plot , d = " + num2str(d) + ".png")

end
%% Run Sanity Checks And Examples
SanityChecks;