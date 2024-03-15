function samples = generateSamp(f,support,M,d,x_j)
%Generate Samples using the acceptance-rejection method:
%generate a sample uniformly on the support of f, and accept the sample
%with probablity f(x)/sup(f) (if sample is rejected generate a new sample)

a = support(1); b = support(2);
supfx = max(f(x_j))*1.1; %upper bound on the maximal value of f 
samples = zeros(M,d);

if(M<10) %for low number of samples, generate sample one by one
    for i = 1:M
        sampleAccepted = false;
        while ~sampleAccepted
            x = rand(1,d)*(b-a)+a;
            u = rand();
            if(u < f(x)/supfx)
                sampleAccepted = true;
                samples(i,:) = x;
            end
        end
    end
else %for high amount of samples,generate M samples, and call the function recursivly with the number of samples that got rejected
    x = rand(M,d)*(b-a)+a;
    u = rand(M,1);
    sample_acceptance = u < f(x)/supfx;
    accepted_samples = sum(sample_acceptance);
    samples(1:accepted_samples,:) = x(sample_acceptance,:);
    samples(accepted_samples+1:end,:) = generateSamp(f,support,M-accepted_samples,d,x_j);
end
    
    