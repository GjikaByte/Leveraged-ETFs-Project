function [Poisson]=PoissonGen( lambda)
    % This function exploits the inverse transform method to generate a Poisson variable

    % F: cumulative distribution, U∼Unif[0,1]
    % Sequential search for the smallest n at which F(n) <= U,
    % F(n) = P( N=0 ) + ... P( N=n )

    % Note that: P( N=k+1) = P( N=k )*lambda/k+1

    Poisson=0;
    U=rand;
    P=exp(-lambda);
    F=P;
    while U > F
        Poisson = Poisson + 1;
        P = P * lambda / Poisson;
        F = F + P;
    end

end