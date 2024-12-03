function [S]= LogNormal(P0, mrf, sigma, n, T, lambda, mrfj, sigmaj)

    % This function generates a simulated Index price through a lognormal stochastic process
    % with jumps

    % mrf:      Index expected return
    % sigma:    Index volatility
    % n:        Number of steps
    % T:        Period of simulation

    % lambda:   Poisson parameter
    % mrfj:     Poisson mean
    % sigmaj:   Poisson standard deviation
    % P0:       Wanted starting Index price

    h= T/n;                                % Magnitude of  temporal steps, delta t
    k= exp(mrfj) - 1;                      % Expected percentage jump
    S(1)= P0;                              % Setting the simulated Index Starting pric to P0 
    rmean= (mrf - lambda*k - (sigma^2)/2); % "Drift" of log returns, consequence of Ito's lemma

    for i= 1 : n
        m= PoissonGen(lambda);             % Generating Poisson distributed variable
        meanj= (mrfj - (sigmaj^2)/2);      % "Drift" of the jumps; consequence of Ito's lemma
        rsdj = sum(randn(1,m));            % Generating m random variables and summing them up

        % We cut the generation of the first exponent argument in two lines
        p= rmean*h/1260;                        
        pp= sigma * randn * sqrt(h/1260) ;

        % we cut the generation of the nsims simulations at step i+1 in two lines
        S(i+1) = S(i) .* exp(p+pp);                         % Log return step
        S(i+1)= S(i+1).* exp(m * meanj + sigmaj .* rsdj);   % Adjusting for extra random jumps
    end

end