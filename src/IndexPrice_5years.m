function [S]= IndexPrice_5years (P0, mrf, sigma, nsims, n, T, lambda, mrfj, sigmaj)

    % This function generates nsims simulated Index prices through the LogNormal function

    % mrf:      Index expected return
    % sigma:    Index volatility
    % nsims:    Number of simulations
    % n:        Number of steps
    % T:        Periods of simulation

    % lambda:   Poisson parameter
    % mrfj:     Poisson mean
    % sigmaj:   Poisson standard deviation
    % P0:       Starting Index price

    for i = 1 : nsims
        S(:,i)= LogNormal_5years(P0, mrf, sigma, n, T, lambda, mrfj, sigmaj);
    end
end