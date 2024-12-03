function [B]= LETFpriceEffective (S, x, A0, R, t)
    %This function calculates the LETF over the holding period "t"
    % S: Index Prices, n x nsims matrix
    % x: Leverage
    % A0: Starting LETF price
    % R: Simple returns of S over the holding period "t"
    % t: Holding period


    A = zeros(t, 1);     % Inizialing the LETF price matrix
    A(1) = A0;

    for i = 1 : t - 1
        A(i+1) = A(i) .* (1 +  x * R(i)); %Calculating the LETF prices
    end

    B = A(t); 
end