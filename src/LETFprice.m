function [A]= LETFprice (S, x, A0)
    %This function calculates the LETF

    % S: Index Prices, n x nsims matrix
    % x: Leverage
    % A0: Starting LETF price

    R = PriceToReturn(S);   % Calculating the simple returns of S

    A = zeros(size(S));     % Inizialing the LETF price matrix
    A(1,:) = A0;

    n = size(S,1);          % Extracting the first dimension of S

    for i = 1 : n-1
        A(i+1,:) = A(i,:) .* (1 + x * R(i,:)); %Calculating the LETF prices
    end

end