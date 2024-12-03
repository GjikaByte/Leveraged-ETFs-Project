function [R] = PriceToReturn (S)
    % This function calculates the daily (period t) returns from the daily (period t) prices
    % R: Returns vector
    % S: Prices vector
    R = ( S(2:end,:)  - S(1:end-1,:) ) ./ S(1:end-1,:)  ; 
end