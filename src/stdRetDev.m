function [out] = stdRetDev(S, S10, S11, S12, S20, S21, t)
    % This function:
    % 1) Calculates the return of the LETFs x2,-x2,-x2,x3,-x3 for holding period of t
    % 2) Calculates the return deviation of the LETFs with naive expected return
    % 3) Returns the standard deviation of the return deviation for each LETF

    % W10: LETFs x2 return over the holding period t
    % W11: LETFs -x2 return over the holding period t
    % W12: LETFs -x1 return over the holding period t
    % W20: LETFs x3 return over the holding period t
    % W21: LETFs -x3 return over the holding period t
    % D10: LETFs x2 return deviations over the holding period t
    % D11: LETFs -x2 return deviations over the holding period t
    % D12: LETFs -x1 return deviations over the holding period t
    % D20: LETFs x3 return deviations over the holding period t
    % D21: LETFs -x3 return deviations over the holding period t

    W= S(t,:)./S(1,:) - 1;
    W10= S10(t,:)./S10(1,:) - 1 ;
    W11= S11(t,:)./S11(1,:) - 1 ;
    W12= S12(t,:)./S12(1,:) - 1 ;
    W20= S20(t,:)./S20(1,:) - 1;
    W21= S21(t,:)./S21(1,:) - 1;
    D10=(- W*2 + W10);
    D11=(- W*-2 + W11);
    D12=(- W*-1 + W12);
    D20=(- W*3 + W20);
    D21=(- W*-3 + W21);

    mD10= std(D10);
    mD11= std(D11);
    mD12= std(D12);
    mD20= std(D20);
    mD21= std(D21);
    out= [mD21 mD11 mD12 mD10 mD20]; % Ordered according to leverage
end