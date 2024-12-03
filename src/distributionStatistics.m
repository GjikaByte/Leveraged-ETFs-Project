function [ avg, med, stdev, minimum, maximum, skew, kurt, excessKurt, ...
    nbPos, nbPercPos, lowerPerc, upperPerc, volatility ] ...
    = distributionStatistics( timeSeries, frequ, lowerPerc, upperPerc )
    
    % mean
    avg = mean(timeSeries);
    % median
    med = median(timeSeries);
    % standard deviation
    stdev = std(timeSeries);
    % minimum
    minimum = min(timeSeries);
    % maximum
    maximum = max(timeSeries);
    % skewness
    skew = skewness(timeSeries);
    % (excess) kurtosis
    kurt = kurtosis(timeSeries);
    excessKurt = kurtosis(timeSeries) - 3;
    % N > 0
    nbPos = sum(timeSeries > 0);
    % N > 0 in %
    nbPercPos = nbPos / size(timeSeries, 1);
    % lower percentile
    lowerPerc = prctile(timeSeries, lowerPerc);
    % upper percentile
    upperPerc = prctile(timeSeries, upperPerc);
    % annualized volatility
    volatility = sqrt(frequ) * stdev;
    
end