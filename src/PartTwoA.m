%% Add personal paths
addpath('/Users/kdallo/Documents/EDHEC/Courses/Master Project/MATLAB Code/Other Functions')
addpath('/Users/kdallo/Documents/EDHEC/Courses/Master Project/MATLAB Code/Code')
addpath('/Users/kdallo/Documents/EDHEC/Courses/Master Project/MATLAB Code/Data')
addpath('/Users/kdallo/Documents/EDHEC/Courses/Master Project/MATLAB Code/Data/Dax - 2X')
addpath('/Users/kdallo/Documents/EDHEC/Courses/Master Project/MATLAB Code/Figures')

%% Some constants and variables
% Set daily (d), weekly (w), monthly (m), quarterly (q) or yearly (y):
frequ = 'd';%           <--------- CHOOSE FREQUENCY
x = -2;         % Multiple (Leverage) = -2
figNb = 1;      % Counter 'i' for figures
nbBins = 100;   % Number of bins for distribution plots
lowerPerc = 10; % lower percentile
upperPerc = 90; % upper percentile
xtickerPerc = [0,0.0855, 0.217, 0.347, 0.478, 0.608, 0.7385, 0.87, 1];

%% Read data from Excel file (has to be csv format)
if frequ == 'd'
    [~, ~, rawData] = FNreadcsv('dailyPrices.csv');
    annualizeFrequ = 252;
    xtickers = xtickerPerc * 2000;
    titleFrequ = 'Daily';
elseif frequ == 'w'
    [~, ~, rawData] = FNreadcsv('weeklyPrices.csv');
    annualizeFrequ = 52;
    xtickers = xtickerPerc * 400;
    titleFrequ = 'Weekly';
elseif frequ == 'm'
    [~, ~, rawData] = FNreadcsv('monthlyPrices.csv');
    annualizeFrequ = 12;
    xtickers = xtickerPerc * 93;
    titleFrequ = 'Monthly';
elseif frequ == 'q'
    [~, ~, rawData] = FNreadcsv('quarterlyPrices.csv');
    annualizeFrequ = 4;
    xtickers = xtickerPerc * 32;
    titleFrequ = 'Quarterly';
elseif frequ == 'y'
    [~, ~, rawData] = FNreadcsv('yearlyPrices.csv');
    annualizeFrequ = 1;
    xtickers = xtickerPerc * 8;
    titleFrequ = 'Yearly';
end

% Historical prices for XSD2 and DAX from 10 May 2010 to 29 December 2017
datesRaw = rawData(2:end,1);
priceXSD2 = cell2mat(rawData(2:end,2));
priceDAX = cell2mat(rawData(2:end,3));

% Compute returns for XSD2, DAX and naive XSD2
returnXSD2 = priceXSD2(2:end) ./ priceXSD2(1:end-1) - 1;
returnDAX = priceDAX(2:end) ./ priceDAX(1:end-1) - 1;
naiveReturnXSD2 = returnDAX * x;

% Compute naive LETF path
naivePriceXSD2 = cumprod([priceXSD2(1); naiveReturnXSD2 + 1]);

%% Price chart for XSD2, naive XSD2 and DAX
figure(figNb)
figNb = figNb + 1;
plot([naivePriceXSD2, priceXSD2, priceDAX])
grid on
title(strcat(titleFrequ, ' Prices for DAX, XSD2 & Naive XSD2 Since Launch of the LETF'))
ylabel('Price')
legend({'naive XSD2','XSD2','DAX'},'Location','northwest', 'FontSize', 14)
xticks(xtickers)
xticklabels({'May','Dec 2010','Dec 2011','Dec 2012','Dec 2013','Dec 2014','Dec 2015','Dec 2016','Dec 2017'})

%% Price chart for XSD2 and naive XSD2
figure(figNb)
figNb = figNb + 1;
plot([naivePriceXSD2, priceXSD2])
grid on
title(strcat(titleFrequ, ' Prices for XSD2 & Naive XSD2 Since Launch of the LETF'))
ylabel('Price')
legend({'naive XSD2', 'XSD2'}, 'FontSize', 14)
xticks(xtickers)
xticklabels({'May','Dec 2010','Dec 2011','Dec 2012','Dec 2013','Dec 2014','Dec 2015','Dec 2016','Dec 2017'})

%% Plots for DAX, XSD2 and naive XSD2 return distribution
% XSD2
figure(figNb)
figNb = figNb + 1;
histfit(returnXSD2, nbBins, 'kernel');
grid on;
title(strcat(titleFrequ, ' Return Distribution for XSD2 Since Launch of the LETF'))
xlabel('Return');
% Naive
figure(figNb)
figNb = figNb + 1;
histfit(naiveReturnXSD2, nbBins, 'kernel');
grid on;
title(strcat(titleFrequ, ' Return Distribution for Naive XSD2 Since Launch of the LETF'))
xlabel('Return');
% DAX
figure(figNb)
figNb = figNb + 1;
histfit(returnDAX, nbBins, 'kernel');
grid on;
title(strcat(titleFrequ, ' Return Distribution for DAX Since Launch of the LETF'))
xlabel('Return');

%% Summary statistics for returns
% DAX
[ meanDAX, medianDAX, stdDAX, minDAX, maxDAX, skewDAX, kurtDAX, ...
    excessKurtDAX, nbPositivesDAX, nbPercPosDAX, tenPercDAX, ...
    nintyPercDAX, volDAX ] = distributionStatistics( returnDAX, ...
    annualizeFrequ, lowerPerc, upperPerc );
% XSD2
[ meanXSD2, medianXSD2, stdXSD2, minXSD2, maxXSD2, skewXSD2, kurtXSD2, ...
    excessKurtXSD2, nbPositivesXSD2, nbPercPosXSD2, tenPercXSD2, ...
    nintyPercXSD2, volXSD2 ] = distributionStatistics( returnXSD2, ...
    annualizeFrequ, lowerPerc, upperPerc );
% Naive XSD2
[ meanNaiveXSD2, medianNaiveXSD2, stdNaiveXSD2, minNaiveXSD2, ...
    maxNaiveXSD2, skewNaiveXSD2, kurtNaiveXSD2, excessKurtNaiveXSD2, ...
    nbPositivesNaiveXSD2, nbPercPosNaiveXSD2, tenPercNaiveXSD2, ...
    nintyPercNaiveXSD2, volNaiveXSD2 ] = distributionStatistics( ...
    naiveReturnXSD2, annualizeFrequ, lowerPerc, upperPerc );

%% Return ratio (XSD2 return / naive XSD2 return)
returnRatio = returnXSD2 ./ naiveReturnXSD2;
% Exclude extreme outliers - !FOR PLOT ONLY!
returnRatioClean = returnRatio;
nbDeleted = size(returnRatioClean(abs(returnRatioClean) > 5), 1);
percDeleted = nbDeleted / size(returnRatioClean, 1);
returnRatioClean(abs(returnRatioClean) > 5) = [];
% Plot distribution
figure(figNb)
figNb = figNb + 1;
histfit(returnRatioClean, 80, 'kernel');
grid on
title(strcat(titleFrequ, ' Return Ratio (XSD2 Return / Naive XSD2 Return) Distribution'))
xlabel('Return Ratio');
% Plot return ratio over time
figure(figNb)
figNb = figNb + 1;
plot(returnRatioClean);
grid on
title(strcat(titleFrequ, ' Return Ratio (XSD2 Return / Naive XSD2 Return)'))
ylabel('Return Ratio')
xticks(xtickers)
xticklabels({'May','Dec 2010','Dec 2011','Dec 2012','Dec 2013','Dec 2014','Dec 2015','Dec 2016','Dec 2017'})
% Compute distribution statistics
[ meanRetRat, medRetRat, stdevRetRat, minimumRetRat, maximuRetRat, ...
    skewRetRat, kurtRetRat, excessKurtRetRat, nbPosRetRat, ...
    nbPercPosRetRat, lowerPercRetRat, upperPercRetRat, volatilityRetRat ...
    ] = distributionStatistics( returnRatioClean, annualizeFrequ, ...
    lowerPerc, upperPerc );

%% Effective leverage (XSD2 return / benchmark return)
timeSeriesEffLev = returnXSD2 ./ returnDAX;
% Compute distribution statistics
[ meanEffLev, medEffLev, stdevEffLev, minimumEffLev, maximuEffLev, ...
    skewEffLev, kurtEffLev, excessKurtEffLev, nbPosEffLev, ...
    nbPercPosEffLev, lowerPercEffLev, upperPercEffLev, volatilityEffLev...
    ] = distributionStatistics( timeSeriesEffLev, annualizeFrequ, ...
    lowerPerc, upperPerc );
% Exclude extreme outliers - !FOR PLOT ONLY!
timeSeriesEffLevClean = timeSeriesEffLev;
nbDeletedEffLev = size(timeSeriesEffLevClean(abs(timeSeriesEffLevClean) > 10), 1);
percDeletedEffLev = nbDeletedEffLev / size(timeSeriesEffLevClean, 1);
timeSeriesEffLevClean(abs(timeSeriesEffLevClean) > 10) = [];
% Plot
figure(figNb)
figNb = figNb + 1;
histfit(timeSeriesEffLevClean,80,'kernel')
grid on
title(strcat(titleFrequ, ' Effective Leverage (XSD2 Return / DAX Return) Distribution'))
xlabel('Return Ratio');
figure(figNb)
figNb = figNb + 1;
plot(timeSeriesEffLevClean)
grid on
title(strcat(titleFrequ, ' Effective Leverage (XSD2 Return / DAX Return)'))
ylabel('Effective Leverage')
xticks(xtickers)
xticklabels({'May','Dec 2010','Dec 2011','Dec 2012','Dec 2013','Dec 2014','Dec 2015','Dec 2016','Dec 2017'})

%% Return deviation
returnDeviation = returnXSD2 - naiveReturnXSD2;
% Summary statistics for return deviation
[ meanDeviation, medianDeviation, stdDeviation, minDeviation, ...
    maxDeviation, skewDeviation, kurtDeviation, excessKurtDeviation, ...
    nbPositivesDeviation, nbPercPosDeviation, tenPercDeviation, ...
    nintyPercDeviation, volDeviation ] = distributionStatistics( ...
    returnDeviation, annualizeFrequ, lowerPerc, upperPerc );
% Plot
figure(figNb)
figNb = figNb + 1;
histfit(returnDeviation, nbBins, 'kernel');
grid on
title(strcat(titleFrequ, ' Return Deviation Distribution'))
xlabel('Return Deviation (Target - Naive)')

%% Random generator of average return deviation per holding period
% If chosen frequency is not daily, read daily dataset. Otherwise use array
% we created in the beginning of this code.
if frequ ~= 'd'
    [~, ~, rawData] = FNreadcsv('dailyPrices.csv');
    % Historical prices for XSD2 and DAX from 10 May 2010 to 29 December 2017
    datesRawDaily = rawData(2:end,1);
    priceXSD2Daily = cell2mat(rawData(2:end,2));
    priceDAXDaily = cell2mat(rawData(2:end,3));
    % Compute returns for XSD2, DAX and naive XSD2
    returnXSD2Daily = priceXSD2Daily(2:end) ./ priceXSD2Daily(1:end-1) - 1;
    returnDAXDaily = priceDAXDaily(2:end) ./ priceDAXDaily(1:end-1) - 1;
    naiveReturnXSD2Daily = returnDAXDaily * x;
    naivePriceXSD2Daily = cumprod([priceXSD2Daily(1); naiveReturnXSD2Daily + 1]);
else
    returnXSD2Daily = returnXSD2;
    returnDAXDaily = returnDAX;
    priceXSD2Daily = priceXSD2;
    priceDAXDaily = priceDAX;
    naiveReturnXSD2Daily = naiveReturnXSD2;
    naivePriceXSD2Daily = naivePriceXSD2;
end

%% Overlapping mean return deviation and effective leverage
% Declaration
days = 1:1:150;
averageRetDev = zeros(size(days,2),1);
averageRetDevAbs = zeros(size(days,2),1);
averageEffLev = zeros(size(days,2),1);
avgLowPerc = zeros(size(days,2),1);
avgUppPerc = zeros(size(days,2),1);
% Go through each period given in the vector 'days' as many times as
% possible with our time series (e.g. for holding period = 1 we can cover
% 1,912 holding periods, for holding period = 252 days we can cover 1,660
% holding periods.
for j=1:size(days, 2)
    retDev = 0; retDevAbs = 0; tmpRetDev = 0; effLev = 0; tmpEffLev = 0;
    increment = days(j);
    % Vector to store return deviation in order to compute lower percentile
    timeSeriesRetDev = zeros(size(priceXSD2Daily,1) - increment,1);
    for i=1:size(priceXSD2Daily,1) - increment
        % Compute returns and sum them up
        tmpRetDev = ( (priceXSD2Daily(i + increment) / ...
            priceXSD2Daily(i) - 1) - ...
            (naivePriceXSD2Daily(i + increment) / ...
            naivePriceXSD2Daily(i) - 1) );
        retDev = retDev + tmpRetDev;
        retDevAbs = retDevAbs + abs(tmpRetDev);
        % Compute effective leverage
        tmpEffLev = ( (priceXSD2Daily(i + increment) / ...
            priceXSD2Daily(i) - 1) / (priceDAXDaily(i + increment) / ...
            priceDAXDaily(i) - 1) );
        if isinf(tmpEffLev)
            % if return of DAX is 0, we have to adjust since we cannot
            % divide by 0
            tmpEffLev = x;
        end
        effLev = effLev + tmpEffLev;
        % Compute lower percentile
        timeSeriesRetDev(i) = tmpRetDev;
    end
    % Take average
    averageRetDev(j) = (100 * retDev / i);
    averageRetDevAbs(j) = (retDevAbs / i);
    averageEffLev(j) = (effLev / i);
    [~,~,~,~,~,~,~,~,~,~,avgLowPerc(j),avgUppPerc(j),~] = distributionStatistics( ...
        timeSeriesRetDev, annualizeFrequ, lowerPerc, upperPerc );
    avgLowPerc(j) = ( avgLowPerc(j) * 100 ) / i;
    avgUppPerc(j) = ( avgUppPerc(j) * 100 ) / i;
end

% Plot average return deviation
% Average return deviation
figure(figNb)
figNb = figNb + 1;
plot(days, averageRetDev)
grid on
title('Average Return Deviation as a Function of the Holding Period')
ylabel('Average Return Deviation (Target - Naive) in %')
xlabel('Holding Period in Days')
% Average absolute return deviation
figure(figNb)
figNb = figNb + 1;
plot(days, averageRetDevAbs)
grid on
title('Average Absolute Return Deviation as a Function of the Holding Period')
ylabel('Average Absolute Return Deviation |Target - Naive|')
xlabel('Holding Period in Days')
% Average effective leverage
figure(figNb)
figNb = figNb + 1;
plot(days, averageEffLev)
grid on
title('Average Effective Leverage as a Function of the Holding Period')
ylabel('Average Effective Leverage |Target / Benchmark|')
xlabel('Holding Period in Days')
% Average lower percentile
figure(figNb)
figNb = figNb + 1;
plot(days,avgLowPerc)
%plot(days(1:600), avgLowPerc(1:600))
grid on
title({'Average 10% Percentile of the Return Deviation','as a Function of the Holding Period'})
ylabel('Average 10% Percentile')
xlabel('Holding Period in Days')
% Average upper percentile
figure(figNb)
figNb = figNb + 1;
plot(days, avgUppPerc)
grid on
title({'Average 90% Percentile of the Return Deviation','as a Function of the Holding Period'})
ylabel('Average 90% Percentile')
xlabel('Holding Period in Days')

%% Additional plots
% The data for these plots is coming from the code above
category = categorical(["daily" "weekly" "monthly" "quarterly" "yearly"]);
category = reordercats(category,{'daily' 'weekly' 'monthly' 'quarterly' 'yearly'});

% Mean and median of return deviation per holding period
rdMeanMedian = [ -0.00195122150757669 -0.0232393117205742; ...
    -0.0129461281342646 -0.00980723872361366; ...
    -0.2994336342329 -0.4685896102131; ...
    -0.2589373127978 0.366267039314; ...
    -4.64042151002 0.3854182637761 ];
figure(figNb)
figNb = figNb + 1;
bar(category, rdMeanMedian)
grid on
title({'Mean and Median per Holding Period'})
legend({'Mean','Median'},'Location','northwest', 'FontSize', 14)
ylabel('%')
xlabel('Holding Period')

% 10 and 90 percentiles of return deviation per holding period
rdPercentiles = [ -0.7493643236195 0.7258468746179; ...
    -1.6046075018999 1.5511248945856; ...
    -3.662599671152 2.8017322853041; ...
    -5.1610462486085 3.29981991952; ...
    -23.2824939317325 10.4140930690162 ];
figure(figNb)
figNb = figNb + 1;
bar(category, rdPercentiles)
grid on
title({'10 and 90 Percentiles per Holding Period'})
legend({'10th Percentile','90th Percentile'},'Location','northwest', 'FontSize', 14)
ylabel('%')
xlabel('Holding Period')

% Minimum and maximum return deviation
rdMinMax = [ -3.1937482610469 7.4504270908775; ...
    -4.6643772632657 5.1708217518666; ...
    -7.0841770218604 8.8875598878199; ...
    -8.3997598578196 4.491510344378; ...
    -25.346555947711 10.7035411067942 ];
figure(figNb)
figNb = figNb + 1;
bar(category, rdMinMax)
grid on
title({'Minimum and Maximum per Holding Period'})
legend({'Minimum','Maximum'},'Location','northwest', 'FontSize', 14)
ylabel('%')
xlabel('Holding Period')