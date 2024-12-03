%% SECOND PART
%% Preliminary Code

%We set the relevant paths
close all
clc 
cd('/Users/andigjika/Desktop/Matlab functions and plots');
addpath('/Users/andigjika/Desktop/Matlab functions and plots');
addpath('/Users/andigjika/Desktop/Matlab functions and plots/Functions');
addpath ('/Users/andigjika/Dropbox/Master Project Thesis/Data');

% The following code extracts the last year prices from the raw prices data
% and calculates the returns
% We extract data from daxcopy
% daxcopy: .csv file with prices of the DAX Index from 30/12/2009 to 17/01/2018
% The function'FNreadcsv' uploads data from .csv files

% We calculate the following values:
% dataAll: Prices of the DAX over 2045 trading days
% dataRaw2: Prices of the DAX over 1 year (252 trading days)
% retAll: Returns of the DAX over 2045 trading days
% ret: Returns of the DAX over 1 year (252 trading days)

% We upload the dax prices
[~,~,dataRaw2]=FNreadcsv('daxcopy.csv');

% We manipulate the data to be able to work on it
Data(1,1) = str2double(dataRaw2{1, 1});
Data(1,2) = str2double(dataRaw2{1, 2});
for i = 2 : size(dataRaw2,1)
    Data(i,1) = str2double(dataRaw2{i, 1});
    Data(i,2) = dataRaw2{i, 2};
end
Dat(:,1) = Data(:,1) + Data(:,2) ./ 100;
clear dataRaw2;

% We select prices and returns of the last year
dataAll= Dat;
dataRaw2= Dat(1:252);
dataAll=flip(dataAll);
dataRaw2=flip(dataRaw2);
ret=PriceToReturn(dataRaw2);
retAll=PriceToReturn(dataAll);

% We clear the stack
clear Dat Data i

%% 2.1 Return Deviation Graphs || Distributions of Return deviations
% In this subsection we do the following:
% 1) Generate 5,000 price paths of the DAX over 5 years
% 2) Elaborate the resulting LETFs -x3, -x2, -x1, x2, x3 prices over the period
% 3) Calculate the return of the LETFs for holding periods of one week, one month,
%    one year
% 4) For each holdig period, we calculate the LETFs returns over the 5,000
%    simulations
% 5) For each holding period, we calculate the return deviation of the
%    LETFs with naive expected return
% 6) For each holding period, we elaborate and plot the distributions
%    of the return deviations

% NOTATION:
% S: 5,000 Simulated DAX prices over the last 5 years
% S10: LETFs x2 on S
% S11: LETFs -x2 on S
% S12: LETFs -x1 on S
% S20: LETFs x3 on S
% S21: LETFs -x3 on S
% W: Return of the 5,000 DAX price simulations over the holding period
% W10: LETFs x2 return over the holding period
% W11: LETFs -x2 return over the holding period
% W12: LETFs -x1 return over the holding period
% W20: LETFs x3 return over the holding period
% W21: LETFs -x3 return over the holding period

% 1) + 2) 
vol5=std(retAll(end-1260:end))*sqrt(1260);
mi5=(dataAll(end)/dataAll(end-1260))^(1/1260)-1;
[S]= IndexPrice_5years (100, mi5, vol5, 5000, 1260, 1260, 0, 0, 0);
S10=LETFprice (S, 2, S(1)); S11=LETFprice(S,-2,S(1)); S12=LETFprice(S,-1,S(1)); S20=LETFprice(S,3,S(1)); S21=LETFprice(S,-3,S(1));

% HOLDING PERIOD: 1 Week || 3) + 4) +5) + 6)
% We use the returns of the first week that we simulated
% NOTATION:
% D10: LETFs x2 return deviations over 1 week holding period
% D11: LETFs -x2 return deviations over 1 week holding period
% D12: LETFs -x1 return deviations over 1 week holding period
% D20: LETFs x3 return deviations over 1 week holding period
% D21: LETFs -x3 return deviations over 1 week holding period
% A: Table with relevant statistics of the return deviations' distribution
W= S(5,:)./S(1,:) - 1;
W10= S10(5,:)./S10(1,:) - 1 ;
W11= S11(5,:)./S11(1,:) - 1 ;
W12= S12(5,:)./S12(1,:) - 1 ;
W20= S20(5,:)./S20(1,:) - 1;
W21= S21(5,:)./S21(1,:) - 1;
D10=(- W*2 + W10);
D11=(- W*-2 + W11);
D12=(- W*-1 + W12);
D20=(- W*3 + W20);
D21=(- W*-3 + W21);
A=[mean(D10) prctile(D10,10) std(D10) skewness(D10) kurtosis(D10);
mean(D11) prctile(D11,10) std(D11) skewness(D11) kurtosis(D11);
mean(D12) prctile(D12,10) std(D12) skewness(D12) kurtosis(D12);
mean(D20) prctile(D20,10) std(D20) skewness(D20) kurtosis(D20);
mean(D21) prctile(D21,10) std(D21) skewness(D21) kurtosis(D21)];

% We plot the table of the return distributions for one week holding period
[D10y, D10x] = hist(D10,300);
[D11y, D11x] = hist(D11,300);
[D12y, D12x] = hist(D12,300);
[D20y, D20x] = hist(D20,300);
[D21y, D21x] = hist(D21,300);
plot(D10x,D10y), hold on, plot(D11x,D11y),plot(D12x,D12y), plot(D20x,D20y), plot(D21x,D21y), legend()
title('Return deviation distribution over one week holding period')
ylabel('Frequency')
xlabel('Return Deviation')
legend({'LETF x2', 'LETF -x2','LETF -x1', 'LETF x3','LETF -x3'},'Location','northeast','FontSize',14)

% HOLDING PERIOD: 1 Month || 3) + 4) +5) + 6)
% We use the returns of the first month that we simulated
% NOTATION:
% D101: LETFs x2 return deviations over 1 month holding period
% D111: LETFs -x2 return deviations over 1 month holding period
% D121: LETFs -x1 return deviations over 1 month holding period
% D201: LETFs x3 return deviations over 1 month holding period
% D211: LETFs -x3 return deviations over 1 month holding period
% NOTE THAT we overwrite W, W10, ... and A
W= S(21,:)./S(1,:) - 1;
W10= S10(21,:)./S10(1,:) - 1 ;
W11= S11(21,:)./S11(1,:) - 1 ;
W12= S12(21,:)./S12(1,:) - 1 ;
W20= S20(21,:)./S20(1,:) - 1;
W21= S21(21,:)./S21(1,:) - 1;
D101=(- W*2 + W10);
D111=(- W*-2 + W11);
D121=(- W*-1 + W12);
D201=(- W*3 + W20);
D211=(- W*-3 + W21);
A=[mean(D101) prctile(D101,10) std(D101) skewness(D101) kurtosis(D101);
mean(D111) prctile(D111,10) std(D111) skewness(D111) kurtosis(D111);
mean(D121) prctile(D121,10) std(D121) skewness(D121) kurtosis(D121);
mean(D201) prctile(D201,10) std(D201) skewness(D201) kurtosis(D201);
mean(D211) prctile(D211,10) std(D211) skewness(D211) kurtosis(D211)];

% We plot the table of the return distributions for one month holding period
[D101y, D101x] = hist(D101,300);
[D111y, D111x] = hist(D111,300);
[D121y, D121x] = hist(D121,300);
[D201y, D201x] = hist(D201,300);
[D211y, D211x] = hist(D211,300);
plot(D101x,D101y), hold on, plot(D111x,D111y), plot(D121x,D121y), plot(D201x,D201y), plot(D211x,D211y), 
title('Return deviation distribution over one month holding period')
ylabel('Frequency')
xlabel('Return Deviation')
legend({'LETF x2', 'LETF -x2','LETF -x1', 'LETF x3','LETF -x3'},'Location','northeast','FontSize',14)

% HOLDING PERIOD: 1 Year || 3) + 4) +5) + 6)
% We use the returns of the first month that we simulated
% NOTATION:
% D102: LETFs x2 return deviations over 1 year holding period
% D112: LETFs -x2 return deviations over 1 year holding period
% D122: LETFs -x1 return deviations over 1 year holding period
% D202: LETFs x3 return deviations over 1 year holding period
% D212: LETFs -x3 return deviations over 1 year holding period
% NOTE THAT we overwrite W, W10, ... and A
W= S(252,:)./S(1,:) - 1;
W10= S10(252,:)./S10(1,:) - 1 ;
W11= S11(252,:)./S11(1,:) - 1 ;
W12= S12(252,:)./S12(1,:) - 1 ;
W20= S20(252,:)./S20(1,:) - 1;
W21= S21(252,:)./S21(1,:) - 1;
D102=(- W*2 + W10);
D112=(- W*-2 + W11);
D122=(- W*-1 + W12);
D202=(- W*3 + W20);
D212=(- W*-3 + W21);
A=[mean(D102) prctile(D102,10) std(D102) skewness(D102) kurtosis(D102);
mean(D112) prctile(D112,10) std(D112) skewness(D112) kurtosis(D112);
mean(D122) prctile(D122,10) std(D122) skewness(D122) kurtosis(D122);
mean(D202) prctile(D202,10) std(D202) skewness(D202) kurtosis(D202);
mean(D212) prctile(D212,10) std(D212) skewness(D212) kurtosis(D212)];

% We plot the table of the return distributions for one year holding period
[D102y, D102x] = hist(D102,300);
[D112y, D112x] = hist(D112,300);
[D122y, D122x] = hist(D122,300);
[D202y, D202x] = hist(D202,300);
[D212y, D212x] = hist(D212,300);
plot(D102x,D102y), hold on, plot(D112x,D112y), plot(D122x,D122y), plot(D202x,D202y), plot(D212x,D212y)
title('Return deviation distribution over one month holding period')
ylabel('Frequency')
xlabel('Return Deviation')
legend({'LETF x2', 'LETF -x2','LETF -x1', 'LETF x3','LETF -x3'},'Location','northeast','FontSize',14)

% We clear the stack
clear D10 D10x D10y D101 D101x D101y D11 D11x D11y D111 D111x D111y
clear D12 D12x D12y D121 D121x D121y D102 D102x D102y D112 D112x D112y
clear D122 D122x D122y D20 D20x D20y D201 D201x D201y D202 D202x D202y
clear D21 D21x D21y D211 D211x D211y D212 D212x D212y W W10 W11 W12 W20 W21

%% 2.1 Return Deviation Graphs || Mean of Return Deviation
% In this subsection we generate a plot of the mean of the return deviation
% against leverage and holding period
% NOTE THAT we use data from our previous simulations
% NOTATION:
% A: Vector of mean of the return deviations
A = zeros(5,25);
for i = 1 : 25
    A(:,i) = meanRetDev(S, S10, S11, S12, S20, S21, i*10);
end
surf([1:25],[1:5],A)
title('Mean of Return Deviation')
ylabel('Leverage')
xlabel('Holding days')
yticklabels({'-x3','-x2','-x1','x2','x3'});
xticklabels({'0','50','100','150','200','250'});
zlabel('Mean Deviation')

% We clear the stack
clear i A
%% 2.1 Return Deviation Graphs || Standard Deviation of Return Deviation
% In this subsection we generate a plot of the standard deviation of the return deviation
% against leverage and holding period
% NOTE THAT we use data from our previous simulations
% NOTATION:
% A: Vector of standard eviation of the return deviations
A = zeros(5,25);
for i = 1 : 25
    A(:,i) = stdRetDev(S, S10, S11, S12, S20, S21, i*10);
end
surf([1:25],[1:5],A)
title('Standard Deviation of Return Deviation')
ylabel('Leverage')
xlabel('Holding days')
yticklabels({'-x3','-x2','-x1','x2','x3'});
xticklabels({'0','50','100','150','200','250'});
zlabel('Standard Deviation')

% We clear the stack
clear i A
%% 2.2 Compounding Effect || 10% Percentiles, different K
% In this subsection we calculate the graph of the 10% percentile of the
% return deviation distribution against leverage and holding days

% We are interested in inferring mrfj, sigmaj and lambda
% We assume the 5 year volatility and daily geometric mean return
% We estimate mrfj, sigmaj and lambda by identifying excess returns
% We use different levels of daily returns volatility
% NOTATION:
% vol: daily volatility
% mi: equal to mi5, geometric mean of returns
% lambda: daily number of excess returns
% excess: vector of excess returns
% mij: mean of excess returns
% volj: standard deviation of excess returns
% S: 5,000 Simulated DAX prices over the last 5 years WITH JUMPS
% S10: LETFs x2 on S
% S11: LETFs -x2 on S
% S12: LETFs -x1 on S
% S20: LETFs x3 on S
% S21: LETFs -x3 on S

vol=std(retAll(end-1260:end));
vol5= vol * sqrt(1260);
mi=(dataAll(end)/dataAll(end-1260))^(1/1260)-1;

% x1 Standard Deviation: 68% || K = 1
lambda = 0;
excess= 0;
j=0;
for i= 1 : 1260
    if retAll(end - i) > mi + vol || retAll(end - i) < mi - vol
        lambda= lambda+1;
        j= j+1;
        excess(j) = retAll(end-i);
    end
end
lambda = lambda / 1260;
mij= mean(excess);
volj= std(excess);

[S]= IndexPrice_5years (100, mi, vol5, 5000, 1260, 1260, lambda, mij, volj);
S10=LETFprice (S, 2, S(1)); S11=LETFprice(S,-2,S(1)); S12=LETFprice(S,-1,S(1)); S20=LETFprice(S,3,S(1)); S21=LETFprice(S,-3,S(1));

% We generate the median/mean of return deviations graph
A = zeros(5,25);
for i = 1 : 25
    A(:,i) = PercentileRetDev(S, S10, S11, S12, S20, S21, i*10,10);
end
surf([1:25],[1:5],A)
title('Excess level at 68%')
ylabel('Leverage')
xlabel('Holding days')
yticklabels({'-x3','-x2','-x1','x2','x3'});
xticklabels({'0','50','100','150','200','250'});
zlabel('10% Percentile')

% We clear the stack
clear i j lambda A excess

% x2 Standard Deviation: 95% || K=2
lambda = 0;
excess= 0;
j=0;
for i= 1 : 1260
    if retAll(end - i) > mi + 2*vol || retAll(end - i) < mi - 2*vol
        lambda= lambda+1;
        j= j+1;
        excess(j) = retAll(end-i);
    end
end
lambda = lambda / 1260;
mij= mean(excess);
volj= std(excess);
[S]= IndexPrice_5years (100, mi, vol5, 5000, 1260, 1260, lambda, mij, volj);
S10=LETFprice (S, 2, S(1)); S11=LETFprice(S,-2,S(1)); S12=LETFprice(S,-1,S(1)); S20=LETFprice(S,3,S(1)); S21=LETFprice(S,-3,S(1));

% We generate the median/mean of return deviations graph
A = zeros(5,25);
for i = 1 : 25
    A(:,i) = PercentileRetDev(S, S10, S11, S12, S20, S21, i*10,10);
end
surf([1:25],[1:5],A)
title('Excess level at 95%')
ylabel('Leverage')
xlabel('Holding days')
yticklabels({'-x3','-x2','-x1','x2','x3'});
xticklabels({'0','50','100','150','200','250'});
zlabel('10% Percentile')

% We clear the stack
clear i j lambda A excess

% 3 vol: 99.7% || K=3
lambda = 0;
excess= 0;
j=0;
for i= 1 : 1260
    if retAll(end - i) > mi + 3*vol || retAll(end - i) < mi - 3*vol
        lambda= lambda+1;
        j= j+1;
        excess(j) = retAll(end-i);
    end
end
lambda = lambda / 1260;
mij= mean(excess);
volj= std(excess);

[S]= IndexPrice_5years (100, mi, vol5, 5000, 1260, 1260, lambda, mij, volj);
S10=LETFprice (S, 2, S(1)); S11=LETFprice(S,-2,S(1)); S12=LETFprice(S,-1,S(1)); S20=LETFprice(S,3,S(1)); S21=LETFprice(S,-3,S(1));

% We generate the median/mean of return deviations graph
A = zeros(5,25);
for i = 1 : 25
    A(:,i) = PercentileRetDev(S, S10, S11, S12, S20, S21, i*10,10);
end
surf([1:25],[1:5],A)
title('Mean of Return Deviations')
ylabel('Leverage')
xlabel('Holding days')
yticklabels({'-x3','-x2','-x1','x2','x3'});
xticklabels({'0','50','100','150','200','250'});
zlabel('Mean Deviation')

% We clear the stack
clear i j lambda A excess

%% 2.2 Compounding Effect || 10% Percentiles, K=2, different volatility
% In this subsection we calculate the graph of the 10% percentile of the
% return deviation distribution against leverage and holding days

% We assume the 5 year volatility and daily geometric mean return
% We estimate mrfj, sigmaj and lambda by using 2 Vol excess returns
% NOTATION:
% lambda: daily number of excess returns
% excess: vector of excess returns
% mij: mean of excess returns
% volj: standard deviation of excess returns
% S: 5,000 Simulated DAX prices over the last 5 years WITH JUMPS
% S10: LETFs x2 on S
% S11: LETFs -x2 on S
% S12: LETFs -x1 on S
% S20: LETFs x3 on S
% S21: LETFs -x3 on S

% We set volatility to 50% of the historical value
vol5= .5 * vol5;
lambda = 0;
excess= 0;
j=0;
for i= 1 : 1260
    if retAll(end - i) > mi + 2*vol || retAll(end - i) < mi - 2*vol
        lambda= lambda+1;
        j= j+1;
        excess(j) = retAll(end-i);
    end
end
lambda = lambda / 1260;
mij= mean(excess);
volj= std(excess);
[S]= IndexPrice_5years (100, mi, vol5, 5000, 1260, 1260, lambda, mij, volj);
S10=LETFprice (S, 2, S(1)); S11=LETFprice(S,-2,S(1)); S12=LETFprice(S,-1,S(1)); S20=LETFprice(S,3,S(1)); S21=LETFprice(S,-3,S(1));

% We generate the 10% percentile of return deviations graph
A = zeros(5,25);
for i = 1 : 25
    A(:,i) = PercentileRetDev(S, S10, S11, S12, S20, S21, i*10,10);
end
surf([1:25],[1:5],A)
title('50% Volatility')
ylabel('Leverage')
xlabel('Holding days')
yticklabels({'-x3','-x2','-x1','x2','x3'});
xticklabels({'0','50','100','150','200','250'});
zlabel('10% Percentile')
vol5 = vol5 * 2;

% We clear the stack
clear i j lambda A excess

% We set volatility to 200% of the historical value
vol5= 2 * vol5;
lambda = 0;
excess= 0;
j=0;
for i= 1 : 1260
    if retAll(end - i) > mi + 2*vol || retAll(end - i) < mi - 2*vol
        lambda= lambda+1;
        j= j+1;
        excess(j) = retAll(end-i);
    end
end
lambda = lambda / 1260;
mij= mean(excess);
volj= std(excess);
[S]= IndexPrice_5years (100, mi, vol5, 5000, 1260, 1260, lambda, mij, volj);
S10=LETFprice (S, 2, S(1)); S11=LETFprice(S,-2,S(1)); S12=LETFprice(S,-1,S(1)); S20=LETFprice(S,3,S(1)); S21=LETFprice(S,-3,S(1));

% We generate the median/mean of return deviations graph
A = zeros(5,25);
for i = 1 : 25
    A(:,i) = PercentileRetDev(S, S10, S11, S12, S20, S21, i*10,10);
end
surf([1:25],[1:5],A)
title('200% Volatility')
ylabel('Leverage')
xlabel('Holding days')
yticklabels({'-x3','-x2','-x1','x2','x3'});
xticklabels({'0','50','100','150','200','250'});
zlabel('10% Percentile')
vol5 = vol5 * .5;

% We clear the stack
clear i j lambda A excess
clear dataAll dataRaw2 ret retAll
clear mi mi5 mij vol vol5 volj
clear S S10 S11 S12 S20 S21
%% THIRD PART
%% Preliminary Code

% In this subsection we upload and work on the data for section 4
% We use data from 09/04/2010 to 12/04/2018
% NOTATION:
% dataAll: Prices of the DAX Index
% dataAllx2: Prices of the XSD2
% retAll: Prices of the DAX Index
% retAllx2: Prices of the XSD2
 
addpath ('/Users/andigjika/Dropbox/Master Project Thesis/Data/part3');

% We extract last year prices from the raw prices data and calculates the returns

% We upload the Dax Index prices
% daxeff2 contains the prices of the DAX Index
[~,~,dataRaw2]=FNreadcsv('daxeff2.csv');

% We manipulate the data to be able to work on it
Data(1,1) = str2double(dataRaw2{1, 1});
Data(1,2) = str2double(dataRaw2{1, 2});
for i = 2 : size(dataRaw2,1)
    Data(i,1) = str2double(dataRaw2{i, 1});
    Data(i,2) = dataRaw2{i, 2};
end
Dat(:,1) = Data(:,1) + Data(:,2) ./ 100;
dataAll= Dat;
dataAll=flip(dataAll);
retAll=PriceToReturn(dataAll);

% We clear the stack
clear dataRaw2;
clear Dat;
clear Data;
clear i;

% We upload the XSD2 prices
% SHRTDAX2eff2 contains the prices of the XSD2
[~,~,dataRaw2]=FNreadcsv('SHRTDAX2eff2.csv');

% We manipulate the data to be able to work on it
Data(1,1) = str2double(dataRaw2{1, 1});
Data(1,2) = str2double(dataRaw2{1, 2});
for i = 2 : size(dataRaw2,1)
    Data(i,1) = str2double(dataRaw2{i, 1});
    Data(i,2) = dataRaw2{i, 2};
end
Dat(:,1) = Data(:,1) + Data(:,2) ./ 100;
dataAllx2= Dat;
dataAllx2=flip(dataAllx2);
retAllx2=PriceToReturn(dataAllx2);

%We clear the stack
clear dataRaw2;
clear Dat;
clear Data;
clear i;

%% 4.1 Holding periods and Gamma
% In this subsection we calculate Gamma over one week and one month holding
% period multiple times
% We then collect relevant statistics for the distributions of gamma on the
% two holding periods (weekly and monthly)

% 1 Month
% NOTATION:
% N: Number of months in the sample
% i: month count
% j: day-count
% r: Realized returns over the holding period
% RLETF: Theoretical returns
% StatisticsM: Relevant statistics for monthly holding period
N = floor(size(retAll,1)/21);
for i = 1 : N
    r = retAllx2(21*(i-1) +1) + 1; % For each month, the first return on the first day
    for j = 2 : 20
        r = (1 + retAllx2(21*(i-1) +j) ) * r;  % Update r to the cumulative return
    end
    R = PriceToReturn(dataAll( (1 + (i-1) * 21) : (21 + (i-1) * 21 ) ));
    RLETF = LETFpriceEffective(dataAll,-2,dataAll((i-1)*21 + 1),R, 21);
    RLETF = RLETF/dataAll((i-1)*21 + 1);
    
    monthDev(i) = r / RLETF;
end 
 MeanM = mean(monthDev);
 VolM = std(monthDev);
 SkwM = skewness(monthDev);
 KrtM = kurtosis(monthDev);
 PrcM = prctile(monthDev,10);
 StatisticsM = [MeanM VolM SkwM KrtM PrcM];
 
 %We clear the stack
 clear MeanM VolM SkwM KrtM PrcM
 clear i j N r R RLETF VolM
 
 
% 1 Week
% NOTATION:
% N: Number of weeks in the sample
% i: week count
% j: day-count
% r: Realized returns over the holding period
% RLETF: Theoretical returns
% StatisticsW: Relevant statistics for weekly holding period
N = floor(size(retAll,1)/5);
for i = 1 : N
    r = retAllx2(5*(i-1) +1) + 1; % For each week, the first return on the first day
    for j = 2 : 4
        r = (1 + retAllx2(5*(i-1) +j) ) * r;  % Update r to the cumulative return
    end
    R = PriceToReturn(dataAll( (1 + (i-1) * 5) : (5 + (i-1) * 5 ) ));
    RLETF = LETFpriceEffective(dataAll,-2,dataAll((i-1)*5 + 1),R, 5);
    RLETF = RLETF/dataAll((i-1)*5 + 1);
    
    weekDev(i) = r / RLETF;
end 
 MeanW = mean(weekDev);
 VolW = std(weekDev);
 SkwW = skewness(weekDev);
 KrtW = kurtosis(weekDev);
 PrcW = prctile(weekDev,10);
 StatisticsW = [MeanW VolW SkwW KrtW PrcW];
 
% We clear the stack
 clear MeanW VolW SkwW KrtW PrcW
 clear i j N r R RLETF VolW weekDev monthDev
 
 %% 4.2 Gamma in Bull & Bear markets
% In this subsection we calculate Gamma over one week holding period
% in bull and bear markets

% NOTATION:
% N: Number of months in the sample
% i: month count
% j: day-count
% r: Realized returns over the holding period
% RLETF: Theoretical returns
% StatisticsMUP: Relevant statistics for bull markets
% StatisticsMDW: Relevant statistics for bear markets
N = floor(size(retAll,1)/5); % Number of months
    k=1;
    p=1;
for i = 1 : N
    r = retAllx2(5*(i-1) +1) + 1; % For each month, the first return on the first day
    rp = retAll(5*(i-1) +1) + 1;
    for j = 2 : 4
        r = (1 + retAllx2(5*(i-1) +j) ) * r; % Update r to the cumulative return
        rp = (1 + retAll(5*(i-1) +j) ) * rp;
    end
    R = PriceToReturn(dataAll( (1 + (i-1) * 5) : (5 + (i-1) * 5 ) ));
    RLETF = LETFpriceEffective(dataAll,-2,dataAll((i-1)*5 + 1),R, 5);
    RLETF = RLETF/dataAll((i-1)*5 + 1);
 
    if rp > 1.02
        monthDevUp(k) = r / RLETF;
        k= k + 1;
    elseif rp < .98
        monthDevDown(p) = r / RLETF;
        p= p + 1;
    end 
end
 MeanUP = mean(monthDevUp);
 VolUP = std(monthDevUp);
 SkwUP = skewness(monthDevUp);
 KrtUP = kurtosis(monthDevUp);
 PrcUP = prctile(monthDevUp,10);
 StatisticsMUP = [MeanUP VolUP SkwUP KrtUP PrcUP];
 MeanDOWN = mean(monthDevDown);
 VolDOWN = std(monthDevDown);
 SkwDOWN = skewness(monthDevDown);
 KrtDOWN = kurtosis(monthDevDown);
 PrcDOWN = prctile(monthDevDown,10);
 StatisticsMDW = [MeanDOWN VolDOWN SkwDOWN KrtDOWN PrcDOWN];
 
% We clear the stack
 clear MeanUP VolUP SkwUP KrtUP PrcUP
 clear MeanDOWN VolDOWN SkwDOWN KrtDOWN PrcDOWN
 clear i j k N p r R rp weekDev RLETF monthDevDown monthDevUp