clear

Signal = readmatrix('Signal.txt');
n=linspace(0,length(Signal)-1,length(Signal));


% Simple RLS
% Load Data: FF.txt includes -1st colum- Estimated Signal - 2nd colum a0 -3rd colum a1
FF=readmatrix('FF.txt');

figure (1)
subplot(2,2,1);
plot(n,Signal)
hold on;
scatter(n,FF(:,1),2.5,'filled') 
grid on;
title("Generalized RLS_Estimation function, forgeting factor=0.99")
xlabel("number of sample")
ylabel("sample's value")
legend("Random Signal", "Estimated Signal")
hold on;
ylim([0,7])

subplot(2,2,3);
plot(n,FF(:,2))
hold on
plot(n,FF(:,3))
hold on
title("Polinomial Coeficients, forgeting factor=0.99")
grid on;
xlabel("number of sample")
legend("a0", "a1")
ylim([-70,10])


% Rectagular Window Aproch - load data
% Load Data: RecWin.txt includes -1st colum- Estimated Signal - 2nd colum a0 -3rd colum a1
RecWin=readmatrix('RecWin.txt');

subplot(2,2,2);
plot(n,Signal)
hold on;
scatter(n,RecWin(:,1),2.5,'filled') 
grid on;
title("Rectangular Window Approach, Window=50 samples")
xlabel("number of sample")
ylabel("sample's value")
legend("Random Signal", "Estimated Signal")
hold on;
ylim([0,7])

subplot(2,2,4);
plot(n,RecWin(:,2))
hold on
plot(n,RecWin(:,3))
hold on
title("Polinomial Coeficients, Window=50 samples")
grid on;
xlabel("number of sample")
legend("a0", "a1")
ylim([-70,10])
