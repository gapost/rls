clear

Signal = dlmread('build/Signal.txt');
n=linspace(0,length(Signal)-1,length(Signal));

% Load Data:
% FF.txt includes -1st colum- Estimated Signal - 2nd colum a0 -3rd colum a1
##FF=dlmread('/home/vd/intership_mat/Intership/rls/examples/filtertest/build/FF.txt');
FF=dlmread('build/FF.txt');
##RecWin.txt includes -1st colum- Estimated Signal - 2nd colum a0 -3rd colum a1
RecWin=dlmread('build/RecWin.txt');
AugChol=dlmread('build/aug_chol.txt');

% Simple RLS

figure (2)
subplot(2,3,1);
plot(n,Signal,'m')
hold on;
scatter(n,FF(:,1),2.5,'filled')
grid on;
title("Generalized RLS_Estimation function, forgeting factor=0.97")
xlabel("number of sample")
ylabel("sample's value")
legend("Random Signal", "Estimated Signal")
hold on;grid on;
ylim([0,10])
xlim([0,length(n)])

% Rectagular Window Aproch - load data

subplot(2,3,2);
plot(n,Signal,'m')
hold on;
plot(RecWin(:,1))
hold on;
##plot(n,ones(1,length(n)),'k')
grid on;
##plot(n,5*ones(1,length(n)),'k')
grid on;
title("Rectangular Window Approach, Window=50 samples")
xlabel("number of sample")
ylabel("sample's value")
legend("Random Signal", "Estimated Signal")
ylim([0,10])
hold on;
xlim([0,length(n)])

% Augmented Cholescy Aproch - load data

subplot(2,3,3);
plot(n,Signal,'m')
hold on;
plot(AugChol(:,1))
hold on;
##plot(n,ones(1,length(n)),'k')
grid on;
##plot(n,5*ones(1,length(n)),'k')
grid on;
title("Augmented Cholescy Aproch, Window=50 samples")
xlabel("number of sample")
ylabel("sample's value")
legend("Random Signal", "Estimated Signal")
ylim([0,10])
hold on;
xlim([0,length(n)])
##
subplot(2,3,4);
plot(n,FF(:,2))
hold on
plot(n,FF(:,3))
hold on
title("Polinomial Coeficients, forgeting factor=0.97")
grid on;
xlabel("number of sample")
legend("a0", "a1")
ylim([min(AugChol(5:end,2))-2,max(AugChol(5:end,2))+2])
xlim([0,length(n)])


##
subplot(2,3,5);
plot(RecWin(:,2))
hold on
plot(RecWin(:,3))
hold on
title("Polinomial Coeficients, Window=50 samples")
grid on;
xlabel("number of sample")
legend("a0", "a1")
ylim([min(AugChol(5:end,2))-2,max(AugChol(5:end,2))+2])
xlim([0,length(n)])

##
subplot(2,3,6);
plot(AugChol(:,2))
hold on
plot(AugChol(:,3))
hold on
title("Polinomial Coeficients, Window=50 samples")
grid on;
xlabel("number of sample")
legend("a0", "a1")
ylim([min(AugChol(5:end,2))-2,max(AugChol(5:end,2))+2])
xlim([0,length(n)])
