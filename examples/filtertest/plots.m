clear
##cd build
Signal = dlmread('/home/vd/rls/examples/filtertest/build/Signal.txt');
##Signal = dlmread('/home/vd/rls/examples/filtertest/build/Signal.txt');
n=linspace(0,length(Signal)-1,length(Signal));

% Load Data:
% FF.txt includes -1st colum- Estimated Signal - 2nd colum a0 -3rd colum a1
##FF=dlmread('/home/vd/intership_mat/Intership/rls/examples/filtertest/build/FF.txt');
FF=dlmread('/home/vd/rls/examples/filtertest/build/FF.txt');
% RecWin.txt includes -1st colum- Estimated Signal - 2nd colum a0 -3rd colum a1
RecWin=dlmread('/home/vd/rls/examples/filtertest/build/RecWin.txt');

% Simple RLS

figure (2)
subplot(2,2,1);
plot(n,Signal,'m')
hold on;
scatter(n,FF(:,1),2.5,'filled')
grid on;
title("Generalized RLS_Estimation function, forgeting factor=0.97")
xlabel("number of sample")
ylabel("sample's value")
legend("Random Signal", "Estimated Signal")
hold on;grid on;
ylim([0,max(RecWin(:,1))+1])
xlim([0,length(n)])
##
subplot(2,2,3);
plot(n,FF(:,2))
hold on
plot(n,FF(:,3))
hold on
title("Polinomial Coeficients, forgeting factor=0.97")
grid on;
xlabel("number of sample")
legend("a0", "a1")
ylim([min(RecWin(:,2))-2,max(RecWin(:,2))+2])
xlim([0,length(n)])


% Rectagular Window Aproch - load data

subplot(2,2,2);
plot(n,Signal,'m')
hold on;
scatter(n,RecWin(:,1),2.5,'filled')
hold on;
##plot(n,ones(1,length(n)),'k')
grid on;
##plot(n,5*ones(1,length(n)),'k')
grid on;
title("Rectangular Window Approach, Window=100 samples")
xlabel("number of sample")
ylabel("sample's value")
legend("Random Signal", "Estimated Signal")
ylim([0,max(RecWin(:,1))+1])
hold on;
xlim([0,length(n)])
##
subplot(2,2,4);
plot(n,RecWin(:,2))
hold on
plot(n,RecWin(:,3))
hold on
title("Polinomial Coeficients, Window=100 samples")
grid on;
xlabel("number of sample")
legend("a0", "a1")
ylim([min(RecWin(:,2))-2,max(RecWin(:,2))+2])
xlim([0,length(n)])

