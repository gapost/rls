clear
Theta = dlmread('/home/vd/rls/testing/chol_testing/build/Theta.txt');
Signal = dlmread('/home/vd/rls/testing/chol_testing/build/Signal.txt');

L=length(Theta);
n=linspace(1,L,L);

figure(1 )

# non-recurcive
plot(n,Theta(:,1),"m")
hold on
plot(n,Theta(:,2),"r")
hold on

# recurcive
plot(n,Theta(:,3),"b")
hold on
plot(n,Theta(:,4),"g")
hold on

grid on
xlabel("Samples")
ylabel("Coeficients")
xlim([0,L])
ylim([-10,20])
legend("a0 - non recurcive","a1 - non recurcive", "a0 - recurcive","a1 - recurcive")

figure(2 )
plot(Signal)
