clear
Theta = dlmread('/home/vd/rls/testing/chol_testing/build/Theta.txt');
Signal = dlmread('/home/vd/rls/testing/chol_testing/build/Signal.txt');

L=length(Theta);
n=linspace(1,L,L);

figure(1 )

# non-recurcive
plot(n,Theta(:,1),"m")
hold on
plot(n,Theta(:,2),"m")
hold on

# recurcive
plot(n,Theta(:,3),"g")
hold on
plot(n,Theta(:,4),"g")
hold on

grid on
xlabel("Samples")
ylabel("Coeficients")
xlim([0,L])
ylim([-2,8])
legend("a0 - non recurcive","a1 - non recurcive", "a0 - recurcive","a1 - recurcive","Signal")
hold on
##figure(2 )
plot(Signal,'y')
##ylim([0,6])
