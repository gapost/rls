clear
Theta = dlmread('/home/vd/rls/testing/chol_testing/build/Aug_Theta.txt');
Signal = dlmread('/home/vd/rls/testing/chol_testing/build/Aug_Signal.txt');

L=length(Theta);
n=Theta(:,1);
l1=4.8; l2=5.3; l3=3.4;
line1=l1*ones(L,1);  line2=l2*ones(L,1); line3=l3*ones(L,1);x=zeros(L,2);


l=length(Signal)-length(Theta(:,1))+1
length(n)
length(Signal(l:end))

figure(1 )
subplot(2,1,1)
#Coefficients
plot(Signal(l:end),'b'); hold on
plot(Theta(:,1),"k"); hold on
plot(Theta(:,2),"k"); hold on
plot(line1,'--c'); hold on
plot(line2,'--c'); hold on
plot(line3,'--c'); hold on
grid on
xlabel("Samples")
ylabel("Coeficients")
ylim([-0.5,6])
legend("Initial Signal","a0","a1")

subplot(2,1,2)
#Reconstractive Signal
plot(Signal(l:end),'b'); hold on
plot(Theta(:,3),"r"); hold on
plot(line1,'--c'); hold on
plot(line2,'--c'); hold on
plot(line3,'--c'); hold on
grid on
xlabel("Samples")
ylabel("Signal")
ylim([-0.5,6])
legend("Initial Signal","Reconstractive Signal")

