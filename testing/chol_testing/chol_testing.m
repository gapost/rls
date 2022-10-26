clear
#Step 1-> check Cholescy Decomposition
l=2500; m=1500; n=800;
l1=4.8; l2=6.3; l3=4;

Phi=ones(l,2);
Y=zeros(l,1);

R=randi(10,n,1); R=(R-mean(R))/max(R); Y(1:n,1)=l1+R;

R=randi(10,m-n,1); R=(R-mean(R))/max(R); Y(n+1:m,1)=l2+R;

R=randi(10,l-m,1); R=(R-mean(R))/max(R); Y(m+1:l,1)=l3+R;

line1=l1*ones(l,1); line2=l2*ones(l,1); line3=l3*ones(l,1);
x=zeros(l,2);
for i=2:l
  Phi(i,2)=i;
  x(i,:)=lschol(Phi(1:i,:),Y(1:i,:));
end

figure(1)
subplot(2,2,1)
plot(Y,"k"); hold on
plot(x(:,1),"m"); hold on
plot(x(:,2),"b"); hold on
plot(line1,'--k'); hold on
plot(line2,'--k'); hold on
plot(line3,'--k'); hold on
ylim([-0.5,7]); xlim([-l/100,l])
title("Data Calculated with Octave")
legend("Initial Signal","Coeff a0","Coeff a1")
grid on

#Step 2 -> Signal Reconstraction

S=zeros(l,1);
for i=1:l
  S(i,1)=Phi(i,:)*x(i,:)';

end

figure(1)
subplot(2,2,3)
plot(Y,"k"); hold on
plot(S,"m"); hold on
plot(line1,'--k'); hold on
plot(line2,'--k'); hold on
plot(line3,'--k'); hold on
legend("Initial Signal","Reconstracted Signal")
ylim([-0.5,7]); xlim([-l/100,l])
grid on

#Step 3 -> Display Data calculate through C++ code

Theta = dlmread('/home/vd/rls/testing/chol_testing/build/Theta.txt');
Signal = dlmread('/home/vd/rls/testing/chol_testing/build/Signal.txt');

l=length(Theta);
st=length(Signal)-length(Theta(:,1))+1

figure(1)
subplot(2,2,2)
# non-recurcive
plot(Theta(:,1),"m"); hold on
plot(Theta(:,2),"m"); hold on

# recurcive
plot(Theta(:,4),"k"); hold on
plot(Theta(:,5),"k"); hold on
plot(Signal(st:end),'b'); hold on
plot(line1,'--k'); hold on
plot(line2,'--k'); hold on
plot(line3,'--k'); hold on
title("Data Calculated with C++")

legend("Initial Signal","Reconstracted Signal")
ylim([-0.5,7]); xlim([-l/100,l])
grid on
xlabel("Time"); ylabel("Coeficients")
##xlim([0,(n)])
ylim([-0.5,7]); xlim([-l/100,l])
legend("a0 - non recurcive","a1 - non recurcive", "a0 - recurcive","a1 - recurcive","Signal")
##ylim([0,6])

subplot(2,2,4)
plot(Signal(st:end),'b'); hold on
plot(Theta(:,3),"m"); hold on
plot(Theta(:,6),"g"); hold on
plot(line1,'--k'); hold on
plot(line2,'--k'); hold on
plot(line3,'--k'); hold on
legend("Initial Signal","Reconstracted Signal without recurtion","Reconstracted Signal with recurtion")
ylim([-0.5,7]); xlim([-l/100,l])
grid on





