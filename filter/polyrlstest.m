clear

% options
M = 100; % block window width
N = 2; % # of poly coefficients
ff = .95; % forgetting factor
sqrt_upd = 1;


[Y, th] = test_signal(1000,10);
noise = randn(size(Y));
dY = 0.1;
Y = Y + dY*noise;

save in.dat -ascii Y

cmd = '../.build/filter/polyrlsfilter ';
if M>0,
  if sqrt_upd,
    opts = ['-n' num2str(N) ' -w' num2str(M) ' -sqrt'];
  else
    opts = ['-n' num2str(N) ' -w' num2str(M)];
  end
else
  opts = ['-n' num2str(N) ' -ff' num2str(ff)];
end

tic
system([cmd opts ' < in.dat > out.dat']);
toc

A = load('out.dat','-ascii');

t = (1:length(Y))';

if M,
  M0 = M;
else
  M0 = 1/(1-ff);
end
A(1:M0,:) = NaN;


figure 1
clf
subplot(2,2,1)
plot(t,[Y A(:,1)],'.-')
xlabel('time')
legend('signal','estimation')
title('Timeseries y(t)')

subplot(2,2,2)
plot(t,[A(:,2) th(:,2)],'.-')
xlabel('time')
title('Estimated rate dy/dt')

subplot(2,2,3)
plotyy(t,A(:,4),t,A(:,5))
xlabel('time')
title('Estimated Parameters')
legend('\theta_0','\theta_1')

subplot(2,2,4)
semilogy(t,A(:,3),'.-',[t(1) t(end)],[1 1]*dY^2*M0)
xlabel('time')
title('Residual')


delete in.dat out.dat




