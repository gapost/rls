clear

% create test signal
N = 10000;
t = (0:N)';
u = round(t/N*10)/10;
i = find(t>N/2);
u(i) = 1 - round(t(i)/N*10)/10;
u = 0.4*u + 0.4;
offset = -15 + (t/N)*2;
slope = 70 - (t/N)*20;
u = u + 0.01*randn(size(t));
dy = 0.1;
y = slope.*u + offset + dy*randn(size(t));


A = [u ones(size(u)) y];
save in.dat -ascii A

% options
M = 100; % block window width
N = 2; % # of poly coefficients
ff = .95; % forgetting factor
sqrt_upd = 1;

cmd = '../.build/filter/rlsfilter ';
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

if M,
  M0 = M;
else
  M0 = 1/(1-ff);
end

A(1:M0,:) = NaN;

figure 1
subplot(2,2,1)
plot(t,u)
title('Input u(t)')
subplot(2,2,2)
plot(t,y,t,A(:,1))
title('Output y(t) = s*u(t) + y_0')
legend('signal','estimation')
xlabel('t')

subplot(2,2,3)
plot(t,[A(:,3) slope])
xlabel('time')
title('Estimated slope s')
subplot(2,2,4)
plot(t,[A(:,4) offset])
xlabel('time')
title('Estimated offset y_0')

delete in.dat out.dat




