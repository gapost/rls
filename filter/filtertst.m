clear

M = 20; % block window width
N = 2; % # of poly coefficients
ff = 1; % forgetting factor

Nstart = 30;
Nmid   = 100;
Nc     = 200;

Y = zeros(Nc,1);

Y(1) = 10;
for i=2:Nstart,
  Y(i) = Y(i-1);
end
for i=Nstart+1:Nmid,
  Y(i) = Y(i-1) + 1;
end
for i=Nmid+1:Nc,
  Y(i) = Y(i-1)*(1-1/20);
end

save in.dat -ascii Y

cmd = '.build/rlsfilter ';
opts = ['-n' num2str(N) ' -w' num2str(M) ' -ff' num2str(ff)];

system([cmd opts ' < in.dat > out.dat']);

A = load('out.dat','-ascii');

figure 1
plot(1:Nc,[Y A(:,1)],'.-')
xlabel('time')
legend('signal','estimation')
title('Output of rlsfilter')

delete in.dat out.dat




