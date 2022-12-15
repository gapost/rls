clear

% options
M = 40; % block window width
N = 2; % # of poly coefficients
ff = 0.98; % forgetting factor
llt = 0;


Nstart = 30;
Nmid   = 100;
Ns     = 200
Nc     = 300;

Y = zeros(Nc,1);

Y(1) = 10;
for i=2:Nstart,
  Y(i) = Y(i-1);
end
for i=Nstart+1:Nmid,
  Y(i) = Y(i-1) + 1;
end
for i=Nmid+1:Ns,
  Y(i) = Y(i-1)*(1-1/20);
end
for i=Ns+1:Nc,
  Y(i) = 50;
end


noise = rand(size(Y))*2-1;
Y = Y + 0.02*max(Y)*noise;

save in.dat -ascii Y

cmd = '../.build/filter/rlsfilter ';
if M>0,
  if llt,
    opts = ['-n' num2str(N) ' -w' num2str(M) ' -llt'];
  else
    opts = ['-n' num2str(N) ' -w' num2str(M)];
  end
else
  opts = ['-n' num2str(N) ' -ff' num2str(ff)];
end

system([cmd opts ' < in.dat > out.dat']);

A = load('out.dat','-ascii');

figure 1
plot(1:Nc,[Y A(:,1)],'.-')
xlabel('time')
legend('signal','estimation')
title('Output of rlsfilter')

figure 2
plot(1:Nc,A(:,2),'.-')
xlabel('time')
title('Estimated rate')

figure 3
plot(1:Nc,A(:,3),'.-')
xlabel('time')
title('Residual')

delete in.dat out.dat




