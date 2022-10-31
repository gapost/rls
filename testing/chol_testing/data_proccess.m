function data_proccess( coeff1, coeff2, in_signal,rec_signal,start, l1, l2, l3,n,m,l)

  line1=l1*ones(l,1);  line2=l2*ones(l,1); line3=l3*ones(l,1);
  clear_signal=zeros(1,length(rec_signal));
  clear_signal(1:n-start)=transpose(line1(1:n-start,1));
  clear_signal(n-start:m-start-1)=transpose(line2(1:m-n,1));
length(clear_signal)
-m+start+1
  clear_signal(m-start:end)=transpose(line3(1:length(clear_signal)-m+start+1,1));
  #Display Data
  figure()
  subplot(3,1,1)
  #Coefficients
  plot(coeff1,"k"); hold on
  plot(coeff2,"k"); hold on
  plot(line1,'--c'); hold on
  plot(line2,'--c'); hold on
  plot(line3,'--c'); hold on
  grid on
  xlabel("Samples")
  ylabel("Coeficients")
  ylim([-0.5,max([l1 l2 l3])+2]); xlim([0,l])
  legend("a0","a1")

  subplot(3,1,2)
  plot(in_signal);hold on
  plot(rec_signal,"r"); hold on
  ##plot(line1,'--c'); hold on
  ##plot(line2,'--c'); hold on
  ##plot(line3,'--c'); hold on
  grid on
  xlabel("Samples")
  ylim([-0.5,max([l1 l2 l3])+2]); xlim([0,l])
  legend("Initial Signal","Reconstractive Signal")

  subplot(3,1,3)
  plot(clear_signal);hold on
  plot(rec_signal,"r"); hold on
  ##plot(line1,'--c'); hold on
  ##plot(line2,'--c'); hold on
  ##plot(line3,'--c'); hold on
  grid on
  xlabel("Samples")
  ylim([-0.5,max([l1 l2 l3])+2]); xlim([0,l])
  legend("Ideal Signal","Reconstractive Signal")
  hold off
  # Calculate snr
  signal_tr=transpose(rec_signal);
##  length(signal_tr)
##  length(clear_signal)
  noise=transpose(rec_signal)-clear_signal;
  snr= 20*log(sum(clear_signal.^2)/sum(noise.^2))
end
