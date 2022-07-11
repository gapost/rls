clear all
obj = recursiveLS(3);
obj.ForgettingFactor = 0.5;
output = [];
input = [];

for i = 1:50
    output(i) = i;
    input(i) = i;
end 

for i = 50:100
    output(i) = 50;
    input(i) = i;
end 

for i = 100:150
    output(i) = 50 + i;
    input(i) = i;
end 

for i = 150:200
    output(i) = 20;
    input(i) = i;
end 

for i = 200:250
    output(i) = (i-198)*(i-198) + 2*(i-198) + 1;
    input(i) = i;
 
end    
for i = 1:numel(input)
    H = [input(i)*input(i) input(i) 1];
    [theta, EstimatedOutput] = obj(output(i),H);
    estimatedOut(i)= EstimatedOutput;
    theta_est(i,:) = max(theta);
    theta_full(i,:) = theta;
    oldInput = input(i);
end

figure(1) 
numSample = 1:numel(input);
plot(numSample,output,'b');
hold on
plot(numSample,estimatedOut,'r-');
legend('Measured Output','Estimated Output');
grid on 

figure(2)
plot(numSample,theta_full(:,1),'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a2')
legend('Estimated Parameters')
grid on 


figure(3)
plot(numSample,theta_full(:,2),'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a1')
legend('Estimated Parameters')
grid on 

figure(4)
plot(numSample,theta_full(:,3),'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a0')
legend('Estimated Parameters')
grid on 