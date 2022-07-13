clear all;
obj = recursiveLS(2);
obj.ForgettingFactor = 0.89;

S1 = load('C:\Users\nicks\rls\MATLAB\MAT_Files\output.mat');
S2 = load('C:\Users\nicks\rls\MATLAB\MAT_Files\input.mat');
S3 = load('C:\Users\nicks\rls\MATLAB\MAT_Files\error1.mat');
S4 = load('C:\Users\nicks\rls\MATLAB\MAT_Files\error2.mat');
S5 = load('C:\Users\nicks\rls\MATLAB\MAT_Files\error3.mat');

output = S1.output();
noise = wgn(1,250,25);
output = noise + output;

input  = S2.input();
error1 = S3.error1();
error2 = S4.error2();
error3 = S5.error3();
    
for i = 1:numel(input)
    H = [input(i) 1];
    [theta, EstimatedOutput] = obj(output(i),H);
    estimatedOut(i)= EstimatedOutput;
    theta_est(i,:) = max(theta);
    theta_full(i,:) = theta;
    error1(i) = error1(i) - theta(1);
    error2(i) = error2(i) - theta(2);
    oldInput = input(i);
end

figure(1);
subplot(2,3,1);
numSample = 1:numel(input);
plot(numSample,output,'b');
hold on 
plot(numSample,estimatedOut,'r-');
title('Output and Measured Output');
xlabel('Time');
ylabel('Value');
legend('Measured Output','Estimated Output');
grid on 

subplot(2,3,2);
plot(numSample,theta_full(:,1),'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a0')
grid on 


subplot(2,3,3);
plot(numSample,theta_full(:,2),'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a1')
grid on 

subplot(2,3,4);
plot(numSample,error1,'b-')
xlabel('Time');
ylabel('Error');
title('Error of Parameter a0')
grid on 

subplot(2,3,5);
plot(numSample,error2,'b-')
xlabel('Time');
ylabel('Error');
title('Error of Parameter a1')
grid on 

print('C:\Users\nicks\rls\MATLAB\images\Test2Check','-dpng');