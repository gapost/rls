clear all;
obj = recursiveLS(3);
obj.ForgettingFactor = 0.9;

S1 = load('..\rls\MATLAB\MAT_Files\output.mat');
S2 = load('..\rls\MATLAB\MAT_Files\input.mat');
S3 = load('..\rls\MATLAB\MAT_Files\error1.mat');
S4 = load('..\rls\MATLAB\MAT_Files\error2.mat');
S5 = load('..\rls\MATLAB\MAT_Files\error3.mat');

output = S1.output();
noise = wgn(1,250,1);
output = noise + output;

input  = S2.input();
true_a0 = S3.true_a0();
true_a1 = S4.true_a1();
true_a2 = S5.true_a2();


for i = 1:numel(input)
    H = [input(i)*input(i) input(i) 1];
    [theta, EstimatedOutput] = obj(output(i),H);
    estimatedOut(i)= EstimatedOutput;
    theta_full(i,:) = theta;
end

figure(1);
subplot(2,2,1);
numSample = 1:numel(input);
plot(numSample,output,'b');
hold on 
plot(numSample,estimatedOut,'r-');
title('Output and Measured Output');
xlabel('Time');
ylabel('Value');
legend('Measured Output','Estimated Output');
grid on 

subplot(2,2,2);
plot(numSample,true_a2,'b-');
hold on 
plot(numSample,theta_full(:,1),'r-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a2');
legend('Estimated Parameter','Measured Parameter');
grid on 


subplot(2,2,3);
plot(numSample,true_a1,'b-');
hold on 
plot(numSample,theta_full(:,2),'r-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a1')
legend('Estimated Parameter','Measured Parameter');
grid on 

subplot(2,2,4);
plot(numSample,true_a0,'b-');
hold on 
plot(numSample,theta_full(:,3),'r-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a0')
legend('Estimated Parameter','Measured Parameter');
grid on  

print('..\rls\MATLAB\images\Test3Check','-dpng');