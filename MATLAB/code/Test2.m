clear all;
obj = recursiveLS(2);
obj.ForgettingFactor = 0.9;

S1 = load('..\rls\MATLAB\MAT_Files\output.mat');
S2 = load('..\rls\MATLAB\MAT_Files\input.mat');
S3 = load('..\rls\MATLAB\MAT_Files\error1.mat');
S4 = load('..\rls\MATLAB\MAT_Files\error2.mat');
S5 = load('..\rls\MATLAB\MAT_Files\error3.mat');

output = S1.output();
output = output(1:200);

noise = wgn(1,200,1);
output = noise + output;

input  = S2.input();
true_a0 = S3.true_a0();
true_a1 = S4.true_a1();
%%true_a2 = S5.true_a2();

input = input(1:200);
true_a0 = true_a0(1:200);
true_a1 = true_a1(1:200);

for i = 1:numel(input)
    H = [input(i) 1];
    [theta, EstimatedOutput] = obj(output(i),H);
    estimatedOut(i)= EstimatedOutput;
    theta_est(i,:) = max(theta);
    theta_full(i,:) = theta;
end

figure(1);
subplot(2,2,1);
numSample = 1:200;
plot(numSample,output,'b');
hold on 
plot(numSample,estimatedOut,'r-');
title('Output and Measured Output');
xlabel('Time');
ylabel('Value');
legend('Estimated Output','Measured Output');
grid on 

subplot(2,2,2);
plot(numSample,true_a1,'r-');
hold on 
plot(numSample,theta_full(:,1),'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a1')
legend('Estimated Parameter','Measured Parameter');
grid on 


subplot(2,2,3);
plot(numSample,true_a0,'r-');
hold on 
plot(numSample,theta_full(:,2),'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a0')
legend('Estimated Parameter','Measured Parameter');
grid on  

print('..\rls\MATLAB\images\Test2Check','-dpng');