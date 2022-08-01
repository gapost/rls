clear all;
obj = recursiveLS(3);
obj.ForgettingFactor = 0.9;

output =   dlmread("..\rls\testing\exp_testing\TXT-Files\PolyRLS\Test_Output.txt");
S2 = load('..\rls\testing\exp_testing\MAT_Files\input.mat');
S3 = load('..\rls\testing\exp_testing\MAT_Files\true_a0.mat');
S4 = load('..\rls\testing\exp_testing\MAT_Files\true_a1.mat');
S5 = load('..\rls\testing\exp_testing\MAT_Files\true_a2.mat');

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

frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

subplot(2,2,1);
numSample = 1:numel(input);
plot(numSample,output,'r-');
hold on 
plot(numSample,estimatedOut,'b-');
title('Output and Measured Output');
xlabel('Time');
ylabel('Value');
legend('Measured Output','Estimated Output');
grid on 

subplot(2,2,2);
plot(numSample,true_a0,'r-');
hold on 
plot(numSample,theta_full(:,3),'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a0')
legend('Real Parameter','Estimated Parameter');
grid on  

subplot(2,2,3);
plot(numSample,true_a1,'r-');
hold on 
plot(numSample,theta_full(:,2),'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a1')
legend('Real Parameter','Estimated Parameter');
grid on 

subplot(2,2,4);
plot(numSample,true_a2,'r-');
hold on 
plot(numSample,theta_full(:,1),'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a2');
legend('Real Parameter','Estimated Parameter');
grid on 

    
print('..\rls\testing\exp_testing\images\Test3Check','-dpng','-r0');