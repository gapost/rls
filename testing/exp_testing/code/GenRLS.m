clear all;
%Take Data from C++ for graphing in testing\exp_testing%
output =   dlmread("..\rls\testing\exp_testing\TXT-Files\PolyRLS\Test_Output.txt");
Test3_a0 = dlmread("..\rls\testing\exp_testing\TXT-Files\GenRLS\TestGen_Param_a0.txt");
Test3_a1 = dlmread("..\rls\testing\exp_testing\TXT-Files\GenRLS\TestGen_Param_a1.txt");
Test3_a2 = dlmread("..\rls\testing\exp_testing\TXT-Files\GenRLS\TestGen_Param_a2.txt");
Test3_est = dlmread("..\rls\testing\exp_testing\TXT-Files\GenRLS\TestGen_Est_Output.txt");
S1 = load('..\rls\testing\exp_testing\MAT_Files\true_a0.mat');
S2 = load('..\rls\testing\exp_testing\MAT_Files\true_a1.mat');
S3 = load('..\rls\testing\exp_testing\MAT_Files\true_a2.mat');

true_a0 = S1.true_a0();
true_a1 = S2.true_a1();
true_a2 = S3.true_a2();


figure(1);

frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);


subplot(2,2,1);
numSample = 1:500;
plot(numSample,output,'r-');
hold on 
plot(numSample,Test3_est,'b-');
title('Output and Measured Output');
xlabel('Time');
ylabel('Value');
legend('Measured Output','Estimated Output');
grid on 

subplot(2,2,2);
numSample = 1:500;
plot(numSample,true_a0,'r-');
hold on 
plot(numSample,Test3_a0,'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a0')
legend('Real Parameter','Estimated Parameter');
grid on  

subplot(2,2,3);
plot(numSample,true_a1,'r-');
hold on 
plot(numSample,Test3_a1,'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a1')
legend('Real Parameter','Estimated Parameter');
grid on  


subplot(2,2,4); 
plot(numSample,true_a2,'r-');
hold on 
plot(numSample,Test3_a2,'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a2')
legend('Real Parameter','Estimated Parameter');
grid on  

print('..\rls\testing\exp_testing\images\GenRLS','-dpng','-r0');