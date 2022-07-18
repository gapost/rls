clear all;
%Take Data from C++ for graphing in MATLAB%
output =   dlmread("..\rls\MATLAB\TXT-Files\Test_Output.txt");
Test2_a0 = dlmread("..\rls\MATLAB\TXT-Files\Test2_Param_a0.txt");
Test2_a1 = dlmread("..\rls\MATLAB\TXT-Files\Test2_Param_a1.txt");
Test3_a0 = dlmread("..\rls\MATLAB\TXT-Files\Test3_Param_a0.txt");
Test3_a1 = dlmread("..\rls\MATLAB\TXT-Files\Test3_Param_a1.txt");
Test3_a2 = dlmread("..\rls\MATLAB\TXT-Files\Test3_Param_a2.txt");
Test2_est = dlmread("..\rls\MATLAB\TXT-Files\Test2_Est_Output.txt");
Test3_est = dlmread("..\rls\MATLAB\TXT-Files\Test3_Est_Output.txt");
S1 = load('..\rls\MATLAB\MAT_Files\true_a0.mat');
S2 = load('..\rls\MATLAB\MAT_Files\true_a1.mat');
S3 = load('..\rls\MATLAB\MAT_Files\true_a2.mat');

true_a0 = S1.true_a0();
true_a1 = S2.true_a1();
true_a2 = S3.true_a2();

figure(1);

frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

subplot(2,2,1);
numSample = 1:400;
plot(numSample,output(1:400),'r-');
hold on 
plot(numSample,Test2_est(1:400),'b-');
title('Output and Measured Output');
xlabel('Time');
ylabel('Value');
legend('Measured Output','Estimated Output');
grid on 

subplot(2,2,2);
plot(numSample,true_a0(1:400),'r-');
hold on 
plot(numSample,Test2_a0(1:400),'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a0')
legend('Real Parameter','Estimated Parameter');
grid on  

subplot(2,2,3);
plot(numSample,true_a1(1:400),'r-');
hold on 
plot(numSample,Test2_a1(1:400),'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Parameter a1')
legend('Real Parameter','Estimated Parameter');
grid on  

print('C:\Users\nicks\rls\MATLAB\images\RLS_Test2','-dpng','-r0');

figure(2);

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

print('..\rls\MATLAB\images\RLS','-dpng',-'r0');