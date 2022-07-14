clear all;
%Take Data from C++ for graphing in MATLAB%
output =   dlmread("..\rls\MATLAB\TXT-Files\Test_Output.txt");
Test2_a0 = dlmread("..\rls\MATLAB\TXT-Files\Test2_Param_a0.txt");
Test2_a1 = dlmread("..\rls\MATLAB\TXT-Files\Test2_Param_a1.txt");
Test3_a0 = dlmread("..\rls\MATLAB\TXT-Files\Test3_Param_a0.txt");
Test3_a1 = dlmread("..\rls\MATLAB\TXT-Files\Test3_Param_a1.txt");
Test3_a2 = dlmread("..\rls\MATLAB\TXT-Files\Test3_Param_a2.txt");


figure(1);
subplot(2,3,1);
numSample = 1:250;
plot(numSample,output,'b');
xlabel('Time');
ylabel('Output');
title('Output')
grid on

subplot(2,3,2);
plot(numSample,Test2_a0,'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Test2 - Parameter a0')
grid on


subplot(2,3,3);
plot(numSample,Test2_a1,'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Test2 - Parameter a1')
grid on

subplot(2,3,4);
numSample = 1:250;
plot(numSample,Test3_a0,'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Test3 - Parameter a0')
grid on

subplot(2,3,5);
plot(numSample,Test3_a1,'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Test3 - Parameter a1')
grid on


subplot(2,3,6); 
plot(numSample,Test3_a2,'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Test3 - Parameter a2')
grid on
print('..\rls\MATLAB\images\RLS','-dpng');