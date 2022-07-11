clear all
%Take Data from C++ for graphing in MATLAB%
Test2_a0 = dlmread("C:\Users\nicks\Downloads\Test2_Param_a0.txt");
Test2_a1 = dlmread("C:\Users\nicks\Downloads\Test2_Param_a1.txt");
Test3_a0 = dlmread("C:\Users\nicks\Downloads\Test3_Param_a0.txt");
Test3_a1 = dlmread("C:\Users\nicks\Downloads\Test3_Param_a1.txt");
Test3_a2 = dlmread("C:\Users\nicks\Downloads\Test3_Param_a2.txt");

figure(1)
numSample = 1:200;
plot(numSample,Test2_a0,'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Test2 - Parameter a0')
legend('Estimated Parameters')
grid on


figure(2)
plot(numSample,Test2_a1,'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Test2 - Parameter a1')
legend('Estimated Parameters')
grid on

figure(3)
numSample = 1:250;
plot(numSample,Test3_a0,'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Test3 - Parameter a0')
legend('Estimated Parameters')
grid on

figure(4)
plot(numSample,Test3_a1,'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Test3 - Parameter a1')
legend('Estimated Parameters')
grid on


figure(5) 
plot(numSample,Test3_a2,'b-')
xlabel('Time');
ylabel('Value of Parameter');
title('Test3 - Parameter a2')
legend('Estimated Parameters')
grid on
