l=2500; m=1500; n=800;
l1=4.8; l2=5.3; l3=3.4;

##data_proccess( coeff1, coeff2, in_signal,rec_signal,start, l1, l2, l3,n,m,l)

# Cholesky Method
RecW_Theta = dlmread('/home/vd/rls/testing/chol_testing/build/RecWinTheta.txt');

# Non Recurcive
data_proccess( RecW_Theta(:,1), RecW_Theta(:,2),RecW_Theta(:,7) ,RecW_Theta(:,3),25, l1, l2, l3,n,m,l)

##Recurcive
data_proccess( RecW_Theta(:,4), RecW_Theta(:,5),RecW_Theta(:,7) ,RecW_Theta(:,6),25, l1, l2, l3,n,m,l)

## Augmented Cholesky Method
Aug_Theta = dlmread('/home/vd/rls/testing/chol_testing/build/Aug_Theta.txt');
data_proccess( Aug_Theta(:,1), Aug_Theta(:,2),Aug_Theta(:,4) ,Aug_Theta(:,3),25, l1, l2, l3,n,m,l)

