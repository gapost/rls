#include <cstdlib>
#include <random>
#include <fstream>
#include "C:\armadillo-11.2.1\examples\RLS_Estimation_Object.h"

using namespace RLS;
using namespace arma;
using namespace std;

int main() {
	//Initializing parameters//
	const int len = 500; //Length of Array of Output/Input
	double init_covar = 100000.;
	//Generate Data for Validation
	double Y[len] = { 0. };
	double X[len] = { 0. };
	

	fstream OutFile("C:/Users/nicks/rls/MATLAB/TXT-Files/PolyRLS/Test_Output.txt", ios::in);
	fstream InFile("C:/Users/nicks/rls/MATLAB/TXT-Files/PolyRLS/Test_Input.txt", ios::in);

	double value;
	int i = 0;
	while (OutFile >> value)
	{
		Y[i] = value;
		i++;
	}

	i = 0;

	while (InFile >> value)
	{
		X[i] = value;
		i++;
	}


	

	//Test the "Set" and "Get" Functions in the class//
	PolyRLS<double, 2> Test_Alg(1., 10000);
	//Do some updates to notice changes//
	for (int i = 0; i < 10; i++) {
		Test_Alg.update_par(Y[i]); // Update parameters 
	}
	cout << "Testing Functions : " <<  endl;
	cout << "Parameter function : " << endl;
	Test_Alg.getEstimatedParameters().print();
	cout << '\n';
	cout << "Covariance function : " << endl;
	Test_Alg.getCovarianceMat().print();
	cout << '\n';
	cout << "Gain function : " << endl;
	Test_Alg.getGains().print();
	cout << '\n';
	cout << "Iteration function : " << Test_Alg.getIterations() << endl;
	cout << '\n';
	cout << "Forgetting factor function : " << Test_Alg.getLambda() << endl;
	cout << '\n';
	cout << "Initial Covariance function : " << Test_Alg.getCovar() << endl;
	cout << '\n';

	//Change some parameters with "Set" and test them again//
	Test_Alg.setLambda(0.5);
	Test_Alg.setCovariance(5000);
	cout << "New forgetting factor : " << Test_Alg.getLambda() << endl;
	cout << '\n';
	cout << "New initial covariancce : " << Test_Alg.getCovar() << endl;
	cout << '\n';

	//Test "reset" function//
	cout << "Testing reset function : " << endl; 
	Test_Alg.reset();
	cout << "Parameter function : " << endl;
	Test_Alg.getEstimatedParameters().print();
	cout << '\n';
	cout << "Covariance function : " << endl;
	Test_Alg.getCovarianceMat().print();
	cout << '\n';
	cout << "Gain function : " << endl;
	Test_Alg.getGains().print();
	cout << '\n';
	cout << "Iteration function : " << Test_Alg.getIterations() << endl;
	cout << '\n';
	cout << "Forgetting factor function : " << Test_Alg.getLambda() << endl;
	cout << '\n';
	cout << "Initial Covariance function : " << Test_Alg.getCovar() << endl;
	cout << '\n';

	//Initialize some different parameters to check differences//
	//Test 1 : 2 Parameters and forgetting factor 1

	PolyRLS<double, 2> Alg_1(1., init_covar);

	cout << "|---------First Test : We expect parameters for each area to have errors as every single input is taken in regard---------|" << endl;
	cout << "|---------and doesn't accurately represent the area.---------|" << endl;

	for (int i = 0; i < 500; i++) {
		cout << "Here is output: " << Y[i] << endl;
		Alg_1.update_par(Y[i]); // Update parameters in respect to Input and Output
		cout << "Here are the estimated parameters at time " << i <<" : " << endl;
		Alg_1.getEstimatedParameters().print(); //Print Matrix 
		cout << '\n';
		
	}
	
	//Test 2 : 2 Parameters and forgetting factor 0.9

	PolyRLS<double, 2> Alg_2(0.9, init_covar);

	cout << "|---------Second Test : We expect parameters to be time-variant and change depending on the output---------|" << endl;
	cout << "|---------In our case , we should see that it resembles the polynomial we put as input.---------|" << endl;

	//Writing Data from Parameters//
	ofstream outFile1;
	ofstream outFile2;
	ofstream outFile_out_2;
	outFile1.open("C:/Users/nicks/rls/MATLAB/TXT-Files/PolyRLS/Test2_Param_a0.txt");
	outFile2.open("C:/Users/nicks/rls/MATLAB/TXT-Files/PolyRLS/Test2_Param_a1.txt");
	outFile_out_2.open("C:/Users/nicks/rls/MATLAB/TXT-Files/PolyRLS/Test2_Est_Output.txt");

	for (int i = 0; i < 500; i++) {
		cout << "Here is estimated output: " << Alg_2.getEstimatedOutput() << endl;
		outFile_out_2 << Alg_2.getEstimatedOutput() << " ";
		cout << "Here is output: " << Y[i] << endl;
		Alg_2.update_par(Y[i]); // Update parameters in respect to Input and Output
		cout << "Here are the estimated parameters at time " << i << " : " << endl;
		Alg_2.getEstimatedParameters().print(); //Print Matrix 
		cout << '\n';
		outFile1 << Alg_2.getEstimatedParameters().row(0) << " ";
		outFile2 << Alg_2.getEstimatedParameters().row(1) << " ";
	}
	outFile1.close();
	outFile2.close();
	outFile_out_2.close();

	//Test 3 : 3 Parameters and forgetting factor 0.9

	PolyRLS<double, 3> Alg_3(0.9, init_covar);

	cout << "|---------Third Test : We are estimating everything as a 2nd order Polynomial---------|"<< endl;

	//Writing Data from Parameters//
	ofstream outFile3;
	ofstream outFile4;
	ofstream outFile5;
	ofstream outFile_out_3;
	outFile3.open("C:/Users/nicks/rls/MATLAB/TXT-Files/PolyRLS/Test3_Param_a0.txt");
	outFile4.open("C:/Users/nicks/rls/MATLAB/TXT-Files/PolyRLS/Test3_Param_a1.txt");
	outFile5.open("C:/Users/nicks/rls/MATLAB/TXT-Files/PolyRLS/Test3_Param_a2.txt");
	outFile_out_3.open("C:/Users/nicks/rls/MATLAB/TXT-Files/PolyRLS/Test3_Est_Output.txt");

	for (int i = 0; i < 500; i++) {
		cout << "Here is estimated output: " << Alg_3.getEstimatedOutput() << endl;
		outFile_out_3 << Alg_3.getEstimatedOutput() << " ";
		cout << "Here is output: " << Y[i] << endl;
		Alg_3.update_par(Y[i]); // Update parameters in respect to Input and Output
		cout << "Here are the estimated parameters at time " << i << " : " << endl;
		Alg_3.getEstimatedParameters().print(); //Print Matrix 
		cout << '\n';
		cout << "Here is cost : " << Alg_3.getCost() << endl;
		outFile3 << Alg_3.getEstimatedParameters().row(0) << " ";
		outFile4 << Alg_3.getEstimatedParameters().row(1) << " ";
		outFile5 << Alg_3.getEstimatedParameters().row(2) << " ";
	}
	outFile3.close();
	outFile4.close();
	outFile5.close();
	outFile_out_3.close();

	cout << "|---------Takes much longer to converge with this method ---------|" << endl;
	cout << "|---------Need to look for big shifts in parameters,to validate the change,ie 50 ->6500 at that moment---------|" << endl;;
	cout << "|--------- or 0.0002 -> -5.8---------|";
	cout << '\n';

	//Test 4 : Testing Generalized RLS_Estimation function
	cout << "|---------Now we are going to test the generalized RLS function with the regressors of the---------|" << endl;
	cout << "|---------polynomial to check the results---------|" << endl;
	cout << '\n';
	vec Reg(3, fill::eye);
	RLS_Estimator<double, 3> Test_Gen(0.9,100000.);
	//Do some updates to notice changes//
	for (int i = 0; i < 10; i++) {
		Reg(1) = i;
		Reg(2) = i * i;
		Test_Gen.update_par(Reg,Y[i]); // Update parameters 
	}

	cout << "Testing Functions : " << endl;
	cout << "Parameter function : " << endl;
	Test_Gen.getEstimatedParameters().print();
	cout << '\n';
	cout << "Covariance function : " << endl;
	Test_Gen.getCovarianceMat().print();
	cout << '\n';
	cout << "Gain function : " << endl;
	Test_Gen.getGains().print();
	cout << '\n';
	cout << "Iteration function : " << Test_Gen.getIterations() << endl;
	cout << '\n';
	cout << "Forgetting factor function : " << Test_Gen.getLambda() << endl;
	cout << '\n';
	cout << "Initial Covariance function : " << Test_Gen.getCovar() << endl;
	cout << '\n';

	//Change some parameters with "Set" and test them again//
	Test_Gen.setLambda(0.5);
	Test_Gen.setCovariance(5000);
	cout << "New forgetting factor : " << Test_Gen.getLambda() << endl;
	cout << '\n';
	cout << "New initial covariancce : " << Test_Gen.getCovar() << endl;
	cout << '\n';

	//Test "reset" function//
	cout << "Testing reset function : " << endl;
	Test_Gen.reset();
	cout << "Parameter function : " << endl;
	Test_Gen.getEstimatedParameters().print();
	cout << '\n';
	cout << "Covariance function : " << endl;
	Test_Gen.getCovarianceMat().print();
	cout << '\n';
	cout << "Gain function : " << endl;
	Test_Gen.getGains().print();
	cout << '\n';
	cout << "Iteration function : " << Test_Gen.getIterations() << endl;
	cout << '\n';
	cout << "Forgetting factor function : " << Test_Gen.getLambda() << endl;
	cout << '\n';
	cout << "Initial Covariance function : " << Test_Gen.getCovar() << endl;
	cout << '\n';

	//Reset back to start
	Test_Gen.reset();
	Test_Gen.setLambda(1.);
	Test_Gen.setCovariance(10000);

	//Test 5: Check output of parameters from the generalized class
	ofstream outFileGen1;
	ofstream outFileGen2;
	ofstream outFileGen3;
	ofstream outFile_out_Gen;
	outFileGen1.open("C:/Users/nicks/rls/MATLAB/TXT-Files/GenRLS/TestGen_Param_a0.txt");
	outFileGen2.open("C:/Users/nicks/rls/MATLAB/TXT-Files/GenRLS/TestGen_Param_a1.txt");
	outFileGen3.open("C:/Users/nicks/rls/MATLAB/TXT-Files/GenRLS/TestGen_Param_a2.txt");
	outFile_out_Gen.open("C:/Users/nicks/rls/MATLAB/TXT-Files/GenRLS/TestGen_Est_Output.txt");

	for (int i = 0; i < 500; i++) {
		Reg(1) = i;
		Reg(2) = i * i;
		cout << "Here is estimated output: " << Test_Gen.getEstimatedOutput(Reg) << endl;
		outFile_out_Gen << Test_Gen.getEstimatedOutput(Reg) << " ";
		cout << "Here is output: " << Y[i] << endl;
		Test_Gen.update_par(Reg,Y[i]); // Update parameters in respect to Input and Output
		cout << "Here are the estimated parameters at time " << i << " : " << endl;
		Test_Gen.getEstimatedParameters().print(); //Print Matrix 
		cout << '\n';
		cout << "Here is cost : " << Test_Gen.getCost() << endl;
		outFileGen1 << Test_Gen.getEstimatedParameters().row(0) << " ";
		outFileGen2 << Test_Gen.getEstimatedParameters().row(1) << " ";
		outFileGen3 << Test_Gen.getEstimatedParameters().row(2) << " ";
	}

	outFileGen1.close();
	outFileGen2.close();
	outFileGen3.close();
	outFile_out_Gen.close();

	//Test 6: Testing Rectangular Window Approach
	BlockRLS<double, 3> Test_Block(1,10, 10000.);
	cout << "|---------Now we are going to test the Block RLS function with the regressors of the---------|" << endl;
	cout << "|---------polynomial to check the results---------|" << endl;
	cout << '\n';
	//Test 6: Check output of parameters from the BlockRLS class
	ofstream outFileBlock1;
	ofstream outFileBlock2;
	ofstream outFileBlock3;
	ofstream outFile_out_Block;
	outFileBlock1.open("C:/Users/nicks/rls/MATLAB/TXT-Files/BlockRLS/TestBlock_Param_a0.txt");
	outFileBlock2.open("C:/Users/nicks/rls/MATLAB/TXT-Files/BlockRLS/TestBlock_Param_a1.txt");
	outFileBlock3.open("C:/Users/nicks/rls/MATLAB/TXT-Files/BlockRLS/TestBlock_Param_a2.txt");
	outFile_out_Block.open("C:/Users/nicks/rls/MATLAB/TXT-Files/BlockRLS/TestBlock_Est_Output.txt");
	for (int i = 0; i < 500; i++) {
		Reg(1) = i;
		Reg(2) = i * i;
		cout << "Here is estimated output: " << Test_Block.getEstimatedOutput(Reg) << endl;
		outFile_out_Block << Test_Block.getEstimatedOutput(Reg) << " ";
		cout << "Here is output: " << Y[i] << endl;
		Test_Block.update_par(Reg, Y[i]); // Update parameters in respect to Input and Output
		cout << "Here are the estimated parameters at time " << i << " : " << endl;
		Test_Block.getEstimatedParameters().print(); //Print Matrix 
		cout << '\n';
		outFileBlock1 << Test_Block.getEstimatedParameters().row(0) << " ";
		outFileBlock2 << Test_Block.getEstimatedParameters().row(1) << " ";
		outFileBlock3 << Test_Block.getEstimatedParameters().row(2) << " ";
	}

	outFileBlock1.close();
	outFileBlock2.close();
	outFileBlock3.close();
	outFile_out_Block.close();

	return 0;
}