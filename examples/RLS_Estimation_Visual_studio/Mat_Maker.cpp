/*
#include <cstdlib>
#include <random>
#include <fstream>
using namespace std;

int main() {
	//Initializing parameters//
	const int len = 500; //Length of Array of Output/Input
	//Generate Data for Validation
	double Y[len] = { 0. };
	double X[len] = { 0. };

	double counter = 0.;

	ofstream outFileOut;
	ofstream outFileIn;
	outFileOut.open("C:/Users/Nikos/rls/testing/exp_testing/TXT-Files/PolyRLS/Test_Output.txt");
	outFileIn.open("C:/Users/Nikos/rls/testing/exp_testing/TXT-Files/PolyRLS/Test_Input.txt");

	for (int i = 0; i < 100; i++) {
		Y[i] = counter + 10; //OUTPUTS OF SYSTEM
		X[i] = counter; //INPUTS OF SYSTEM
		outFileOut << Y[i] << " ";
		outFileIn << X[i] << " ";
		counter += 1.;
	}
	for (int i = 100; i < 200; i++) {
		Y[i] = 50.; //OUTPUTS OF SYSTEM
		X[i] = counter; //INPUTS OF SYSTEM
		outFileOut << Y[i] << " ";
		outFileIn << X[i] << " ";
		counter += 1.;
	}

	for (int i = 200; i < 300; i++) {
		Y[i] = counter - 50.; //OUTPUTS OF SYSTEM
		X[i] = counter; //INPUTS OF SYSTEM
		outFileOut << Y[i] << " ";
		outFileIn << X[i] << " ";
		counter += 1.;
	}

	for (int i = 300; i < 400; i++) {
		Y[i] = 20.; //OUTPUTS OF SYSTEM
		X[i] = counter; //INPUTS OF SYSTEM
		outFileOut << Y[i] << " ";
		outFileIn << X[i] << " ";
		counter += 1.;
	}

	for (int i = 400; i < 500; i++) {
		Y[i] = 0.4*(i-398.)*(i-398.) - 2.*(i-398.) + 1.; //OUTPUTS OF SYSTEM
		X[i] = counter; //INPUTS OF SYSTEM
		outFileOut << Y[i] << " ";
		outFileIn << X[i] << " ";
		counter += 1.;
	}
	outFileOut.close();
	outFileIn.close();
}*/