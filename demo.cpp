#include <boost/math/tools/roots.hpp>
#include "cantera/thermo.h"
#include <iostream>
#include<string.h>
#include<fstream>
#include<math.h>

using namespace Cantera;
using namespace std;

/*calculation of redlich constant a */
double constantCal_a(double Tc, double Pc, double R) {
	double a;
	a = (0.42748 * (pow(R, 2) * pow(Tc, 2))) / Pc;
	return a;
}

/*calculation of redlich constant b */
double constantCal_b(double Tc, double Pc, double R) {
	double b;
	b = (0.08664 * (R * Tc)) / Pc;
	return b;
}

/*Partial derivative of dv with respect to dT derived from the original Redlich-Kwong equation*/
double dv_dT(double Tc, double Pc, double R, double T, double P, double V) {
	double Tr = T / Tc;
	double a = constantCal_a(Tc, Pc, R);
	double b = constantCal_b(Tc, Pc, R);
	double num = (1.5 * R * V * (V + b) * (pow(Tr, 0.5))) - (P * 0.5 * (pow(V, 3) - (pow(b, 2) * V)) / ((pow(Tc, 0.5)) * (pow(T, 0.5))));
	double den = ((pow(Tr, 0.5) * P) * ((3 * pow(V, 2)) - pow(b, 2))) + a - ((R * pow(T, 1.5) / (pow(Tc, 0.5))) * ((2 * V) + b));
	return (num / den);
}

/*Partial derivative of dv with respect to dP derived from the original Redlich-Kwong equation*/
double dv_dP(double Tc, double Pc, double R, double T, double P, double V) {
	double Tr = T / Tc;
	double a = constantCal_a(Tc, Pc, R);
	double b = constantCal_b(Tc, Pc, R);
	double num = 1;
	double den = (-R * T / (pow(V - b, 2))) + ((a / (pow((V * V) + (b * V), 2) * pow(Tr, 0.5))) * ((2 * V) + b));
	return (num / den);
}

/*calculation of Cv used in Redlich Kwong EOS*/
double CvCal_Redlich(double Tc, double Pc, double R, double T, double P, double V, double Cp, double a, double b) {
	double Tr = T / Tc;
	double Z = (P * V) / (R * T);
	double dV_dP = dv_dP(Tc, Pc, R, T, P, V);
	double dV_dT = dv_dT(Tc, Pc, R, T, P, V);

	double dZ_dT = dV_dT * ((-b / (pow((V - b), 2))) + ((a * pow(Tc, 0.5)) / (pow(T, 1.5) * R * (pow(V + b, 2))))) + (1.5 * a * pow(Tc, 0.5) / (R * (V + b) * pow(T, 2.5)));
	double dZ_dP = ((V / (V - b)) - (a / ((V + b) * pow(Tr, 0.5) * R * T)) - Z + (dV_dP * (-P * b / (pow(V - b, 2)) + (a * P / (R * T * pow(Tr, 0.5) * pow(V + b, 2)))))) / P;
	double Cv = Cp - ((R * pow((Z + (T * dZ_dT)), 2)) / (Z - (P * dZ_dP)));
	return Cv;
}

/*calculation of Cp used in Redlich Kwong EOS*/
double CpCal_Redlich(double Tc, double Pc, double R, double T, double P, double V) {
	double a = constantCal_a(Tc, Pc, R);
	double b = constantCal_b(Tc, Pc, R);
	double Tr = T / Tc;
	double dV = dv_dT(Tc, Pc, R, T, P, V);
	double Cp = (P * dV) - R + (((a / b) * (0.75) * (log((V + b) / V) * (pow(Tc, 0.5)))) / pow(T, 1.5)) + (a * 1.5 * (pow(Tr, -0.5) * dV) / ((V + b) * V));
	return Cp;
}

/*function which calculates the value of enthalpy equation for a given temperature*/
double functionValue(double T,double R) {
	double H_O2, H_N2, H_AR;
	double H_air;
	double Xo2, Xn2, Xar;
	
	Xo2 = 0.21;
	Xn2 = 0.78;
	Xar = 0.01;
	
	/* H(T) = (a0 + (a1/2)T + (a2/3)T^2 + (a3/4)T^3 + (a4/5)T^4 + (a5/T))RT */
	H_O2 = (3.78245636 - (((2.99673416e-03) / 2) * (T)) + (((9.84730201e-06) / 3) * (pow(T, 2))) - (((9.68129509e-09) / 4) * (pow(T, 3))) + ((3.24372837e-12 / 5) * (pow(T, 4))) - (1063.94356 / T)) * R * T;
	H_N2 = (3.298677 + (((1.4082404e-03) / 2) * (T)) - (((3.963222e-06) / 3) * (pow(T, 2))) + (((5.641515e-09) / 4) * (pow(T, 3))) - (((2.444854e-12) / 5) * (pow(T, 4))) - (1020.8999) / T) * R * T;
	H_AR = (2.5 - (745.375 / T)) * R * T;

	/*final enthalpy equation as a fucntion of temperature of the gas*/
	H_air = (Xo2 * H_O2) + (Xn2 * H_N2) + (Xar * H_AR);
	return H_air;
}

/*function which calcualtes the derivative value of enthalpy equation for a given temperature*/
double firstDer(double T,double R) {
	double dH_O2, dH_N2, dH_AR;
	double dH_air, fun_val;
	double Xo2, Xn2, Xar;

	//composition of gases in air in molefractions
	Xo2 = 0.21;
	Xn2 = 0.78;
	Xar = 0.01;

	/* dH(T)/dT = a0 + (a1)T + (a2)T^2 + (a3)T^3 + (a4)T^4 */
	dH_O2 = (3.78245636 - ((2.99673416e-03) * (T)) + ((9.84730201e-06) * (pow(T, 2))) - ((9.68129509e-09) * (pow(T, 3))) + ((3.24372837e-12) * (pow(T, 4)))) * R;
	dH_N2 = (3.298677 + ((1.4082404e-03) * (T)) - ((3.963222e-06) * (pow(T, 2))) + ((5.641515e-09) * (pow(T, 3))) - ((2.444854e-12) * (pow(T, 4)))) * R;
	dH_AR = (2.5) * R;

	/*final enthalpy equation as a fucntion of temperature of the gas*/
	dH_air = (Xo2 * dH_O2) + (Xn2 * dH_N2) + (Xar * dH_AR);
	return fun_val = dH_air;
}

/*function which implements Newton Raphson method*/
double newtonRaphson() {
	auto sol = newSolution("air.yaml");
	auto gas = sol->thermo();
	int iterNo = 0; //variable to store Number of iterations
	double x = 527; //initial guess value of temperature
	double xPre = 0; //variable to store x previous value
	double error = 1; //variable to store error 
	double const tol = 1e-5;//variable to store tolerance 
	cout << endl << "Newton raphson method implementation" << endl;
	//when error is less than tolerance loop terminates
	while (error > tol || iterNo < 1) {
		iterNo += 1; //variable to count number of iterations
		xPre = x; //variable to store previous value of root
		double R = gas->RT() / (x * gas->meanMolecularWeight());

		//Xn+1 ​= Xn​ − (f(x)/f′(x))​.
		x = x - (functionValue(x,R) / firstDer(x,R));
		cout << "R value in Newton method  = " << R << "  when temp = " << x << endl;
		error = abs(x - xPre);  // error = |x-xpre|
		cout << "iteration No = " << iterNo << "  current root = " << x << "  previous root = " << xPre << "  error = " << error << endl;
	}
	return x;
}

/*Calcualtion of Cp with the formula used in NASA7 co-efficients method */
double Cp_cal(double T,double R) {
	double Cp_O2,Cp_N2, Cp_AR;
	double Cp_air;
	double Xo2, Xn2, Xar;

	//composition of gases in air in molefractions
	Xo2 = 0.21;
	Xn2 = 0.78;
	Xar = 0.01;

	/* Cp(T)/R = a0 + a1T + a2T^2 + a3T^3 + a4T^4  */
	Cp_O2 = (3.78245636 - ((2.99673416e-03) * (T)) + ((9.84730201e-06) * (pow(T, 2))) - ((9.68129509e-09) * (pow(T, 3))) + ((3.24372837e-12) * (pow(T, 4))))*R;
	Cp_N2 = (3.298677 + ((1.4082404e-03) * (T)) - ((3.963222e-06) * (pow(T, 2))) + ((5.641515e-09) * (pow(T, 3))) - ((2.444854e-12) * (pow(T, 4))))*R;
	Cp_AR = 2.5*R;

	// Cp = ∑(Xi*Cpi)
	Cp_air = (Xo2 * Cp_O2 + (Xn2 * Cp_N2) + (Xar*Cp_AR));
	return Cp_air;
}

/*Calcualtion of enthalpy with the formula used in NASA7 co-efficients method */
double H_cal(double T, double R) {
	double H_O2, H_N2, H_AR;
	double H_air;
	double Xo2, Xn2, Xar;

	//composition of gases in air in molefractions
	Xo2 = 0.21;
	Xn2 = 0.78;
	Xar = 0.01;

	/* H(T)/RT = a0 + (a1/2)T + (a2/3)T^2 + (a3/4)T^3 + (a4/5)T^4 + (a5/T) */
	H_O2 = (3.78245636 - (((2.99673416e-03) / 2) * (T)) + (((9.84730201e-06) / 3) * (pow(T, 2))) - (((9.68129509e-09) / 4) * (pow(T, 3))) + ((3.24372837e-12 / 5) * (pow(T, 4))) - (1063.94356 / T)) * R * T;
	H_N2 = (3.298677 + (((1.4082404e-03) / 2) * (T)) - (((3.963222e-06) / 3) * (pow(T, 2))) + (((5.641515e-09) / 4) * (pow(T, 3))) - (((2.444854e-12) / 5) * (pow(T, 4))) - (1020.8999) / T) * R * T;
	H_AR = (2.5 - (745.375 / T)) * R * T;

	//RT() * mean_X(enthalpy_RT_ref());
	// H = ∑(Xi*Hi)
	H_air = (Xo2 * H_O2) + (Xn2 * H_N2) + (Xar * H_AR);
	return H_air; //return the calculated enthalpy of air
}

/*Calcualtion of entropy with the formula used in NASA7 co-efficients method */
double S_cal(double T, double R) {
	double S_O2, S_N2, S_AR;
	double S_air;
	double Xo2, Xn2, Xar;
	//composition of gases in air in molefractions
	Xo2 = 0.21;
	Xn2 = 0.78;
	Xar = 0.01;

	/* S(T)/R = a0lnT + a1T + (a2/2)T^2 + (a3/3)T^3 + (a4/4)T^4 + a6 */
	S_O2 = ((3.78245636 * log(T)) - ((2.99673416e-03) * (T)) + (((9.84730201e-06) / 2) * (pow(T, 2))) - (((9.68129509e-09) / 3) * (pow(T, 3))) + (((3.24372837e-12) / 4) * (pow(T, 4))) + (3.65767573)) * R;
	S_N2 = ((3.298677 * log(T)) + ((1.4082404e-03) * (T)) - (((3.963222e-06) / 2) * (pow(T, 2))) + (((5.641515e-09) / 3) * (pow(T, 3))) - (((2.444854e-12) / 4) * (pow(T, 4))) + (3.950372)) * R;
	S_AR = ((2.5 * log(T)) + (4.366)) * R;

	//Small number to compare differences of mole fractions against.
	double Tiny = 1.e-20;

	//∑XlogX
	double Sum_xlogx = (Xo2 * log(Xo2 + Tiny)) + (Xn2 * log(Xn2 + Tiny)) + (Xar * log(Xar + Tiny));

	//S = ∑(Xi * Si) − Rlog(P / P0) - R∑XlogX;
	S_air = (Xo2 * S_O2) + (Xn2 * S_N2) + (Xar * S_AR) - (R*log(538657.7/OneAtm)) - (R*Sum_xlogx);
	return S_air; //return the calculated entropy of air

}

//function which calculates thermodynamic properties of a gas with the EOS model specified
std::map<std::string, double> gastherm(double T, double P, double R_fluid, std::string model_eos)
{

	double C_sound, zcf, gamma, Cp;
	
	if (model_eos == "tc_ideal")
	{
		Cp = 1005; // Remove and take from input
		//h0 = h_ref + (Cp*(T0-T_ref));  Add these also if possiblr
		//s0 = s_ref + (Cp*log(T0/ T_ref)) - (R_fluid * log(P0 / P_ref)); Add these also if possible
		gamma = Cp / (Cp - R_fluid);    // Mayer's Relation only applicable in this case
		C_sound = sqrt(gamma * R_fluid * T);
		zcf = 1;   // Compressibility factor (unity for thermally ideal gas)
		cout << "cp = " << Cp << "  gamma = " << gamma << " C_sound  = " << C_sound << endl;
	}
	
	//if the model specified is cantera then this block is implemented
	else if (model_eos == "cantera") {
		double Cv, H, S, cp, cv, H_calc, S_calc;		

		zcf = 1;   // Compressibility factor (unity for thermally ideal gas)

		// Create a new phase
		auto sol = newSolution("air.yaml");
		auto gas = sol->thermo();

		//set the temperature and pressure of the gas
		gas->setState_TP(T,P);
		double R = gas->RT() / (T * gas->meanMolecularWeight());

		/*print out T,P,V values of gas to know it's state*/
		cout << "temp = " << gas->temperature() << "  pressure = " << gas->pressure() << "  R = " << R << endl;

		//gas->setState_SV(gas->entropy_mass(), gas->molarVolume() / gas->meanMolecularWeight());
		//gas->setState_SH(gas->entropy_mass(),gas->enthalpy_mass());

		//cout << "temp = " << gas->temperature() << "  pressure = " << gas->pressure() << "  volume = " << gas->molarVolume() << endl;
	
		/*Thermodynamic properties calculated using inbuilt functions in cantera*/
		Cp = gas->cp_mass(); //inbuilt function in cantera to calucalte cp
		Cv = gas->cv_mass(); //inbuilt function in cantera to calucalte cv
		H = gas->enthalpy_mass(); //inbuilt function in cantera to calucalte enthalpy
		S = gas->entropy_mass(); //inbuilt function in cantera to calucalte entropy
		gamma = Cp / Cv; 
		C_sound = sqrt(gamma * R_fluid * T);
		
		/*Thermodynamics properties calculated using NASA7 co-efficient method */
		cp = Cp_cal(T, R); //call the function Cp_cal
		cv = cp - R; //for an ideal gas cv = cp - R
		H_calc = H_cal(T, R); //call the function H_cal
		S_calc = S_cal(T, R); //call the function S_cal

		cout << endl << "cp value actually given by cantera = " << cp <<"  ,cp calucalted = " << Cp << "  ,gamma = " << gamma << endl;
		cout << "enthalpy value actually given by cantera = " << H << "  ,enthalpy calculated = " << H_calc << endl;
		cout << "entropy value actually given by cantera = " << S << "  ,entropy calculated = " << S_calc << endl;
		printf("Cv value actually given by cantera = %.10lf, Cv calucalted %.10lf \n ", Cv,cv);

		double T = newtonRaphson();
		cout << "R = " << R << endl;
		
	}
	//if the model specified is Redlich Kwong then this block is implemented
	else if (model_eos == "Redlich Kwong") {
		double Cp_air;
		double Cp_O2, Cp_N2, Cp_AR;
		double Xo2, Xn2, Xar;
		Xo2 = 0.21;
		Xn2 = 0.78;
		Xar = 0.01;
		Cp_O2 = CpCal_Redlich(154.6, 5050000, 259.8, T, P, 0.25469);
		Cp_N2 = CpCal_Redlich(126.2, 3390000, 296.8, T, P, 0.291169);
		Cp_AR = CpCal_Redlich(150.8, 4870000, 208.1, T, P, 0.204034);

		// Cp = ∑(Xi*Cpi)
		Cp_air = (Xo2 * Cp_O2) + (Xn2 * Cp_N2) + (Xar * Cp_AR);
		cout << "Cp of O2 = " << Cp_O2 << " ,Cp of N2 = " << Cp_N2 << " ,Cp of AR = " << Cp_AR << " ,Cp of Air = " << Cp_air << endl;

		double a_O2 = constantCal_a(154.6, 5050000, 259.8);
		double b_O2 = constantCal_b(154.6, 5050000, 259.8);
		double a_N2 = constantCal_a(126.2, 3390000, 296.8);
		double b_N2 = constantCal_b(126.2, 3390000, 296.8);
		double a_AR = constantCal_a(150.8, 4870000, 208.1);
		double b_AR = constantCal_b(150.8, 4870000, 208.1);

		// Molecular weight of mixture = ∑(Xi*Mi)
		double avg_M = (Xo2 * 0.031999) + (Xn2 * 0.028) + (Xar * 0.039948);

		// √a = ∑(Xi*√ai)    =>   a = (∑(Xi*√ai))^2
		double avg_a = pow((Xo2 * pow(a_O2, 0.5)) + (Xn2 * pow(a_N2, 0.5)) + (Xar * pow(a_AR, 0.5)), 2);

		// b = ∑(Xi*bi)
		double avg_b = (Xo2 * b_O2) + (Xn2 * b_N2) + (Xar * b_AR);

		//R of mixture  = R_universal / M_mix
		double avg_R = 8.3144626 / avg_M;

		//Tc = (a/b)(0.20268/R)
		double Tc = ((avg_a / avg_b) * 0.20268) / avg_R;

		//Pc = (0.08664 R Tc)/b
		double Pc = (0.08664 * avg_R * Tc) / avg_b;

		double V = 0.00815534 / avg_M;

		cout << "Tc = " << Tc << "   Pc = " << Pc << "volume = " << V << endl;

		double Cv_air = CvCal_Redlich(Tc, Pc, avg_R, T, P, V, Cp_air, avg_a, avg_b);
	}

	//return the calculated Cp, gamma, c_sound and Zcf values
	return { {"Cp",Cp},{"gamma",gamma},{"C_sound",C_sound},{"zcf",zcf} };
}


// the main program just calls function simple_demo within
// a 'try' block, and catches CanteraError exceptions that
// might be thrown
int main()
{
    try {
		gastherm(527.69, 538657.7,287.00264,"cantera");
	}
    catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
    }
}




