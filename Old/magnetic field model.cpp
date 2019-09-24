#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
using namespace std;

double magnitude(double array[]);
const double PI = atan(1) * 4;
const double GYRO_RATIO = 1.832 * (pow(10, 8)); // neutron gyromagnetic ratio
const double HBAR = 6.582 * (pow(10, -16)); // in Ev*s
const double MAG_THICK = 0.0254; //thickness of layer of PM array
const double MAG_SPACE = 0.0508; //characteristic spacing of magnets
const double B_REM = 1.4; //remnant magnet field strength
const double N_TERMS = 5; //number of terms used from infinite series for b field
const double GRAV = 9.80665; //gravitational constant
const double MASS_N = 1.674927471 * (pow(10, -27)); //mass of neutron
const double MU_N = -9.6623647 * (pow(10, -27)); //magnetic moment of neutron
const double JTONEV = 6.241509 * (pow(10, 27)); //joules to neV

int main()
{
	double magField[3] = {0, 0, 0};
	double position[3] = {0, 0, 0.001};
	double k_n;
	ofstream outFile;
	
	outFile.open("spindata.txt");
	for (position[0] = 0; position[0] <= 0.05; position[0] += 0.00001)
	{
		for (int i = 0; i < 3; i++)
		{
			magField[i] = 0.0;
		}
		
		for (int n = 1; n <= N_TERMS; ++n)
		{
			k_n = 2 * PI * (4*n - 3) / MAG_SPACE;
			if (n % 2 == 1)
			{
				magField[0] += 1.0/(4*n - 3) * (1-exp(-k_n*MAG_THICK))*exp(-k_n*position[2])*sin(k_n*position[0]);
				magField[2] += 1.0/(4*n - 3) * (1-exp(-k_n*MAG_THICK))*exp(-k_n*position[2])*cos(k_n*position[0]);
			}
		}
		for (int i = 0; i < 3; i++)
		{
			magField[i] *= 4 * B_REM / (PI * sqrt(2));
		}
		
		outFile << position[0] << " " << magnitude(magField) << endl;
	}
	return 0;
}

double magnitude(double array[])
{
	return sqrt(array[0] * array[0] + array[1] * array[1] + array[2] * array[2]);
}