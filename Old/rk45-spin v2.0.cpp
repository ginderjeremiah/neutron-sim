#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
using namespace std;

void magFieldCalc(double position[], double magField[]);
void magForceCalc(double position[], double spin[], double force[]);
void crossProduct(double array1[], double array2[], double array3[]);
double dotProduct(double array1[], double array2[]);
double magnitude(double array[]);
const double PI = atan(1) * 4;
const double GYRO_RATIO = -1.832 * (pow(10, 8)); // neutron gyromagnetic ratio
const double HBAR = 1.055 * (pow(10, -34)); // in J*s
const double MAG_THICK = 0.0254; //thickness of layer of PM array
const double MAG_SPACE = 0.0508; //characteristic spacing of magnets
const double B_REM = 1.4; //remnant magnet field strength
const double N_TERMS = 3; //number of terms used from infinite series for b field
const double GRAV = 9.80665; //gravitational constant
const double MASS_N = 1.674927471 * (pow(10, -27)); //mass of neutron
const double MU_N = -9.6623647 * (pow(10, -27)); //magnetic moment of neutron
const double JTONEV = 6.241509 * (pow(10, 27)); //joules to neV

int main()
{
	double coefficients[5][6] = { {1.0/4, 1.0/4, 0, 0, 0, 0},  //for Runge-Kutta-Fehlberg method
							   {3.0/8, 3.0/32, 9.0/32, 0, 0, 0},
							   {12.0/13, 1932.0/2197, -7200.0/2197, 7296.0/2197, 0, 0},
							   {1, 439.0/216, -8, 3680.0/513, -845.0/4104, 0},
							   {0.5, -8.0/27, 2, -3544.0/2565, 1859.0/4104, -11.0/40} };
	double order5[6] = {16.0/135, 0, 6656.0/12825, 28561.0/56430, -9.0/50, 2.0/55};  //for 5th order solution
	double order4[6] = {25.0/216, 0, 1408.0/2565, 2197.0/4104, -0.2, 0};  //for 4th order solution
	double increment[6][3]; //mini-steps
	//double spin[3] = {HBAR / 2 / sqrt(3), HBAR / 2 / sqrt(3), HBAR / 2 / sqrt(3)}; // has magnitude hbar / 2
	//double spin[3] = {-HBAR / 4, 0, -HBAR * sqrt(3) / 4};
	double spin[3];
	double magField[3] = {0, 0, 0};
	double spinCrossField[3];
	double time = 0;
	double maxTime = 5,
		   timeGoal = 0,
		   position[3] = {0, 0, .005},
		   tempPosition[3], //used for calculation new B-field
		   tempSpin[3], //used for calculating spinCrossField at different increments for equation
		   velocity[3] = {3, 3, 0},
		   tempVelocity[3] = {0, 0, 0},
		   nextStepVelocity[3] = {0, 0, 0},
		   force[3] = {0, 0, 0},
		   adaptiveStep = 0.0001, //starting stepSize for spin
		   adaptiveStep2 = 0.0001, //starting stepSize for position 
		   stepScalar[3] = {1, 1, 1},
		   errorTolerance = 1 * (pow(10, -37)),
		   errorTolerance2 = 5 * (pow(10, -12)),
		   approximation5[3], //5th order approximation
		   approximation4[3]; //4th order approximation
	long int numberOfSteps = 0,
	         numberOfSteps2 = 0;
	bool reset = false;
	double fieldStrength;
	ofstream outFile;
	string tempcin;
		   
	//equations solving for are dv/dt = F / MASS_N and dS/dt = (1.832*10^8 s^-1 T^-1)*(S x B)
	cout << "starting..." << endl;
	outFile.open("spindata.txt");
	while (time < maxTime)
	{
		for (int i = 0; i < 3; ++i)
		{
			tempVelocity[i] = velocity[i];
			tempPosition[i] = position[i];
		}
		
		if (numberOfSteps2 == 0)
		{
			magFieldCalc(tempPosition, magField);
			fieldStrength = magnitude(magField);
			for (int l = 0; l < 3; ++l)
			{
				spin[l] = magField[l] * (HBAR / 2) / fieldStrength;
			}
		}
		
		for (int i = 0; i < 5; ++i)
		{
			magForceCalc(tempPosition, spin, force);
			for (int j = 0; j < 3; ++j)
			{
				increment[i][j] = adaptiveStep2 * force[j] / MASS_N;
				tempPosition[j] = position[j] + tempVelocity[j] * adaptiveStep2 * coefficients[0][i];
				
				tempVelocity[j] = velocity[j];
				for (int k = 0; k <= i; ++k)
				{
					tempVelocity[k] += increment[k][j] * coefficients[i][k+1];
				}
			}
		}
		
		magForceCalc(tempPosition, spin, force);
		for (int i = 0; i < 3; ++i)
		{
			increment[5][i] = adaptiveStep2 * force[i] / MASS_N;
		}
		
		for (int i = 0; i < 3; ++i)
		{
			approximation4[i] = velocity[i];
			approximation5[i] = velocity[i];
			for (int j = 0; j < 6; ++j)
			{
			approximation4[i] += increment[j][i] * (order4[j]);
			approximation5[i] += increment[j][i] * (order5[j]);
			}
		}
		
		for(int i = 0; i < 3; ++i)
		{
			if (approximation4[i] - approximation5[i] == 0)
			{
				stepScalar[i] = 100000000;
			}
			else
			{
			stepScalar[i] = pow((errorTolerance2 * adaptiveStep2 / (2 * abs(approximation4[i] - approximation5[i]))), 0.25);
			}
		}
		if (stepScalar[0] > stepScalar[1])
		{
			stepScalar[0] = stepScalar[1];
		}
		if (stepScalar[0] > stepScalar[2])
		{
			stepScalar[0] = stepScalar[2];
		}
		if (stepScalar[0] > 2)
		{
			stepScalar[0] = 2;
		}
		if (numberOfSteps2 == 0)
		{
			numberOfSteps2++;
			adaptiveStep2 *= stepScalar[0];
			
			if (adaptiveStep2 > 1 * (pow(10, -4)))
			{
				adaptiveStep2 = 1 * (pow(10, -4));
			}
			
			if (adaptiveStep2 < 5 * (pow(10, -6)))
			{
				adaptiveStep2 = 5 * (pow(10, -6));
			}
			
			continue;
		}
		timeGoal += adaptiveStep;
		numberOfSteps2++;
		for (int i = 0; i < 3; ++i)
		{
			nextStepVelocity[i] = approximation4[i];
		}
		
		adaptiveStep2 *= stepScalar[0];
		
		if (adaptiveStep2 > 1 * (pow(10, -4)))
		{
			adaptiveStep2 = 1 * (pow(10, -4));
		}
		
		if (adaptiveStep2 < 5 * (pow(10, -6)))
		{
			adaptiveStep2 = 5 * (pow(10, -6));
		}
		if (adaptiveStep2 > maxTime - timeGoal)
		{
			adaptiveStep2 = maxTime - timeGoal;
		}
		
		while(time < timeGoal)
		{
			if (reset == false)
			{
				adaptiveStep = 0.0001;
			}
			//cout << time << " dumb" << endl << timeGoal << endl;
			magFieldCalc(position, magField);
			//Probably a way to put the following stuff into a loop instead of typing it a bunch, but my brain is tired
			crossProduct(spin, magField, spinCrossField);
			for (int i = 0; i < 3; ++i) //calculate increment 1 for each dimension and new position for increment 2
			{
				increment[0][i] = adaptiveStep * GYRO_RATIO * spinCrossField[i];
				tempPosition[i] = position[i] + velocity[i] * adaptiveStep * coefficients[0][0];
			}
			magFieldCalc(tempPosition, magField); //calculate magnetic field at new position
			for (int i = 0; i < 3; ++i) //calculate new spin vector
			{
				tempSpin[i] = spin[i] + increment[0][i] * coefficients[0][1];
			}
			crossProduct(tempSpin, magField, spinCrossField);
			for (int i = 0; i < 3; ++i) //calculate increment 2 for each dimension and new position for increment 3
			{
				increment[1][i] = adaptiveStep * GYRO_RATIO * spinCrossField[i];
				tempPosition[i] = position[i] + velocity[i] * adaptiveStep * coefficients[1][0];
			}
			magFieldCalc(tempPosition, magField); //calculate magnetic field at new position
			for (int i = 0; i < 3; ++i) //calculate new spin vector
			{
				tempSpin[i] = spin[i] + increment[0][i] * coefficients[1][1] + increment[1][i] * coefficients[1][2];
			}
			crossProduct(tempSpin, magField, spinCrossField);
			for (int i = 0; i < 3; ++i) //calculate increment 3 for each dimension and new position for increment 4
			{
				increment[2][i] = adaptiveStep * GYRO_RATIO * spinCrossField[i];
				tempPosition[i] = position[i] + velocity[i] * adaptiveStep * coefficients[2][0];
			}
			magFieldCalc(tempPosition, magField); //calculate magnetic field at new position
			for (int i = 0; i < 3; ++i) //calculate new spin vector
			{
				tempSpin[i] = spin[i] + increment[0][i] * coefficients[2][1] + increment[1][i] * coefficients[2][2] + increment[2][i] * coefficients[2][3];
			}
			crossProduct(tempSpin, magField, spinCrossField);
			for (int i = 0; i < 3; ++i) //calculate increment 4 for each dimension and new position for increment 5
			{
				increment[3][i] = adaptiveStep * GYRO_RATIO * spinCrossField[i];
				tempPosition[i] = position[i] + velocity[i] * adaptiveStep * coefficients[3][0];
			}
			magFieldCalc(tempPosition, magField); //calculate magnetic field at new position
			for (int i = 0; i < 3; ++i) //calculate new spin vector
			{
				tempSpin[i] = spin[i] + increment[0][i] * coefficients[3][1] + increment[1][i] * coefficients[3][2] + increment[2][i] * coefficients[3][3] + increment[3][i] * coefficients[3][4];
			}
			crossProduct(tempSpin, magField, spinCrossField);
			for (int i = 0; i < 3; ++i) //calculate increment 5 for each dimension and new position for increment 6
			{
				increment[4][i] = adaptiveStep * GYRO_RATIO * spinCrossField[i];
				tempPosition[i] = position[i] + velocity[i] * adaptiveStep * coefficients[4][0];
			}
			magFieldCalc(tempPosition, magField); //calculate magnetic field at new position
			for (int i = 0; i < 3; ++i) //calculate new spin vector
			{
				tempSpin[i] = spin[i] + increment[0][i] * coefficients[4][1] + increment[1][i] * coefficients[4][2] + increment[2][i] * coefficients[4][3] + increment[3][i] * coefficients[4][4] + increment[4][i] * coefficients[4][5];
			}
			crossProduct(tempSpin, magField, spinCrossField);
			for (int i = 0; i < 3; ++i) //calculate increment 6 for each dimension
			{
				increment[5][i] = adaptiveStep * GYRO_RATIO * spinCrossField[i];
			}
			//finally done with getting the increments now look at error in solution
			for (int i = 0; i < 3; ++i)
			{
				approximation4[i] = spin[i];
				approximation5[i] = spin[i];
				for (int j = 0; j < 6; ++j)
				{
				approximation4[i] += increment[j][i] * (order4[j]);
				approximation5[i] += increment[j][i] * (order5[j]);
				}
			}
			//find what the optimal step size should be
			for(int i = 0; i < 3; ++i)
			{
				if (approximation4[i] - approximation5[i] == 0)
				{
					stepScalar[i] = 100000000;
				}
				else
				{
				stepScalar[i] = pow((errorTolerance * adaptiveStep / (2 * abs(approximation4[i] - approximation5[i]))), 0.25);
				}
			}
			if (stepScalar[0] > stepScalar[1])
			{
				stepScalar[0] = stepScalar[1];
			}
			if (stepScalar[0] > stepScalar[2])
			{
				stepScalar[0] = stepScalar[2];
			}
			if (stepScalar[0] > 2)
			{
				stepScalar[0] = 2;
			}
			if (reset == false)
			{
				numberOfSteps++;
				adaptiveStep *= stepScalar[0];
				
				if (adaptiveStep > 1 * (pow(10, -8)))
				{
					adaptiveStep = 1 * (pow(10, -8));
				}
				if (adaptiveStep < 1 * (pow(10, -11)))
				{
					adaptiveStep = 1 * (pow(10, -11));
				}
				if (adaptiveStep > timeGoal - time)
				{
					adaptiveStep = timeGoal - time;
				}
				if (adaptiveStep2 < adaptiveStep && adaptiveStep2 != 0)
				{
					adaptiveStep = adaptiveStep2;
				}
				
				reset = true;
				continue;
			}
			//cout << adaptiveStep << " sd" << endl;
			/*
			if (time / maxTime > .8)
			{
				outFile << time << " " << timeGoal << " " << adaptiveStep << " " << adaptiveStep2 << endl;
			}
			*/
			if (numberOfSteps % 1000000 == 1)
			{
				magFieldCalc(position, magField);
				magForceCalc(position, spin, force);
				//outFile << time << " " << position[2] << " " << MASS_N * JTONEV * (GRAV * position[2] + (pow(magnitude(velocity), 2)) / 2) - dotProduct(spin, magField) * GYRO_RATIO * JTONEV << " " << dotProduct(spin, magField) * GYRO_RATIO * JTONEV << endl;
				outFile << time << " " << position[2] << " " << force[2] << " " << spin[2] << " " << magnitude(magField) << " " << dotProduct(spin, magField) / magnitude(spin) / magnitude(magField) << endl;
				if ((numberOfSteps) % 100000 == 1)
				{
					cout << "time and time goal: " << time << " : " << timeGoal << endl;
					cout << "spin magnitude ratio: " << magnitude(spin) / (HBAR / 2) << endl;
					//cout << "x-spin: " << spin[0] << endl;
					//cout << "y-spin: " << spin[1] << endl;
					//cout << "z-spin: " << spin[2] << endl;
					cout << "adaptivestep vel and spin: " << adaptiveStep2 << " : " << adaptiveStep << endl;
					cout << time / maxTime * 100 << "% complete" << endl;
				}
			}
			
			time += adaptiveStep;
			numberOfSteps++;
			for (int i = 0; i < 3; ++i)
			{
				spin[i] = approximation4[i];
				position[i] += adaptiveStep * velocity[i];
			}
			adaptiveStep *= stepScalar[0];
			
			if (adaptiveStep > 1 * (pow(10, -8)))
			{
				adaptiveStep = 1 * (pow(10, -8));
			}
			if (adaptiveStep < 1 * (pow(10, -11)))
			{
				adaptiveStep = 1 * (pow(10, -11));
			}
			if (adaptiveStep > timeGoal - time)
			{
				adaptiveStep = timeGoal - time;
			}
			if (adaptiveStep2 < adaptiveStep && adaptiveStep2 != 0)
			{
				adaptiveStep = adaptiveStep2;
			}
			/*
			if (spin[2] < 0)
			{
				cout << position[0];
				return 0;
			}
			*/
		}
		reset = false;
		for (int i = 0; i < 3; ++i)
		{
			velocity[i] = nextStepVelocity[i];
		}
	}
	
	cout << "Total steps taken: " << numberOfSteps + numberOfSteps2 << endl;
	cout << "Final x-spin: " << spin[0] << endl;
	cout << "Final y-spin: " << spin[1] << endl;
	cout << "Final z-spin: " << spin[2] << endl;
	
	return 0;
}
		
void magFieldCalc(double position[], double magField[])
{ //Magnetic Flux Density in Teslas
	magField[0] = 0.0;
	magField[1] = 0.0;
	magField[2] = 0.0;
	double k_n;
	
	for (int n = 1; n <= N_TERMS; n += 1)
		{
			k_n = 2*PI*(4.0*n-3.0)/MAG_SPACE;
			magField[0] += (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*exp(-k_n*position[2])*sin(k_n*position[0]);
			magField[2] += (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*exp(-k_n*position[2])*cos(k_n*position[0]);
		}
		
	magField[0] *= 4 * B_REM / (PI * sqrt(2));
	magField[1] = 0.01;
	magField[2] *= 4 * B_REM / (PI * sqrt(2));
	
	return;
}
void magForceCalc(double position[], double spin[], double force[])
{
	force[0] = 0.0;
	force[2] = 0.0;
	double xTerm,
		   zTerm,
		   xSum,
		   zSum,
		   xDerXSum,
		   xDerZSum,
		   zDerXSum,
		   zDerZSum;
	double k_n,
		   constant1,
		   constant2;
	
	for (int n = 1; n <= N_TERMS; n += 1)
	{
		k_n = 2*PI*(4.0*n-3.0)/MAG_SPACE;
		force[0] += (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*(spin[0]*exp(-k_n*position[2])*k_n*cos(k_n*position[0])+spin[2]*exp(-k_n*position[2])*-k_n*sin(k_n*position[0]));
		force[2] += (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*(spin[0]*-k_n*exp(-k_n*position[2])*sin(k_n*position[0])+spin[2]*-k_n*exp(-k_n*position[2])*cos(k_n*position[0]));
	}
		
	force[0] *= GYRO_RATIO * 4 * B_REM / (PI * sqrt(2));
	force[2] *= GYRO_RATIO * 4 * B_REM / (PI * sqrt(2));
	force[2] += MASS_N * (-GRAV);
	
	/*
	for (int n = 1; n <= N_TERMS; n += 1)
	{
		k_n = 2*PI*(4.0*n-3.0)/MAG_SPACE;
		constant1 = (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*exp(-k_n*position[2]);
		xTerm = constant1*sin(k_n*position[0]);
		zTerm = constant1*cos(k_n*position[0]);
		xSum += xTerm;
		zSum += zTerm;
		xDerXSum += zTerm*k_n;
		xDerZSum += xTerm*-k_n;
	}
	constant2 = 4 * B_REM / (PI * sqrt(2));
	xSum *= constant2;
	zSum *= constant2;
	xDerXSum *= constant2;
	xDerZSum *= constant2;
	zDerXSum = xDerZSum;
	zDerZSum = -xDerXSum;
	
	force[0] = MU_N * (xSum*xDerXSum + zSum*zDerXSum) / sqrt(xSum*xSum + zSum*zSum);
	force[2] = MU_N * (xSum*xDerZSum + zSum*zDerZSum) / sqrt(xSum*xSum + zSum*zSum);
	force[2] += MASS_N * (-GRAV);
	*/
	return;
}
/*
void magForceCalc(double position[], double spin[], double force[])
{
	force[0] = 0.0;
	force[2] = 0.0;
	double k_n;
	
	for (int n = 1; n <= N_TERMS; n += 1)
		{
			k_n = 2*PI*(4.0*n-3.0)/MAG_SPACE;
			force[0] += (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*(spin[0]*exp(-k_n*position[2])*k_n*cos(k_n*position[0])+spin[2]*exp(-k_n*position[2])*-k_n*sin(k_n*position[0]));
			force[2] += (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*(spin[0]*-k_n*exp(-k_n*position[2])*sin(k_n*position[0])+spin[2]*-k_n*exp(-k_n*position[2])*cos(k_n*position[0]));
		}
		
	force[0] *= GYRO_RATIO * 4 * B_REM / (PI * sqrt(2));
	force[2] *= GYRO_RATIO * 4 * B_REM / (PI * sqrt(2));
	//force[0] *= 4 * B_REM / (PI * sqrt(2));
	//force[2] *= 4 * B_REM / (PI * sqrt(2));
	force[2] += MASS_N * (-GRAV);
	
	return;
}
*/
void crossProduct(double array1[], double array2[], double array3[])
{
	array3[0] = array1[1] * array2[2] - array1[2] * array2[1];
	array3[1] = array1[2] * array2[0] - array1[0] * array2[2];
	array3[2] = array1[0] * array2[1] - array1[1] * array2[0];
	return;
}

double dotProduct(double array1[], double array2[])
{
	return array1[0]*array2[0] + array1[1]*array2[1] + array1[2]*array2[2];
}

double magnitude(double array[])
{
	return sqrt(array[0] * array[0] + array[1] * array[1] + array[2] * array[2]);
}