/*
 * Changes in test 1.3 vs. NEW INTEGRATION METHOD:
 * 
 * -changed dynamics integrator from rk45 to symplectic method
 * 
 */

#include <iostream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
using namespace std;

void createRotationMatrix(double rotationMatrix[][3], double magField[], double angleRotated);
void magFieldCalc(double position[], double magField[]);
void magForceCalc(double position[], double spin[], double force[]);
void crossProduct(double array1[], double array2[], double array3[]);
double dotProduct(double array1[], double array2[]);
double magnitude(double array[]);
double zetaCoordinate(double position[]);
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
const double GTT = 0.0001; // Gauss to Tesla
const double MAXTIME = 1.3; //time to simulate neutron (sim time not real time)
const double maxSpinStep = 1 * pow(10, -9);
const double minSpinStep = 1 * pow(10, -9);
const double eulerStep = 1 * pow(10, -6);
const double dynamicsStep = 1 * pow(10, -5);
const int INTEGRATION_MODE = 1; //set 0 for normal, 1 for spin normalization every step
const int GEOMETRY_MODE = 0; //set 0 for infinite plane, 1 for array
const int OUTPUT_MODE = 1; //set 0 lifetime, 1 for detailed, 2 for bounce location
const int SLURM_MODE = 0; //set 1 to enable randomized starting conditions from file
double holdField;
const double SYMTERM = pow(2, 1.0/3); //used often in symplectic coefficients

int main(int argCount, char** argv)
{
    double coefficients[5][6] = { {1.0/4, 1.0/4, 0, 0, 0, 0},  //for Runge-Kutta-Fehlberg method
                                  {3.0/8, 3.0/32, 9.0/32, 0, 0, 0},
                                  {12.0/13, 1932.0/2197, -7200.0/2197, 7296.0/2197, 0, 0},
                                  {1, 439.0/216, -8, 3680.0/513, -845.0/4104, 0},
                                  {0.5, -8.0/27, 2, -3544.0/2565, 1859.0/4104, -11.0/40} },
		   symCoefficients[2][4] = { {.5153528374311229364, -.085782019412973646, .4415830236164665242, .1288461583653841854},
								   {.1344961992774310892, -.2248198030794208058, .7563200005156682911, .3340036032863214255} },
           order5[6] = {16.0/135, 0, 6656.0/12825, 28561.0/56430, -9.0/50, 2.0/55},  //for 5th order solution
           order4[6] = {25.0/216, 0, 1408.0/2565, 2197.0/4104, -0.2, 0},  //for 4th order solution
           increment[6][3], //mini-steps
           spin[3],
           magField[3] = {0, 0, 0},
           spinCrossField[3],
           t = 0,
           timeGoal = 0,
           position[3] = {0, 0, .02},
           tempPosition[3], //used for calculation new B-field
		   tempPosition2[3],
           tempSpin[3],
		   tempSpin2[3], //used for calculating spinCrossField at different increments
           velocity[3] = {0, 0, 0},
           tempVelocity[3] = {0, 0, 0},
           nextStepVelocity[3] = {0, 0, 0},
           force[3] = {0, 0, 0},
           adaptiveStep = 0.0001, //starting stepSize for spin
           stepScalar[3] = {1, 1, 1},
           errorTolerance = 1 * (pow(10, -33)),
           //errorTolerance2 = 1 * (pow(10, -12)),
           approximation5[3], //5th order approximation
           approximation4[3], //4th order approximation
           fieldStrength,
           polarBeforeB,
           polarAfterB,
           totalDepolar = 0,
           netDepolar,
           vel,
           phi,
           theta,
           bounce[7] = {0, 0, 0, 0, 0, 0},
		   previousPolarization = 1.0,
		   polarization = 1,
		   outputZetaGoal = .5,
		   magDirection[3],
		   prevMagDirection[3] = {0, 0, 0},
		   outputTimeGoal = 0,
		   angleRotated,
		   larmorFreq,
		   rotationMatrix[3][3],
		   transverseDrop,
		   dropHeight,
		   transverseVel,
		   tempField[3],
		   determinant;
	long int numberOfSteps = 0,
	         numberOfSteps2 = 0,
                 numberOfBounces = 0,
                 totalBounces = 0;
	bool reset = false,
             isBouncing = false;
	int presetID;
	ofstream outFile;
	ifstream inFile;
        char fileName[30];
	string tempString;

	//equations solving for are dv/dt = F / MASS_N and dS/dt = (-1.832*10^8 s^-1 T^-1)*(S x B)
	presetID = atoi(argv[1]);
	dropHeight = atoi(argv[2]) / 100.0;
	transverseVel = atoi(argv[3]) / 100.0;
	transverseDrop = atoi(argv[4]) / 100.0;
    //INTEGRATION_MODE = atoi(argv[2]);
	if (SLURM_MODE)
	{
		inFile.open("./randomConditions.txt");
		for (int i = 0; i < presetID % 100; i++)
		{
			getline(inFile, tempString);
		}
		for (int i = 0; i < 3; i++)
		{
			inFile >> position[i] >> velocity[i];
		}
	}
	else
	{
		position[0] = 0;
		position[1] = transverseDrop;
		position[2] = -1.5 + dropHeight;
		velocity[0] = 0;
		velocity[1] = transverseVel;
		velocity[2] = 0;
		/*
		position[0] = 0;
		position[1] = 0;
		position[2] = -1.4;
		velocity[0] = 1;
		velocity[1] = 1;
		velocity[2] = 0;
		*/
		
	}

	strcpy(fileName, "./spindata ");
	strcat(fileName, argv[1]);
	strcat(fileName, " ");
	strcat(fileName, argv[2]);
	strcat(fileName, " ");
	strcat(fileName, argv[3]);
	strcat(fileName, " ");
	strcat(fileName, argv[4]);
	strcat(fileName, " symI.txt");
	outFile.open(fileName);
	
	outFile << scientific;
	outFile.precision(15);

	holdField = (presetID / 1000) * 0.4 * GTT;  //presetID of 1000 = .4 G or 0.00004 T
	
	while (((GEOMETRY_MODE == 1) ? zetaCoordinate(tempPosition2) : tempPosition2[2] + 1.5) > 0 && timeGoal < MAXTIME)
	{
		for (int i = 0; i < 3; ++i)
		{
			tempPosition2[i] = position[i];
		}
		
		if (numberOfSteps2 == 0)
		{
			magFieldCalc(tempPosition2, magField);
			fieldStrength = magnitude(magField);
			for (int l = 0; l < 3; ++l)
			{
				spin[l] = magField[l] * (HBAR / 2) / fieldStrength;
			}
		}

		if (numberOfSteps2 != 0)
		{
			for (int i = 0; i < 4; i++)
			{
				magForceCalc(position, spin, force);
				for (int j = 0; j < 3; j++)
				{
					velocity[j] += dynamicsStep * force[j] / MASS_N * symCoefficients[1][i];
				}
				for (int j = 0; j < 3; j++)
				{
					position[j] += velocity[j] * dynamicsStep * symCoefficients[0][i];
				}
			}
		}
	
		for (int i = 0; i < 3; ++i)
		{
			tempVelocity[i] = (position[i] - tempPosition2[i]) / dynamicsStep;
		}
		timeGoal += dynamicsStep;
		numberOfSteps2++;
		
		while((t < timeGoal || numberOfSteps2 == 0) && ((GEOMETRY_MODE == 1) ? zetaCoordinate(tempPosition2) : tempPosition2[2] + 1.5) > 0)
		{
			for (int i = 0; i < 3; ++i)
			{
				tempPosition2[i] += tempVelocity[i] * eulerStep;
			}
			
			magFieldCalc(tempPosition2, magField);
			
			for (int i = 0; i < 3; ++i)
			{
				magField[i] *= -eulerStep * GYRO_RATIO;
			}
			
			determinant = 1 + pow(magnitude(magField), 2);
			
			tempSpin[0] = spin[0]*(1 + pow(magField[0], 2))/determinant - spin[1]*(magField[0]*magField[1] + magField[2])/determinant + spin[2]*(magField[0]*magField[2] + magField[1])/determinant;
			tempSpin[1] = spin[0]*(magField[0]*magField[1] + magField[2])/determinant + spin[1]*(1 + pow(magField[1], 2))/determinant + spin[2]*(magField[1]*magField[2] - magField[0])/determinant;
			tempSpin[2] = spin[0]*(magField[0]*magField[2] - magField[1])/determinant - spin[1]*(magField[1]*magField[2] - magField[0])/determinant + spin[2]*(1 + pow(magField[2], 2))/determinant;
			
			
			t += eulerStep;
			numberOfSteps++;
			for (int i = 0; i < 3; ++i)
			{
				if (INTEGRATION_MODE == 0)
				{
					spin[i] = tempSpin[i];
				}
				else if (INTEGRATION_MODE == 1)
				{
					spin[i] = tempSpin[i] / magnitude(tempSpin) * (HBAR / 2);
				}
				else
				{
					// INTEGRATION_MODE 2 is currently unavailable
				}
			}
			if (OUTPUT_MODE == 2)
			{
				if (numberOfSteps % 100000 == 0)
				{
					cout << t << endl;
				}
				if (magnitude(force) / MASS_N > 50)
				{
					if (bounce[0] < magnitude(force) / MASS_N)
					{
						bounce[0] = magnitude(force) / MASS_N;
						bounce[1] = t;
						bounce[2] = tempPosition2[0];
						bounce[3] = tempPosition2[1];
						bounce[4] = tempPosition2[2] + 1.5;
						bounce[5] = MASS_N * JTONEV * (magnitude(velocity) * magnitude(velocity) / 2 + GRAV * (tempPosition2[2] + 1.5)) - JTONEV * GYRO_RATIO * dotProduct(spin, magField);
						bounce[6] = dotProduct(spin, magField) / magnitude(spin) / magnitude(magField);
					}
					isBouncing = true;
				}
				if (magnitude(force) / MASS_N < 10 && isBouncing == true)
				{
					outFile << bounce[1] << "    " << bounce[2]*100 << "       " << bounce[3]*100 << "      " << bounce[4]*100 << endl;
					numberOfBounces++;
					bounce[0] = 0;
					isBouncing = false;
				}
			}
			else if (OUTPUT_MODE)
			{
				//if ((tempPosition2[2] + 1.5) <= outputZetaGoal)
				//if (numberOfSteps % 1000 == 0)
				if (t >= outputTimeGoal)
				{
					//outputZetaGoal -= 0.0001;
					magFieldCalc(tempPosition2, magField);
					polarization = dotProduct(spin, magField) / magnitude(spin) / magnitude(magField);
					outFile << t << " " << tempPosition2[0] << " "  << tempPosition2[1] << " "  <<  1000*(tempPosition2[2] + 1.5) << " " << abs(polarization - previousPolarization) << " " << 1.0 - polarization << " " << magnitude(magField) << " " << adaptiveStep << endl;
					//outFile << t << " " << tempPosition2[0] << " "  << tempPosition2[1] << " "  <<  1000 * (tempPosition2[2] + 1.5) << " " << 1.0 - polarization << " " << magnitude(magField) << endl;
					previousPolarization = polarization;
					
					//outFile << presetID << " " << holdField << " " << t << " " << tempPosition2[0] << " "  << tempPosition2[1] << " "  <<  tempPosition2[2] + 1.5 << " " << dotProduct(spin, magField) / magnitude(spin) / magnitude(magField) << " " << MASS_N * JTONEV * (magnitude(velocity) * magnitude(velocity) / 2 + GRAV * (tempPosition2[2] + 1.5)) - JTONEV * GYRO_RATIO * dotProduct(spin, magField) <<  endl;
					
					
					outputTimeGoal += 0.0002;
					//magFieldCalc(tempPosition2, magDirection);
					//outFile << t << " " << tempPosition2[0] << " "  << tempPosition2[1] << " "  <<  tempPosition2[2] + 1.5 << " " << 1.0 - dotProduct(magDirection, prevMagDirection) / magnitude(magDirection) / magnitude(prevMagDirection) << " " << MASS_N * JTONEV * (magnitude(velocity) * magnitude(velocity) / 2 + GRAV * (tempPosition2[2] + 1.5)) - JTONEV * GYRO_RATIO * dotProduct(spin, magField) << " " << adaptiveStep << endl;
					//magFieldCalc(tempPosition2, prevMagDirection);
					
				}
			}
			else
			{
				if (numberOfSteps % 1000000000 == 0)
				{
					magFieldCalc(tempPosition2, magField);
					outFile << presetID << " " << holdField << " " << t << " " << tempPosition2[0] << " "  << tempPosition2[1] << " "  <<  tempPosition2[2] + 1.5 << " " << dotProduct(spin, magField) / magnitude(spin) / magnitude(magField) << " " << MASS_N * JTONEV * (magnitude(velocity) * magnitude(velocity) / 2 + GRAV * (tempPosition2[2] + 1.5)) - JTONEV * GYRO_RATIO * dotProduct(spin, magField) <<  endl;
				}
			}
			
			if (((GEOMETRY_MODE == 1) ? zetaCoordinate(tempPosition2) : tempPosition2[2] + 1.5) < 0)
			{
				magFieldCalc(tempPosition2, magField);
				//outFile << presetID << " " << holdField << " " << t << " " << tempPosition2[0] << " "  << tempPosition2[1] << " "  <<  tempPosition2[2] + 1.5 << " " << dotProduct(spin, magField) /magnitude(spin) / magnitude(magField) << " " << MASS_N * JTONEV * (magnitude(velocity) * magnitude(velocity) / 2 + GRAV * (tempPosition2[2] + 1.5)) - JTONEV * GYRO_RATIO * dotProduct(spin, magField) <<  " below array" <<  endl;
				outFile.close();
				return 0;
			}
		}
	}
	//magFieldCalc(tempPosition2, magField);
    //outFile << presetID << " " << holdField << " " << t << " " << tempPosition2[0] << " "  << tempPosition2[1] << " "  <<  tempPosition2[2] + 1.5 << " " << dotProduct(spin, magField) / magnitude(spin) / magnitude(magField) << " " << MASS_N * JTONEV * (magnitude(velocity) * magnitude(velocity) / 2 + GRAV * (tempPosition2[2] + 1.5)) - JTONEV * GYRO_RATIO * dotProduct(spin, magField) <<  endl;
    outFile.close();
	return 0;
}

void createRotationMatrix(double rotationMatrix[][3], double magField[], double angleRotated)
{
	double sintheta = sin(angleRotated),
		   costheta = cos(angleRotated),
		   unitAxis[3],
		   determinant;
		   
	for (int i = 0; i < 3; i++)
	{
		unitAxis[i] = magField[i] / magnitude(magField);
	}
	
	rotationMatrix[0][0] = costheta + unitAxis[0] * unitAxis[0] * (1 - costheta);
	rotationMatrix[0][1] = unitAxis[0] * unitAxis[1] * (1 - costheta) - unitAxis[2] * sintheta;
	rotationMatrix[0][2] = unitAxis[0] * unitAxis[2] * (1 - costheta) + unitAxis[1] * sintheta;
	rotationMatrix[1][0] = unitAxis[0] * unitAxis[1] * (1 - costheta) + unitAxis[2] * sintheta;
	rotationMatrix[1][1] = costheta + unitAxis[1] * unitAxis[1] * (1 - costheta);
	rotationMatrix[1][2] = unitAxis[1] * unitAxis[2] * (1 - costheta) - unitAxis[0] * sintheta;
	rotationMatrix[2][0] = unitAxis[0] * unitAxis[2] * (1 - costheta) - unitAxis[1] * sintheta;
	rotationMatrix[2][1] = unitAxis[1] * unitAxis[2] * (1 - costheta) + unitAxis[0] * sintheta;
	rotationMatrix[2][2] = costheta + unitAxis[2] * unitAxis[2] * (1 - costheta);
	
	//determinant = rotationMatrix[0][0] * (rotationMatrix[1][1] * rotationMatrix[2][2] - rotationMatrix[1][2] * rotationMatrix[2][1]) - rotationMatrix[0][1] * (rotationMatrix[1][0] * rotationMatrix[2][2] - rotationMatrix[1][2] * rotationMatrix[2][0]) + rotationMatrix[0][2] * (rotationMatrix[1][0] * rotationMatrix[2][1] - rotationMatrix[1][1] * rotationMatrix[2][0]);
	
	//cout << determinant << " " << angleRotated << endl;
	
	return;
}

void magFieldCalc(double position[], double magField[])
{ //Magnetic Flux Density in Teslas
	magField[0] = 0.0;
	magField[1] = 0.0;
	magField[2] = 0.0;
	double k_n,
		   toroidal[3], // xi, eta, zeta
		   theta,
		   phi,
		   rho,
		   R,
		   r,
		   sinphi,
		   cosphi,
		   sintheta,
		   costheta;
	
	if (GEOMETRY_MODE)
	{
		if (position[1] > 0.0)
		{
			R = 0.5;
			r = 1.0;
		}
		else
		{
			R = 1.0;
			r = 0.5;
		}
		
		rho = sqrt(pow(position[0], 2) + pow(position[2], 2));
		theta = atan2(position[0], -position[2]);
		phi = atan2(position[1], (rho - R));
		
		sinphi = sin(phi);
		cosphi = cos(phi);
		sintheta = sin(theta);
		costheta = cos(theta);
		
		//toroidal[0] = rho * theta;
		toroidal[1] = r * phi;
		toroidal[2] = r - sqrt(pow((rho - R), 2) + pow(position[1], 2));
		
		for (int n = 1; n <= N_TERMS; n += 1)
		{
			k_n = 2*PI*(4.0*n-3.0)/MAG_SPACE;
			magField[1] += (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*exp(-k_n*toroidal[2])*sin(k_n*toroidal[1]);
			magField[2] += (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*exp(-k_n*toroidal[2])*cos(k_n*toroidal[1]);
		}
		
		magField[0] = holdField;
		magField[1] *= 4 * B_REM / (PI * sqrt(2));
		magField[2] *= 4 * B_REM / (PI * sqrt(2));
		
		toroidal[0] = (-costheta) * magField[0] + (-sintheta * sinphi) * magField[1] +  (-sintheta * cosphi) * magField[2];
		toroidal[1] =                                      (cosphi) * magField[1] +             (-sinphi) * magField[2];
		toroidal[2] = (-sintheta) * magField[0] +  (costheta * sinphi) * magField[1] + (costheta * cosphi) * magField[2];

		for (int i = 0; i < 3; i++)
		{
			magField[i] = toroidal[i];
		}
	}
	else
	{
		for (int n = 1; n <= N_TERMS; n += 1)
		{
			k_n = 2*PI*(4.0*n-3.0)/MAG_SPACE;
			magField[1] += (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*exp(-k_n*(position[2]+1.5))*sin(k_n*position[1]);
			magField[2] += (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*exp(-k_n*(position[2]+1.5))*cos(k_n*position[1]);
		}
		
		magField[0] = holdField;
		magField[1] *= 4 * B_REM / (PI * sqrt(2));
		magField[2] *= 4 * B_REM / (PI * sqrt(2));
	}
	
	return;
}

void magForceCalc(double position[], double spin[], double force[])
{
	force[0] = 0.0;
	force[1] = 0.0;
	force[2] = 0.0;
	double k_n,
		   toroidal[3], // xi, eta, zeta
		   toroidalspin[3],
		   theta,
		   phi,
		   rho,
		   R,
		   r,
		   sinphi,
		   cosphi,
		   sintheta,
		   costheta;
	
	if (GEOMETRY_MODE)
	{
		if (position[1] > 0.0)
		{
			R = 0.5;
			r = 1.0;
		}
		else
		{
			R = 1.0;
			r = 0.5;
		}

		rho = sqrt(pow(position[0], 2) + pow(position[2], 2));
		theta = atan2(position[0], -position[2]);
		phi = atan2(position[1], (rho - R));

		//toroidal[0] = rho * theta;
		toroidal[1] = r * phi;
		toroidal[2] = r - sqrt(pow((rho - R), 2) + pow(position[1], 2));

		sinphi = sin(phi);
		cosphi = cos(phi);
		sintheta = sin(theta);
		costheta = cos(theta);
				
		toroidalspin[0] =          (-costheta) * spin[0] +                               (-sintheta) * spin[2];
		toroidalspin[1] = (-sintheta * sinphi) * spin[0] + (cosphi) * spin[1] + (costheta * sinphi) * spin[2];
		toroidalspin[2] =  (-sintheta * cosphi) * spin[0] + (-sinphi) * spin[1] + (costheta * cosphi) * spin[2];
				
		for (int n = 1; n <= N_TERMS; n += 1)
		{
			k_n = 2*PI*(4.0*n-3.0)/MAG_SPACE;
			force[1] += (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*(toroidalspin[1]*exp(-k_n*toroidal[2])*k_n*cos(k_n*toroidal[1])+toroidalspin[2]*exp(-k_n*toroidal[2])*-k_n*sin(k_n*toroidal[1]));
			force[2] += (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*(toroidalspin[1]*-k_n*exp(-k_n*toroidal[2])*sin(k_n*toroidal[1])+toroidalspin[2]*-k_n*exp(-k_n*toroidal[2])*cos(k_n*toroidal[1]));
		}
			
		force[1] *= GYRO_RATIO * 4 * B_REM / (PI * sqrt(2));
		force[2] *= GYRO_RATIO * 4 * B_REM / (PI * sqrt(2));

		toroidal[0] = (-costheta) * force[0] + (-sintheta * sinphi) * force[1] +  (-sintheta * cosphi) * force[2];
		toroidal[1] =                                      (cosphi) * force[1] +             (-sinphi) * force[2];
		toroidal[2] = (-sintheta) * force[0] +  (costheta * sinphi) * force[1] + (costheta * cosphi) * force[2] + MASS_N * (-GRAV);

		for (int i = 0; i < 3; i++)
		{
			force[i] = toroidal[i];
		}
	}
	else
	{
		for (int n = 1; n <= N_TERMS; n += 1)
		{
				k_n = 2*PI*(4.0*n-3.0)/MAG_SPACE;
				force[1] += (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*(spin[1]*exp(-k_n*(position[2]+1.5))*k_n*cos(k_n*position[1])+spin[2]*exp(-k_n*(position[2]+1.5))*-k_n*sin(k_n*position[1]));
				force[2] += (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*(spin[1]*-k_n*exp(-k_n*(position[2]+1.5))*sin(k_n*position[1])+spin[2]*-k_n*exp(-k_n*(position[2]+1.5))*cos(k_n*position[1]));
		}
		force[1] *= GYRO_RATIO * 4 * B_REM / (PI * sqrt(2));
		force[2] *= GYRO_RATIO * 4 * B_REM / (PI * sqrt(2));
		force[2] += MASS_N * (-GRAV);
	}
	
	return;
}

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

double zetaCoordinate(double position[])
{
    double r,
           R;
    
    if (position[1] > 0.0)
    {
        R = 0.5;
        r = 1.0;
    }
    else
    {
        R = 1.0;
        r = 0.5;
    }
    
    return r - sqrt(pow(((sqrt(pow(position[0], 2) + pow(position[2], 2))) - R), 2) + pow(position[1], 2));
}
