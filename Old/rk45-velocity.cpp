#include <iostream>
#include <string>
#include <cmath>
using namespace std;

int main()
{
	double coefficients[5][6] = { {1.0/4, 1.0/4, 0, 0, 0, 0},  //coefficients for Runge-Kutta-Fehlberg method
							   {3.0/8, 3.0/32, 9.0/32, 0, 0, 0},
							   {12.0/13, 1932.0/2197, -7200.0/2197, 7296.0/2197, 0, 0},
							   {1, 439.0/216, -8, 3680.0/513, -845.0/4104, 0},
							   {0.5, -8.0/27, 2, -3544.0/2565, 1859.0/4104, -11.0/40} };
	double order5[6] = {16.0/135, 0, 6656.0/12825, 28561.0/56430, -9.0/50, 2.0/55};  //for 5th order solution
	double order4[6] = {25.0/216, 0, 1408.0/2565, 2197.0/4104, -0.2, 0};  //for 4th order solution
	double increment[6][3];
	double time = 0,
		   maxTime = 25,
		   position[3] = {4, 2, 3},
		   tempPosition[3],
//		   velocity[3] = {2, 1, 3},
		   adaptiveStep = 0.01,
		   stepScalar = 1,
//		   errorTolerance = .000000001,
		   errorTolerance = .001,
		   approximation5[3],
		   approximation4[3];
	int numberOfStepsTaken = 0;
	
	//equation solving for is dx/dt = 0.32x
	
	while (time < maxTime) //keep going until you reach the end
	{
		for (int i = 0; i < 3; ++i) //loop 3 times for x, y, and z dimensions
		{
			increment[0][i] = adaptiveStep * 0.32 * position[i]; // first increment based off orginal equation dx/dt = 0.32 * x
			
			tempPosition[i] = position[i] + increment[0][i] * coefficients[0][1]; //extrapolate new position based off of previous increment
			increment[1][i] = adaptiveStep * 0.32 * tempPosition[i]; //find second increment using new position
			
			tempPosition[i] = position[i] + increment[0][i] * coefficients[1][1] + increment[1][i] * coefficients[1][2]; // extrapolate using previous 2 increments... 
			increment[2][i] = adaptiveStep * 0.32 * tempPosition[i]; //find third increment using extrapolation...
			//the rest of this stuff is the same concept as previous lines
			tempPosition[i] = position[i] + increment[0][i] * coefficients[2][1] + increment[1][i] * coefficients[2][2] + increment[2][i] * coefficients[2][3];
			increment[3][i] = adaptiveStep * 0.32 * tempPosition[i];
			
			tempPosition[i] = position[i] + increment[0][i] * coefficients[3][1] + increment[1][i] * coefficients[3][2] + increment[2][i] * coefficients[3][3] + increment[3][i] * coefficients[3][4];
			increment[4][i] = adaptiveStep * 0.32 * tempPosition[i];
			
			tempPosition[i] = position[i] + increment[0][i] * coefficients[4][1] + increment[1][i] * coefficients[4][2] + increment[2][i] * coefficients[4][3] + increment[3][i] * coefficients[4][4] + increment[4][i] * coefficients[4][5];
			increment[5][i] = adaptiveStep * 0.32 * tempPosition[i];
		}
		// now that we have found all of our increments we can use them to help us find the solution
		for (int i = 0; i < 3; ++i)
		{
			approximation4[i] = position[i]; 
			approximation5[i] = position[i]; 
			for (int j = 0; j < 6; ++j)
			{
			approximation4[i] += increment[j][i] * (order4[j]); //as you can see we have two approximations, one 4th and one 5th order solution
			approximation5[i] += increment[j][i] * (order5[j]); //we'll be using the 4th order solution as our answer and using the difference between the 4th and 5th to adjust our step size to gain more speed or accuracy
			}
		}
		// and this next line is the formula used to find our optimal adjustment to the step size; right now it only uses the x-position but I think it can be adjusted to account for all of them in some fashion
		stepScalar = pow((errorTolerance * adaptiveStep / (2 * abs(approximation4[0] - approximation5[0]))), 0.25);
		time += adaptiveStep; //advance one step forward
		numberOfStepsTaken++;
		for (int i = 0; i < 3; ++i) //update position
		{
			position[i] = approximation4[i];
		}
		if (adaptiveStep * stepScalar < maxTime - time) //is the step size is larger than the remaining time left this will make sure it doesnt go past it
		{
			adaptiveStep *= stepScalar;
		}
		else
		{
			adaptiveStep = maxTime - time;
		}
	}
	
	cout << "Total steps taken: " << numberOfStepsTaken << endl;
	cout << "Final x-position: " << position[0] << endl;
	cout << "Final y-position: " << position[1] << endl;
	cout << "Final z-position: " << position[2] << endl;
	
}
		
		
		
	
	