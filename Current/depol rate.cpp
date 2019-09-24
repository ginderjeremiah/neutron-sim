#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
//#include <string>
using namespace std;

int main()
{
	ifstream inFile;
	long double temp,
				temp2,
				depolSum = 0;
	string tempString;
	stringstream stream;
	int numRuns = 10;
	
	for (int i = 0; i < numRuns; i++)
	{
		//stream.clear();
		stream.str("");
		stream << "spindata " << i << " symN.txt";
		inFile.open(stream.str().c_str());

		inFile >> tempString;
		
		while (!inFile.eof())
		{
			inFile >> temp >> temp >> temp >> temp >> temp2 >> temp >> tempString;
			//cout << i << " " << temp << endl;
			depolSum += temp2;
		}
		inFile.close();
	}
	
	cout << depolSum / (numRuns * 10) << endl;
	
	return 0;
}