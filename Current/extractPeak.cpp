#include <iostream>
#include <cstring>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
using namespace std;

int main ()
{
	ofstream outFile;
	ifstream inFile;
	stringstream filestream;
    string fileName;
	string tempString;
	double tempNumber,
		   peak;
	
	outFile.open("./xd.txt");
	int j = 0;
	for (int i = 5; i <= 100; i += 5)
	{
		peak = 0;
		filestream << "./spindata 250000 25 " << i << " 0 symN.txt";
		fileName = filestream.str();
		filestream.clear();
		filestream.str("");
		inFile.open(fileName.c_str());
		while(!inFile.eof())
		{
			inFile >> tempNumber >> tempNumber >> tempNumber >> tempNumber >> tempNumber >> tempNumber;
			if (tempNumber > peak)
			{
				peak = tempNumber;
			}
			getline(inFile, tempString);
		}
		
		outFile << i << " " << peak << endl;
		
		inFile.close();
	}
	
	outFile.close();
	
	
	
	
	
	return 0;
}