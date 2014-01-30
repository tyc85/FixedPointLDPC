#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include "ArrayLDPC.h"
using std::endl;
using std::cout;
using std::cin;

int main(int argc, char *argv[])
{
	char Filename[30];
	double db_start, db_end, db_step;
	int i = 0;
	if(argc == 5)
	{
		printf("%.1f, %.1f, %.1f, %s\n", atof(argv[1]),atof(argv[2]),atof(argv[3]),argv[4]);
		ArrayLDPC(atof(argv[1]), atof(argv[2]), atof(argv[3]), argv[4]);
	}
	else if(argc != 1)
	{
		printf("argument error, argc = %d \n manual input\n", argc);
		cout << "snr from:";
		cin >> db_start;
		cout << "snr end:";
		cin >> db_end;
		cout << "snr step:";
		cin >> db_step;
		cout << "csv filname:";
		cin >> Filename;
		ArrayLDPC(db_start, db_end, db_step, Filename);
	}
	else
	{
		//cout << "snr from:";
		//cin >> db_start;
		//cout << "snr end:";
		//cin >> db_end;
		//cout << "snr step:";
		//cin >> db_step;
		//cout << "csv filname:";
		//cin >> Filename;
		ArrayLDPC(1, 2, 1, "test.csv");
	}
	system("pause");
	return 0;
}