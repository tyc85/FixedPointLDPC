#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include "ArrayLDPCMacro.h"
#include "ArrayLDPC.h"
#include "PerfTest.h"
#include <windows.h>



using std::endl;
using std::cout;
using std::cin;

int main(int argc, char *argv[])
{
	//char Filename[30];
	double db_start, db_end, db_step;
	double Rate, EbN0, snr;
	int i = 0;
	if(argc == 5)
	{
		//printf("%.1f, %.1f, %.1f, %s\n", atof(argv[1]),atof(argv[2]),atof(argv[3]),argv[4]);
		//ArrayLDPC_PerfTest(atof(argv[1]), atof(argv[2]), atof(argv[3]), argv[4]);
		ArrayLDPC_PerfTest(2, 2, 1, "test.csv");
	}
	else if(1) // debugging
	{
		//ArrayLDPC_Debug_Shorten(36);
		ArrayLDPC_Debug_Wifi();
		return 0;
		long MaxPckNum;
		printf("argc = %d \n time trial with: manual input\n", argc);
		cout << "EbN0 in dB: ";
		cin >> EbN0;
		cout << "simulate how many packets? ( " << 2209 << " is the codeword length ): ";
		cin >> MaxPckNum;

		char InfoStream[248] = "OMG how long should this string be to make it 248, just imagine that. \
								   I guess it's still not long enough. Let's see. This is a testing string \
									for a lot of characters so that we have some random bit stream that's\
									correct";
		//---- code timging for windows

		//---- encoder trial
		//EncodeTrial(InfoStream, MaxPckNum);
		//---- decoder trial
 		//DecodeTrial(EbN0, MaxPckNum);
		//---- test error rate
		//ArrayLDPC_Debug();

		
	}
	else
	{
		long MaxPckNum;
		printf("argc = %d \n time trial with: manual input\n", argc);
		cout << "EbN0 in dB: ";
		cin >> db_start;
		cout << "simulate how many packets? ( " << 2209 << " is the codeword length ): ";
		cin >> MaxPckNum;

		//---- code timging for windows
		//---- should create a linux version too
		__int64 ctr1 = 0, ctr2 = 0, freq = 0;
      // Start timing the code.
      if(QueryPerformanceCounter((LARGE_INTEGER *) &ctr1) != 0)
		{
			// Do what ever you do, what ever you need to time...
			
			// Finish timing the code.
			ArrayLDPC_TimeTrial(db_start, MaxPckNum, "timing.txt");
         QueryPerformanceCounter((LARGE_INTEGER *) &ctr2);
         QueryPerformanceFrequency((LARGE_INTEGER *) &freq);

         // Print the time spent in microseconds to the console.

         std::cout << ((ctr2 - ctr1) * 1.0 / freq) << "  seconds" << std::endl;
			std::cout << 2209*MaxPckNum/((ctr2 - ctr1) * 1.0 / freq) << " bits per second" << std::endl;
      }
		//cout << "csv filname:";
		//cin >> Filename;
		
	}
	//else
	//{
	//	//cout << "snr from:";
	//	//cin >> db_start;
	//	//cout << "snr end:";
	//	//cin >> db_end;
	//	//cout << "snr step:";
	//	//cin >> db_step;
	//	//cout << "csv filname:";
	//	//cin >> Filename;
	//	ArrayLDPC_PerfTest(1, 2, 1, "test.csv");
	//}
	system("pause");
	return 0;
}

//---------- PerfTest Main
//int main(int argc, char *argv[])
//{
//	//char Filename[30];
//	double db_start, db_end, db_step;
//	int i = 0;
//	if(argc == 5)
//	{
//		printf("%.1f, %.1f, %.1f, %s\n", atof(argv[1]),atof(argv[2]),atof(argv[3]),argv[4]);
//		ArrayLDPC_PerfTest(atof(argv[1]), atof(argv[2]), atof(argv[3]), argv[4]);
//	}
//	else
//	{
//		long MaxPckNum;
//		printf("argc = %d \n time trial with: manual input\n", argc);
//		cout << "EbN0 in dB: ";
//		cin >> db_start;
//		cout << "simulate how many packets? ( " << 2209 << " is the codeword length ): ";
//		cin >> MaxPckNum;
//
//		//---- code timging for windows
//		//---- should create a linux version too
//		__int64 ctr1 = 0, ctr2 = 0, freq = 0;
//      // Start timing the code.
//      if(QueryPerformanceCounter((LARGE_INTEGER *) &ctr1) != 0)
//		{
//			// Do what ever you do, what ever you need to time...
//			
//			// Finish timing the code.
//			ArrayLDPC_TimeTrial(db_start, MaxPckNum, "timing.txt");
//         QueryPerformanceCounter((LARGE_INTEGER *) &ctr2);
//         QueryPerformanceFrequency((LARGE_INTEGER *) &freq);
//
//         // Print the time spent in microseconds to the console.
//
//         std::cout << ((ctr2 - ctr1) * 1.0 / freq) << "  seconds" << std::endl;
//			std::cout << 2209*MaxPckNum/((ctr2 - ctr1) * 1.0 / freq) << " bits per second" << std::endl;
//      }
//		//cout << "csv filname:";
//		//cin >> Filename;
//		
//	}
//	//else
//	//{
//	//	//cout << "snr from:";
//	//	//cin >> db_start;
//	//	//cout << "snr end:";
//	//	//cin >> db_end;
//	//	//cout << "snr step:";
//	//	//cin >> db_step;
//	//	//cout << "csv filname:";
//	//	//cin >> Filename;
//	//	ArrayLDPC_PerfTest(1, 2, 1, "test.csv");
//	//}
//	//system("pause");
//	return 0;
//}