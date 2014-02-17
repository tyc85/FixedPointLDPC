#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <windows.h>
#include "rvgs.h"
#include "rngs.h"
#include "ArrayLDPCMacro.h"
#include "ArrayLDPC.h"
#include "PerfTest.h"



using std::endl;
using std::cout;
using std::cin;
using std::ofstream;
//---------------------------- Non in-class code
void noMoreMemory()
{
	cerr << "Unable to satisfy request for memory\n";
	abort();
}

int ArrayLDPC_Debug()
{
	char InfoStream[29] = "This is a testing string for";
	int InfoBit[INFO_LENGTH];
	int i, j;
	int CurrentInd = 0;
	int LLR_fp[CWD_LENGTH];
	//int PckError = 0;
	double BitError = 0;
	long Counter = 0;
	long Clk = 0;
	int Codeword[CWD_LENGTH];
	double Receive[CWD_LENGTH];
	double LLR[CWD_LENGTH];
	double snr;
	double db = 1.0;
	double Rate = 0;
	double EbN0_dB = 4;
	double EbN0 = 0;
	double Sigma;
	double BlkError = 0; 
	double PckError = 0;
	double **MemBank; //Memory bank
	bool per_flag;
	
	class FP_Decoder Decoder;
	class FP_Encoder Encoder("code\H_array_p47_r5g.txt", 0);
	class ControlFSM FSM;
	
	
	//cout << "SNR? ";
	//cin >> EbN0_dB;
	EbN0_dB = 5;
	Rate = Decoder.getRate();
	EbN0 = 2*pow(10.0, EbN0_dB/10)*Rate;			// Eb/N0	
	Sigma = sqrt(1/EbN0);
	BitError = 0;
	Clk = 0;
	while(PckError < 100)
	{
		BlkError = 0;
		Encoder.encode(InfoStream, Codeword, 231);
		for(i = 0; i < CWD_LENGTH; i++)
		{
			//-- All zero codeword
			// Using Box-Muller method to generate Gaussian noise
			LLR[i] = 2*EbN0*(1 + Normal(0, Sigma));
			// Using Wallace method to generate Gaussian noise
			//LLR[i] = 2*EbN0*(1 + Wallace(0, Sigma));
			LLR_fp[i] = int(LLR[i]*(1<<FRAC_WIDTH));
			//Debug
			//LLR_fp[i] = i % 47;
		}
		Decoder.setState(PCV);
		//BlkError = Decoder.decode(LLR);
		BlkError = Decoder.decode_fixpoint(LLR_fp);
		if(BlkError > 0)
			PckError++;
		BitError += BlkError;
		Counter++;
	}
	cout << BitError << " " << PckError << " " << Counter << endl
		<< " FER: " << PckError/Counter << " BER: " << BitError/Counter/CWD_LENGTH << endl;
	//Decoder.getInfoBit();
	return 0;
}


int ArrayLDPC_PerfTest(double db_start, double db_end, double db_step, char* Filename)
{
	char Filename2[40]="";
	int InfoBit[INFO_LENGTH];
	int i, j;
	int CurrentInd = 0;
	int LLR_fp[CWD_LENGTH];
	//int PckError = 0;
	double BitError = 0;
	long Counter = 0;
	long Clk = 0;

	double Receive[CWD_LENGTH];
	double LLR[CWD_LENGTH];
	double snr;
	double db = 1.0;
	double Rate = 0;
	double EbN0_dB = 4;
	double EbN0 = 0;
	double Sigma;
	double BlkError = 0; 
	double PckError = 0;
	double **MemBank; //Memory bank
	bool per_flag;

	class FP_Decoder Decoder;
	class ControlFSM FSM;

	ofstream FilePtr(Filename);
	if(!FilePtr)
	{
		cerr << "failed to open " << Filename << endl;
		system("pause");
		exit(0);
	}
	FilePtr.close();
	//strcat(Filename2, "log_");
	//strcat(Filename2, Filename);
	//strcat(Filename2, ".txt");
	sprintf(Filename2, "%s_log.txt", Filename);
	ofstream LogPtr(Filename2);
	if(!LogPtr)
	{
		cerr << "failed to open " << Filename2 << endl;
		system("pause");
		exit(0);
	}
	LogPtr.close();

	
	//cout << "SNR? ";
	//cin >> EbN0_dB;
	EbN0_dB = db_start;
	Rate = Decoder.getRate();
	EbN0 = 2*pow(10.0, EbN0_dB/10)*Rate;			// Eb/N0	
	Sigma = sqrt(1/EbN0);
	BitError = 0;
	Clk = 0;
	while(PckError < 100)
	{
		BlkError = 0;
		for(i = 0; i < CWD_LENGTH; i++)
		{
			//-- All zero codeword
			// Using Box-Muller method to generate Gaussian noise
			LLR[i] = 2*EbN0*(1 + Normal(0, Sigma));
			// Using Wallace method to generate Gaussian noise
			//LLR[i] = 2*EbN0*(1 + Wallace(0, Sigma));
			LLR_fp[i] = int(LLR[i]*(1<<FRAC_WIDTH));
			//Debug
			//LLR_fp[i] = i % 47;
		}
		Decoder.setState(PCV);
		//BlkError = Decoder.decode(LLR);
		BlkError = Decoder.decode_fixpoint(LLR_fp);
		if(BlkError > 0)
			PckError++;
		BitError += BlkError;
		Counter++;
	}
	cout << BitError << " " << PckError << " " << Counter << endl
		<< " FER: " << PckError/Counter << " BER: " << BitError/Counter/CWD_LENGTH << endl;
	//Decoder.getInfoBit();
	return 0;
}


int ArrayLDPC_TimeTrial(double db, int MaxPckNum, char* Filename)
{
	char Filename2[40]="";
	int InfoBit[INFO_LENGTH];
	int i, j;
	int CurrentInd = 0;
	int LLR_fp[CWD_LENGTH];
	//int PckError = 0;
	double BitError = 0;
	long Counter = 0;
	long Clk = 0;

	double Receive[CWD_LENGTH];
	double LLR[CWD_LENGTH];
	double snr;
	//double db = 1.0;
	double Rate = 0;
	double EbN0_dB = 4;
	double EbN0 = 0;
	double Sigma;
	double BlkError = 0; 
	double PckError = 0;
	double **MemBank; //Memory bank
	bool per_flag;

	class FP_Decoder Decoder;
	class ControlFSM FSM;

	ofstream FilePtr(Filename);
	if(!FilePtr)
	{
		cerr << "failed to open " << Filename << endl;
		system("pause");
		exit(0);
	}
	FilePtr.close();
	//strcat(Filename2, "log_");
	//strcat(Filename2, Filename);
	//strcat(Filename2, ".txt");
	sprintf(Filename2, "%s_log.txt", Filename);
	ofstream LogPtr(Filename2);
	if(!LogPtr)
	{
		cerr << "failed to open " << Filename2 << endl;
		system("pause");
		exit(0);
	}
	LogPtr.close();

	
	//cout << "SNR? ";
	//cin >> EbN0_dB;
	EbN0_dB = db;
	Rate = Decoder.getRate();
	EbN0 = 2*pow(10.0, EbN0_dB/10)*Rate;			// Eb/N0	
	Sigma = sqrt(1/EbN0);
	BitError = 0;
	Clk = 0;
	while(Counter < MaxPckNum)
	{
		BlkError = 0;
		for(i = 0; i < CWD_LENGTH; i++)
		{
			//-- All zero codeword
			// Using Box-Muller method to generate Gaussian noise
			LLR[i] = 2*EbN0*(1 + Normal(0, Sigma));
			// Using Wallace method to generate Gaussian noise
			//LLR[i] = 2*EbN0*(1 + Wallace(0, Sigma));
			LLR_fp[i] = int(LLR[i]*(1<<FRAC_WIDTH));
			//Debug
			//LLR_fp[i] = i % 47;
		}
		Decoder.setState(PCV);
		//BlkError = Decoder.decode(LLR);
		BlkError = Decoder.decode_fixpoint(LLR_fp);
		if(BlkError > 0)
			PckError++;
		BitError += BlkError;
		Counter++;
	}
	cout << BitError << " " << PckError << " " << Counter << endl
		<< " FER: " << PckError/Counter << " BER: " << BitError/Counter/CWD_LENGTH << endl;
	//Decoder.getInfoBit();
	
	return 0;
}


