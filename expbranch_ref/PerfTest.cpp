#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <bitset>
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


int ArrayLDPC_Debug_Wifi()
{
	//2209 - 231 = 1978
	//1978 / 8 = 247, 1978 % 8 = 2
	//char InfoStream[248] = "OMG how long should this string be to make it 248, just imagine that. \
								   I guess it's still not long enough. Let's see. This is a testing string \
									for a lot of characters so that we have some random bit stream that's\
									correct";
	//--- unused debugging varibles  
	//char InfoStream[243] = "";// all zero info, exactly 243 bytes for wifi code
	char InfoStream[122] = "OMG  how long   dd   should this string be to make it 243";
	//int InfoBit[INFO_LENGTH];
	//int Codewordtemp[CWD_LENGTH];
	//double Receive[CWD_LENGTH];
	//--- unused debugging varibles
	int i, j;
	int LLR_fp[CWD_LENGTH];
	//int pckerror = 0;
	double biterror = 0;
	long Counter = 0;

	double LLR[CWD_LENGTH];
	double db = 1.0;
	double Rate = 0;
	double EbN0_dB = 4;
	double snr = 0;
	double sigma;
	double blkerror = 0; 
	double pckerror = 0;
	int info_indx[INFO_LENGTH];
	
	class FP_Decoder Decoder;
	class FP_Encoder Encoder("H_802.11_IndZerog.txt", 0);
	
	EbN0_dB = 3;
	//Rate = Decoder.getRate();
	
	cout << "EbNo in dB? ";
	cin >> EbN0_dB;
	snr = 2*pow(10.0, EbN0_dB/10)*0.5;			// Eb/N0	to snr
	sigma = sqrt(1/snr);
	cout << "SNR is " << 10*log10(snr) << " dB" << endl;
	biterror = 0;

	//---- debugging
	//char out[277] = "";
	//out = new char[553];
	// working now
	//Encoder.encode(InfoStream, out, 248);
	//---- 

	Decoder.setInfoBit(InfoStream, 122);
	//for(i = 0; i < INFO_LENGTH; i++)
	//{
	//	//info_indx[i] = Encoder.getInfoIndex(i);
	//	// for wifi code it's systematic due to staircase structure
	//	info_indx[i] = i;
	//}
	// NEED IMPROVEMENT! SHOULD BE ABLE TO HAVE SYSTEMATIC BITS
	for(i = 0; i < INFO_LENGTH; i++)
	{
		info_indx[i] = Encoder.getInfoIndex(i);
	}
	//Decoder.setInfoIndex(info_indx);
	Decoder.setInfoIndex(info_indx);
	//for(i = 0; i < CWD_LENGTH; i++)
	//		Codewordtemp[i] = Encoder.getCodeword(i);
	//Decoder.setCodeword(Codewordtemp);
	//cout << "true codeword checksum:" << Decoder.check() <<endl;

	Decoder.ReadH();
	while(pckerror < 100)
	{
		blkerror = 0;

		//for(i = 247; i >= 0; i--)
		//	cout << std::bitset<8>(InfoStream[i]) << ", ";
		//cout <<endl;
		Encoder.encode(InfoStream, 122);//WIFI code
		
		//cout << "codeword checksum:" << Decoder.check_fp(Codewordtemp) <<endl;
		//cout << "codeword checksum: (0 is pass, 1 is not pass):" << Decoder.check_fp(Codeword)<<endl;
		for(i = 0; i < CWD_LENGTH; i++)
		{
			//-- BPSK modulation => 0 -> 1, 1 -> -1
			// Using Box-Muller method to generate Gaussian noise
			LLR[i] = 2*snr*(1 - 2*Encoder.getCodeword(i) + Normal(0, sigma));
			// -- all zero
			//LLR[i] = 2*snr*(1 + Normal(0, sigma));
			//-- noiseless case
			//LLR[i] = 2*snr*(1 - 2*Encoder.getCodeword(i));
			// Using Wallace method to generate Gaussian noise
			//LLR[i] = 2*EbN0*(1 + Wallace(0, sigma));
			//LLR_fp[i] = int(LLR[i]*(1<<FRAC_WIDTH));
		}
		Decoder.setState(PCV);
		
		
		//cout << "checking hard output before decoding:" << Decoder.hardDecision(LLR_fp) << endl;
		//blkerror = Decoder.decode(LLR);
		//cout << Decoder.decode_fixpoint(LLR_fp) << " total iterations" <<endl;
		Decoder.decode(LLR);
		Decoder.resetBER();
		blkerror = Decoder.calculateBER();
		if(blkerror > 0)
			pckerror++;
		biterror += blkerror;
		Counter++;
	}
	cout << biterror << " " << pckerror << " " << Counter << endl
		<< " FER: " << pckerror/Counter << " BER: " << biterror/Counter/CWD_LENGTH << endl;
	//Decoder.getInfoBit();
	return 0;
}

//---------------------------- Non in-class code
void noMoreMemory()
{
	cerr << "Unable to satisfy request for memory\n";
	abort();
}
int DecodeTrial(double EbN0_dB, int MaxPacket)
{
	double LLR[CWD_LENGTH];
	int LLR_fp[100][CWD_LENGTH];
	int i, j; 
	double snr, sigma, Rate = 0;; 
	class FP_Decoder Decoder;
	__int64 ctr1 = 0, ctr2 = 0, freq = 0;

	//EbN0_dB = 4.5;
	Rate = Decoder.getRate();
	snr = 2*pow(10.0, EbN0_dB/10)*Rate;			// Eb/N0	to snr
	sigma = sqrt(1/snr);
   cout << "equivalent SNR is: " << 10*log10(snr) << endl;
	for(j = 0; j < 100; j++)
	{
		for(i = 0; i < CWD_LENGTH; i++)
		{
			//-- BPSK modulation => 0 -> 1, 1 -> -1 assume all zero codeword
			// Using Box-Muller method to generate Gaussian noise
			LLR[i] = 2*snr*(1 + Normal(0, sigma));
			LLR_fp[j][i] = int(LLR[i]*(1<<FRAC_WIDTH));
		}
	}
	// Start timing the code.for(i = 0; i < CWD_LENGTH; i++)
   if(QueryPerformanceCounter((LARGE_INTEGER *) &ctr1) != 0)
	{
		// Do what ever you do, what ever you need to time...
		
		// Finish timing the code.
		for(i = 0; i < MaxPacket; i++)
		{
			Decoder.setState(PCV);
			Decoder.decode_fixpoint(LLR_fp[i%100]);
		}
      QueryPerformanceCounter((LARGE_INTEGER *) &ctr2);
      QueryPerformanceFrequency((LARGE_INTEGER *) &freq);

      // Print the time spent in microseconds to the console.

      std::cout << ((ctr2 - ctr1) * 1.0 / freq) << "  seconds" << std::endl;
		std::cout << 2209*MaxPacket/((ctr2 - ctr1) * 1.0 / freq) << " bits per second for decoder" << std::endl;
   }
	return 0;
}
int EncodeTrial(char *info, int MaxPacket)
{
	class FP_Encoder Encoder("G_array_forward.txt", 0);
	int i;

	__int64 ctr1 = 0, ctr2 = 0, freq = 0;
      // Start timing the code.
   if(QueryPerformanceCounter((LARGE_INTEGER *) &ctr1) != 0)
	{
		for(i = 0; i < MaxPacket; i++)
		{
			Encoder.encode(info, 248);
		}
		QueryPerformanceCounter((LARGE_INTEGER *) &ctr2);
      QueryPerformanceFrequency((LARGE_INTEGER *) &freq);

         // Print the time spent in microseconds to the console.

      cout << ((ctr2 - ctr1) * 1.0 / freq) << "  seconds" << endl;
		cout << 2209*MaxPacket/((ctr2 - ctr1) * 1.0 / freq) << " bits per second for encoder" << endl;
	}
	return 0;
}

int ArrayLDPC_Debug()
{
	//2209 - 231 = 1978
	//1978 / 8 = 247, 1978 % 8 = 2
	char InfoStream[248] = "OMG how long should this string be to make it 248, just imagine that. \
								   I guess it's still not long enough. Let's see. This is a testing string \
									for a lot of characters so that we have some random bit stream that's\
									correct";
	//--- unused debugging varibles  
	//char InfoStream[248] = "";
	//int InfoBit[INFO_LENGTH];
	//int Codewordtemp[CWD_LENGTH];
	//double Receive[CWD_LENGTH];
	//--- unused debugging varibles
	int i, j;
	int LLR_fp[CWD_LENGTH];
	//int pckerror = 0;
	double biterror = 0;
	long Counter = 0;
	long Clk = 0;

	double LLR[CWD_LENGTH];
	double db = 1.0;
	double Rate = 0;
	double EbN0_dB = 4;
	double snr = 0;
	double sigma;
	double blkerror = 0; 
	double pckerror = 0;
	int info_indx[INFO_LENGTH];
	
	class FP_Decoder Decoder;
	class FP_Encoder Encoder("G_array_forward.txt", 0);
	
	EbN0_dB = 4.5;
	Rate = Decoder.getRate();
	snr = 2*pow(10.0, EbN0_dB/10)*Rate;			// Eb/N0	to snr
	sigma = sqrt(1/snr);
	cout << "SNR is " << 10*log10(snr) << " dB" << endl;
	biterror = 0;
	Clk = 0;

	//---- debugging
	char out[277] = "";
	//out = new char[553];
	// working now
	Encoder.encode(InfoStream, out, 248);
	//---- 

	Decoder.setInfoBit(InfoStream, 248);
	for(i = 0; i < INFO_LENGTH; i++)
	{
		info_indx[i] = Encoder.getInfoIndex(i);
	}
	Decoder.setInfoIndex(info_indx);
	//for(i = 0; i < CWD_LENGTH; i++)
	//		Codewordtemp[i] = Encoder.getCodeword(i);
	//Decoder.setCodeword(Codewordtemp);
	//cout << "true codeword checksum:" << Decoder.check() <<endl;
	while(pckerror < 100)
	{
		blkerror = 0;

		//for(i = 247; i >= 0; i--)
		//	cout << std::bitset<8>(InfoStream[i]) << ", ";
		//cout <<endl;
		Encoder.encode(InfoStream, 248);
		
		//cout << "codeword checksum:" << Decoder.check_fp(Codewordtemp) <<endl;
		//cout << "codeword checksum: (0 is pass, 1 is not pass):" << Decoder.check_fp(Codeword)<<endl;
		for(i = 0; i < CWD_LENGTH; i++)
		{
			//-- BPSK modulation => 0 -> 1, 1 -> -1
			// Using Box-Muller method to generate Gaussian noise
			LLR[i] = 2*snr*(1 - 2*Encoder.getCodeword(i) + Normal(0, sigma));
			//-- noiseless case
			//LLR[i] = 2*snr*(1 - 2*Encoder.getCodeword(i));
			// Using Wallace method to generate Gaussian noise
			//LLR[i] = 2*EbN0*(1 + Wallace(0, sigma));
			LLR_fp[i] = int(LLR[i]*(1<<FRAC_WIDTH));
		}
		Decoder.setState(PCV);
		
		
		//cout << "checking hard output before decoding:" << Decoder.hardDecision(LLR_fp) << endl;
		//blkerror = Decoder.decode(LLR);
		//cout << Decoder.decode_fixpoint(LLR_fp) << " total iterations" <<endl;
		Decoder.decode_fixpoint(LLR_fp);
		Decoder.resetBER();
		blkerror = Decoder.calculateBER();
		if(blkerror > 0)
			pckerror++;
		biterror += blkerror;
		Counter++;
	}
	cout << biterror << " " << pckerror << " " << Counter << endl
		<< " FER: " << pckerror/Counter << " BER: " << biterror/Counter/CWD_LENGTH << endl;
	//Decoder.getInfoBit();
	return 0;
}

int ArrayLDPC_Debug_Shorten(int short_len)
{
	//2209 - 231 = 1978
	//1978 / 8 = 247, 1978 % 8 = 2
	char InfoStream[248] = "OMG how long should this string be to make it 248, just imagine that. \
								   I guess it's still not long enough. Let's see. This is a testing string \
									for a lot of characters so that we have some random bit stream that's\
									correct";
	//--- unused debugging varibles  
	//char InfoStream[248] = "";
	//int InfoBit[INFO_LENGTH];
	//int Codewordtemp[CWD_LENGTH];
	//double Receive[CWD_LENGTH];
	//--- unused debugging varibles
	int i, j;
	int LLR_fp[CWD_LENGTH];
	//int pckerror = 0;
	double biterror = 0;
	long Counter = 0;
	long Clk = 0;

	double LLR[CWD_LENGTH];
	double db = 1.0;
	double Rate = 0;
	double EbN0_dB = 4;
	double snr = 0;
	double sigma;
	double blkerror = 0; 
	double pckerror = 0;
	int info_indx[INFO_LENGTH];
	
	class FP_Decoder Decoder;
	class FP_Encoder Encoder("G_array_forward.txt", 0);
	
	EbN0_dB = 4.5;
	Rate = Decoder.getRate();
	snr = 2*pow(10.0, EbN0_dB/10)*Rate;			// Eb/N0	to snr
	snr = 2*pow(10.0, EbN0_dB/10)*(1978.0 - 976.0)/2209.0;
	sigma = sqrt(1/snr);
	cout << "SNR is " << 10*log10(snr) << " dB" << endl;
	biterror = 0;
	Clk = 0;

	//---- shortening, try 976 out of 1978 since it corresondes to 122 char
	for(i = 0; i < short_len; i ++)
	{
		//--- first short_len number of string becomes 
		InfoStream[i] = 0;
	}

	//---- debugging
	//char out[277] = "";
	//out = new char[553];
	// working now
	//Encoder.encode(InfoStream, out, 248);
	//---- 
	//Decoder.setShorten(short_len);
	Decoder.setInfoBit(InfoStream, 248);
	for(i = 0; i < INFO_LENGTH; i++)
	{
		info_indx[i] = Encoder.getInfoIndex(i);
	}
	Decoder.setInfoIndex(info_indx);
	//for(i = 0; i < CWD_LENGTH; i++)
	//		Codewordtemp[i] = Encoder.getCodeword(i);
	//Decoder.setCodeword(Codewordtemp);
	//cout << "true codeword checksum:" << Decoder.check() <<endl;
	while(pckerror < 100)
	{
		blkerror = 0;

		//for(i = 247; i >= 0; i--)
		//	cout << std::bitset<8>(InfoStream[i]) << ", ";
		//cout <<endl;
		Encoder.encode(InfoStream, 248);
		
		//cout << "codeword checksum:" << Decoder.check_fp(Codewordtemp) <<endl;
		//cout << "codeword checksum: (0 is pass, 1 is not pass):" << Decoder.check_fp(Codeword)<<endl;
		for(i = 0; i < CWD_LENGTH; i++)
		{
			//-- BPSK modulation => 0 -> 1, 1 -> -1
			// Using Box-Muller method to generate Gaussian noise
			LLR[i] = 2*snr*(1 - 2*Encoder.getCodeword(i) + Normal(0, sigma));
			//-- noiseless case
			//LLR[i] = 2*snr*(1 - 2*Encoder.getCodeword(i));
			// Using Wallace method to generate Gaussian noise
			//LLR[i] = 2*EbN0*(1 + Wallace(0, sigma));
			LLR_fp[i] = int(LLR[i]*(1<<FRAC_WIDTH));
		}
		Decoder.setState(PCV);
		
		//----- shortening part
		for(i = 0; i < short_len; i ++)
		{
			// set to max integer value
			LLR_fp[info_indx[i]] = 7*(1<<FRAC_WIDTH );
		}
		
		//cout << "checking hard output before decoding:" << Decoder.hardDecision(LLR_fp) << endl;
		//blkerror = Decoder.decode(LLR);
		//cout << Decoder.decode_fixpoint(LLR_fp) << " total iterations" <<endl;
		cout << Decoder.decode_fixpoint(LLR_fp) << ", ";
		Decoder.resetBER();
		blkerror = Decoder.calculateBER();
		if(blkerror > 0)
			pckerror++;
		biterror += blkerror;
		Counter++;
	}
	cout << biterror << " " << pckerror << " " << Counter << endl
		<< " FER: " << pckerror/Counter << " BER: " << biterror/Counter/CWD_LENGTH << endl;
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
	//int pckerror = 0;
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
	double sigma;
	double blkerror = 0; 
	double pckerror = 0;
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
	sigma = sqrt(1/EbN0);
	BitError = 0;
	Clk = 0;
	while(pckerror < 100)
	{
		blkerror = 0;
		for(i = 0; i < CWD_LENGTH; i++)
		{
			//-- All zero codeword
			// Using Box-Muller method to generate Gaussian noise
			LLR[i] = 2*EbN0*(1 + Normal(0, sigma));
			// Using Wallace method to generate Gaussian noise
			//LLR[i] = 2*EbN0*(1 + Wallace(0, sigma));
			LLR_fp[i] = int(LLR[i]*(1<<FRAC_WIDTH));
			//Debug
			//LLR_fp[i] = i % 47;
		}
		Decoder.setState(PCV);
		//blkerror = Decoder.decode(LLR);
		blkerror = Decoder.decode_fixpoint(LLR_fp);
		if(blkerror > 0)
			pckerror++;
		BitError += blkerror;
		Counter++;
	}
	cout << BitError << " " << pckerror << " " << Counter << endl
		<< " FER: " << pckerror/Counter << " BER: " << BitError/Counter/CWD_LENGTH << endl;
	//Decoder.getInfoBit();
	return 0;
}


int ArrayLDPC_TimeTrial(double db, int MaxPckNum, char* Filename)
{
	char Filename2[40]="";
	char InfoStream[248] = "OMG how long should this string be to make it 248, just imagine that. \
								   I guess it's still not long enough. Let's see. This is a testing string \
									for a lot of characters so that we have some random bit stream that's\
									correct";
	int i, j;
	int CurrentInd = 0;
	int LLR_fp[CWD_LENGTH];
	double BitError = 0;
	long Counter = 0;
	long Clk = 0;

	//double Receive[CWD_LENGTH];
	double LLR[CWD_LENGTH];
	double snr;
	//double db = 1.0;
	double Rate = 0;
	double EbN0_dB = 4;
	double EbN0 = 0;
	double sigma;
	double blkerror = 0; 
	double pckerror = 0;
	//double **MemBank; //Memory bank
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
	sigma = sqrt(1/EbN0);
	BitError = 0;
	Clk = 0;
	while(Counter < MaxPckNum)
	{
		blkerror = 0;
		for(i = 0; i < CWD_LENGTH; i++)
		{
			//-- All zero codeword
			// Using Box-Muller method to generate Gaussian noise
			LLR[i] = 2*EbN0*(1 + Normal(0, sigma));
			// Using Wallace method to generate Gaussian noise
			//LLR[i] = 2*EbN0*(1 + Wallace(0, sigma));
			LLR_fp[i] = int(LLR[i]*(1<<FRAC_WIDTH));
			//Debug
			//LLR_fp[i] = i % 47;
		}
		Decoder.setState(PCV);
		//blkerror = Decoder.decode(LLR);
		blkerror = Decoder.decode_fixpoint(LLR_fp);
		if(blkerror > 0)
			pckerror++;
		BitError += blkerror;
		Counter++;
	}
	cout << BitError << " " << pckerror << " " << Counter << endl
		<< " FER: " << pckerror/Counter << " BER: " << BitError/Counter/CWD_LENGTH << endl;
	//Decoder.getInfoBit();
	
	return 0;
}


