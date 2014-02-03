#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rngs.h"
#include "rvgs.h"
#include "MemoryMacro.h"
#include "Memory.h"
#include "ArrayLDPC.h"
#include "ArrayLDPCMacro.h"
using namespace std;

int FP_Decoder::sxor(int x, int y){		
	int v1, v2;
	int sum;
	int diff;
	int part1;
	int part2;
	int i;
	v1 = abs(x);
	v2 = abs(y);
	sum = (v1+v2) & WIDTH_MASK;
	diff = abs(v1-v2) & WIDTH_MASK;
	part1 = Constant - (sum >> 2); //divided by 4
	part1 = (part1 > 0) ? part1:0;
	part2 = Constant - (diff >> 2); 
	part2 = (part2 > 0) ? part2:0;
	//cout << sgn(x)*sgn(y)*(min(v1, v2) + part1 - part2) << endl;
	return sgn(x)*sgn(y)*(min(v1, v2) + part1 - part2);
}

int FP_Decoder::decode_fixpoint(const int *LLR)
{
	int Counter, Iteration; 
	int i, j, k;
	//int VarShift[NUM_VGRP];
	static int VarAddr;
	static int ChkAddr;
	static int Shift;
	static int Accum[CIR_SIZE];
	static int BankSelect;
	static int Reg_v2c_pre[CIR_SIZE][CHK_DEG];
	static int Reg_c2v[CIR_SIZE][CHK_DEG];
	static int Reg_c2v_pre[CIR_SIZE][VAR_DEG];
	static int Reg_v2c[CHK_DEG];
	static int Forward[CIR_SIZE][CHK_DEG];
	static int Backward[CIR_SIZE][CHK_DEG];
	//Initialize the Edge RAM with channel values from variable node 
	int temp;
	
	if(FSM.getState() == PCV) // If prepare channel value state is on
	{
		for(i = 0; i < NUM_VGRP; i++)
		{
			for(k = 0; k < VAR_DEG; k++)
			{
				for(j = 0; j < CIR_SIZE; j++)	//-- Parallel Process
				{
					VarAddr = i*CIR_SIZE + j;
					Shift = CodeROM.getCirShift(k, i);
					BankSelect = (j + Shift) % CIR_SIZE;	//Select which bank of RAM is active
					ChkAddr = (i + k*CHK_DEG);	//The Address within the selected bank of RAM
					EdgeRAM[BankSelect].setAddress(ChkAddr);	
					EdgeRAM[BankSelect].wrtData(LLR[VarAddr]);
				}
			}
		}
		FSM.setState(C2V);	// go to check to variable state
	}
	Iteration = 0;
	while(Iteration < MAX_ITER && FSM.getState() == C2V)
	{
		//-- Process each check group at a time. 
		for(i = 0; i < NUM_CGRP; i++)
		{
			//-- CIR_SIZE banks of RAM and serially retrieve CHK_DEG messages from the memory
			for(k = 0; k < CHK_DEG; k++)
			{
				//-- Check node parallel process denoted in for loop (loop through all check nodes in one group)
				for(j = 0; j < CIR_SIZE; j++)		
				{
					//-- Prepare all the V2C message
					VarAddr = k + i*CHK_DEG;
					EdgeRAM[j].setAddress(VarAddr); 
					Reg_v2c_pre[j][k] = EdgeRAM[j].rdData();	
				}
			}
			//-- Do binary operation of the sxor: forward backward computation
			for(j = 0; j < CIR_SIZE; j++)	//-- Parallel process
			{
				Forward[j][0] = Reg_v2c_pre[j][0];
				Backward[j][CHK_DEG-1] = Reg_v2c_pre[j][CHK_DEG-1];
			}
			for(k = 1; k < CHK_DEG-1; k++) //-- Serial process
			{
				for(j = 0; j < CIR_SIZE; j++)	//-- Parallel process
				{
					Forward[j][k] = sxor(Forward[j][k-1], Reg_v2c_pre[j][k]);
					Backward[j][CHK_DEG - k - 1] = sxor(Backward[j][CHK_DEG - k], Reg_v2c_pre[j][CHK_DEG - 1 - k]);
				}
			}
			////-- Compute MC2V and accumulate the posteriori:
			//for(j = 0; j < CIR_SIZE; j++)//-- Parallel process, initialize the accumulate register
			//{
			//	Accum[j] = 0;
			//}
			//-- 1. Compute the first and last MC2V and write back to RAM
			for(j = 0; j < CIR_SIZE; j++)	//-- Parallel process
			{
				k = 0; // 1st
				Reg_c2v[j][k] = Backward[j][k+1];
				VarAddr = k + i*CHK_DEG;
				EdgeRAM[j].setAddress(VarAddr); 
				EdgeRAM[j].wrtData(Reg_c2v[j][k]);

				//Reg_c2v[CHK_DEG-1] = Forward[CHK_DEG-1];
				k = CHK_DEG-1; // last
				Reg_c2v[j][k] = Forward[j][k-1];
				VarAddr = k + i*CHK_DEG;
				EdgeRAM[j].setAddress(VarAddr); 
				EdgeRAM[j].wrtData(Reg_c2v[j][k]);
			}
			//-- 2. Compute rest of the MC2V message and write back to the RAM
			for(k = 1; k < CHK_DEG-1; k++)
			{							
				for(j = 0; j < CIR_SIZE; j++) //-- Parallel process
				{
					Reg_c2v[j][k] = sxor(Forward[j][k-1], Backward[j][k+1]);
					VarAddr = k + i*CHK_DEG;
					EdgeRAM[j].setAddress(VarAddr); 
					EdgeRAM[j].wrtData(Reg_c2v[j][k]);
				}
			}
		}
		FSM.setState(V2C);
		//-- Check node group process complete

		//-- Prepare MC2V and compute the posteriori and new MV2C
		//if(FSM.getState() == V2C)	
		for(i = 0; i < NUM_VGRP; i++)
		{
			for(j = 0; j < CIR_SIZE; j++)//-- Parallel process
			{
				Accum[j] = 0;
			}
			for(k = 0; k < VAR_DEG; k++)	//-- Serially accumulate
			{
				// Accumulate all the MC2V message
				for(j = 0; j < CIR_SIZE; j++)	//-- Parallel process
				{
					VarAddr = i*CIR_SIZE + j;
					Shift = CodeROM.getCirShift(k, i);
					BankSelect = (j + Shift) % CIR_SIZE;	//Select which bank of RAM is active
					ChkAddr = (i + k*CHK_DEG);	//The Address within the selected bank of RAM
					EdgeRAM[BankSelect].setAddress(ChkAddr);	
					Reg_v2c_pre[j][k] = EdgeRAM[BankSelect].rdData();
					Accum[j] = Accum[j] + Reg_v2c_pre[j][k];
				}
			}
			for(j = 0; j < CIR_SIZE; j++)	//-- Parallel process
			{
				VarAddr = i*CIR_SIZE + j;
				// Add in the channel value (fixed point)
				Accum[j] = Accum[j] + LLR[VarAddr];
				// Update the Posteriori
				wrtPost(VarAddr, Accum[j]);
			}
				
			for(k = 0; k < VAR_DEG; k++) //-- Write back to RAM the MV2C (serial)
			{
				for(j = 0; j < CIR_SIZE; j++) //-- Parallel process
				{
					Shift = CodeROM.getCirShift(k, i);
					BankSelect = (j + Shift) % CIR_SIZE;
					ChkAddr = (i + k*CHK_DEG);
					EdgeRAM[BankSelect].setAddress(ChkAddr);	
					EdgeRAM[BankSelect].wrtData(Accum[j] - Reg_v2c_pre[j][k]);
				}
			}
		}
		Iteration++;
		// If there is no error then break (cannot detect codeword error) (fixed point)
		if(!checkPost_fp())
		{
			return 0;
		}
		else
		{
			FSM.setState(C2V);
		}
	}// While Iteration loop end
	checkPost_fp();
	return BitError;
}

inline int FP_Decoder::fmax(int x, int y){		
	if(x >= y)
		return x;
	else
		return y;
}
double FP_Decoder::sxor(double x, double y){		
	double v1, v2;
	double sum_abs, diff_abs;
	v1 = fabs(x);
	v2 = fabs(y);
	sum_abs = v1+v2;
	diff_abs = fabs(v1-v2);
	return sgn(x)*sgn(y)*(min(v1, v2) + log(1 + exp(-sum_abs)) - log(1 + exp(-diff_abs)) );
}





int FP_Decoder::decode(const double *LLR)
{
	//int Counter; 
	int i, j, k;
	//int VarShift[NUM_VGRP];
	int VarAddr;
	int ChkAddr;
	int Shift;
	int Iteration = 0;
	double MV2C[CHK_DEG];
	double MC2V[VAR_DEG];
	double Forward[CHK_DEG];
	double Backward[CHK_DEG];
	//Initialize the Edge RAM with channel values from variable node 
	double temp;
	for(i = 0; i < NUM_VGRP; i++)
	{
		for(j = 0; j < CIR_SIZE; j++)
		{
			VarAddr = i*CIR_SIZE + j;
			for(k = 0; k < VAR_DEG; k++)
			{
				Shift = CodeROM.getCirShift(k, i);
				ChkAddr = (j + Shift) % CIR_SIZE;
				EdgeRAM[k].setAddress(ChkAddr + i*CIR_SIZE);	
				EdgeRAM[k].wrtData(LLR[VarAddr]);
				//temp = EdgeRAM[k].getData();
			}
		}
	}
	while(Iteration < MAX_ITER)
	{
		//-- Process each check group at a time. 
		for(i = 0; i < NUM_CGRP; i++)
		{
			//-- Check node parallel process denoted in for loop (loop through all check nodes in one group)
			for(j = 0; j < CIR_SIZE; j++)
			{
				
				for(k = 0; k < CHK_DEG; k++)
				{
					//This is simply for better read of the code. Not necessary.
					Shift = CodeROM.getCirShift(i, k);
					//-- Prepare all the V2C message
					//VarAddr = k*47 + Shift;
					VarAddr = j + k*CHK_DEG;
					EdgeRAM[i].setAddress(VarAddr); 
					MV2C[k] = EdgeRAM[i].getData();	
				}
				//-- Do binary operation of the sxor: forward backward computation
				Forward[0] = MV2C[0];
				Backward[CHK_DEG-1] = MV2C[CHK_DEG-1];
				for(k = 1; k < CHK_DEG-1; k++)
				{
					Forward[k] = sxor(Forward[k-1], MV2C[k]);
					Backward[CHK_DEG - k - 1] = sxor(Backward[CHK_DEG - k], MV2C[CHK_DEG - 1 - k]);
				}
				//-- make MV2C temp memory to store all the computed MC2V message from one check node
				//-- Compute the first and last MC2V and write back to RAM
				k = 0; // 1st
				MV2C[k] = Backward[k+1];
				VarAddr = j + k*CHK_DEG;
				EdgeRAM[i].setAddress(VarAddr); 
				EdgeRAM[i].wrtData(MV2C[k]);

				//MV2C[CHK_DEG-1] = Forward[CHK_DEG-1];
				k = CHK_DEG-1; // last
				MV2C[k] = Forward[k-1];
				VarAddr = j + k*CHK_DEG;
				EdgeRAM[i].setAddress(VarAddr); 
				EdgeRAM[i].wrtData(MV2C[k]);
				
				//-- Compute the MV2C message and write back to the RAM
				for(k = 1; k < CHK_DEG-1; k++)
				{
					MV2C[k] = sxor(Forward[k-1], Backward[k+1]);
					VarAddr = j + k*CHK_DEG;
					EdgeRAM[i].setAddress(VarAddr); 
					EdgeRAM[i].wrtData(MV2C[k]);
				}

			}
		}
		//-- Check node group process complete

		//-- Prepare MC2V and compute the posteriori and new MV2C
		for(i = 0; i < NUM_VGRP; i++)
		{
			for(j = 0; j < CIR_SIZE; j++)
			{
				VarAddr = i*CIR_SIZE + j;
				temp = 0;
				// Accumulate all the MC2V message
				for(k = 0; k < VAR_DEG; k++)
				{
					Shift = CodeROM.getCirShift(k, i);
					ChkAddr = (j + Shift) % CIR_SIZE;
					EdgeRAM[k].setAddress(ChkAddr + i*CIR_SIZE);	
					MC2V[k] = EdgeRAM[k].getData();
					temp = temp + MC2V[k];
					//temp = EdgeRAM[k].getData();
				}
				// Add in the channel value
				temp = temp + LLR[VarAddr];
				wrtPost(VarAddr, temp);
				for(k = 0; k < VAR_DEG; k++)
				{
					Shift = CodeROM.getCirShift(k, i);
					ChkAddr = (j + Shift) % CIR_SIZE;
					EdgeRAM[k].setAddress(ChkAddr + i*CIR_SIZE);	
					EdgeRAM[k].wrtData(temp - MC2V[k]);
					//temp = EdgeRAM[k].getData();
				}
			}
		}
		Iteration++;
		// If there is no error then break (cannot detect codeword error)
		// early termination part
		//if(!checkPost())
		//{
		//	return 0;
		//}
	}
	checkPost();
	return BitError;
}

int FP_Decoder::checkPost()
{
	int i = 0;
	BitError = 0;
	double temp = 0;
	for(i = 0; i < CWD_LENGTH; i++)
	{
		temp = getPost(i);
		if(getPost(i) < 0)
			BitError++;
	}
	if(BitError != 0)
		return 1;
	else
		return 0;
}

//---------------------------- old stuff=> moved to PerfTest.cpp
//void noMoreMemory()
//{
//	cerr << "Unable to satisfy request for memory\n";
//	abort();
//}


//int ArrayLDPC_PerfTest(double db_start, double db_end, double db_step, char* Filename)
//{
//	char Filename2[40]="";
//	int InfoBit[INFO_LENGTH];
//	int i, j;
//	int CurrentInd = 0;
//	int LLR_fp[CWD_LENGTH];
//	//int PckError = 0;
//	double BitError = 0;
//	long Counter = 0;
//	long Clk = 0;
//
//	double Receive[CWD_LENGTH];
//	double LLR[CWD_LENGTH];
//	double snr;
//	double db = 1.0;
//	double Rate = 0;
//	double EbN0_dB = 4;
//	double EbN0 = 0;
//	double Sigma;
//	double BlkError = 0; 
//	double PckError = 0;
//	double **MemBank; //Memory bank
//	bool per_flag;
//
//	class FP_Decoder Decoder;
//	class ControlFSM FSM;
//
//	ofstream FilePtr(Filename);
//	if(!FilePtr)
//	{
//		cerr << "failed to open " << Filename << endl;
//		system("pause");
//		exit(0);
//	}
//	FilePtr.close();
//	//strcat(Filename2, "log_");
//	//strcat(Filename2, Filename);
//	//strcat(Filename2, ".txt");
//	sprintf(Filename2, "%s_log.txt", Filename);
//	ofstream LogPtr(Filename2);
//	if(!LogPtr)
//	{
//		cerr << "failed to open " << Filename2 << endl;
//		system("pause");
//		exit(0);
//	}
//	LogPtr.close();
//
//	
//	//cout << "SNR? ";
//	//cin >> EbN0_dB;
//	EbN0_dB = db_start;
//	Rate = Decoder.getRate();
//	EbN0 = 2*pow(10.0, EbN0_dB/10)*Rate;			// Eb/N0	
//	Sigma = sqrt(1/EbN0);
//	BitError = 0;
//	Clk = 0;
//	while(PckError < 100)
//	{
//		BlkError = 0;
//		for(i = 0; i < CWD_LENGTH; i++)
//		{
//			//-- All zero codeword
//			// Using Box-Muller method to generate Gaussian noise
//			LLR[i] = 2*EbN0*(1 + Normal(0, Sigma));
//			// Using Wallace method to generate Gaussian noise
//			//LLR[i] = 2*EbN0*(1 + Wallace(0, Sigma));
//			LLR_fp[i] = int(LLR[i]*(1<<FRAC_WIDTH));
//			//Debug
//			//LLR_fp[i] = i % 47;
//		}
//		Decoder.setState(PCV);
//		//BlkError = Decoder.decode(LLR);
//		BlkError = Decoder.decode_fixpoint(LLR_fp);
//		if(BlkError > 0)
//			PckError++;
//		BitError += BlkError;
//		Counter++;
//	}
//	cout << BitError << " " << PckError << " " << Counter << endl
//		<< " FER: " << PckError/Counter << " BER: " << BitError/Counter/CWD_LENGTH << endl;
//	//Decoder.getInfoBit();
//	return 0;
//}
//
//int ArrayLDPC_TimeTrial(double db, int MaxPckNum, char* Filename)
//{
//	char Filename2[40]="";
//	int InfoBit[INFO_LENGTH];
//	int i, j;
//	int CurrentInd = 0;
//	int LLR_fp[CWD_LENGTH];
//	//int PckError = 0;
//	double BitError = 0;
//	long Counter = 0;
//	long Clk = 0;
//
//	double Receive[CWD_LENGTH];
//	double LLR[CWD_LENGTH];
//	double snr;
//	//double db = 1.0;
//	double Rate = 0;
//	double EbN0_dB = 4;
//	double EbN0 = 0;
//	double Sigma;
//	double BlkError = 0; 
//	double PckError = 0;
//	double **MemBank; //Memory bank
//	bool per_flag;
//
//	class FP_Decoder Decoder;
//	class ControlFSM FSM;
//
//	ofstream FilePtr(Filename);
//	if(!FilePtr)
//	{
//		cerr << "failed to open " << Filename << endl;
//		system("pause");
//		exit(0);
//	}
//	FilePtr.close();
//	//strcat(Filename2, "log_");
//	//strcat(Filename2, Filename);
//	//strcat(Filename2, ".txt");
//	sprintf(Filename2, "%s_log.txt", Filename);
//	ofstream LogPtr(Filename2);
//	if(!LogPtr)
//	{
//		cerr << "failed to open " << Filename2 << endl;
//		system("pause");
//		exit(0);
//	}
//	LogPtr.close();
//
//	
//	//cout << "SNR? ";
//	//cin >> EbN0_dB;
//	EbN0_dB = db;
//	Rate = Decoder.getRate();
//	EbN0 = 2*pow(10.0, EbN0_dB/10)*Rate;			// Eb/N0	
//	Sigma = sqrt(1/EbN0);
//	BitError = 0;
//	Clk = 0;
//	while(Counter < MaxPckNum)
//	{
//		BlkError = 0;
//		for(i = 0; i < CWD_LENGTH; i++)
//		{
//			//-- All zero codeword
//			// Using Box-Muller method to generate Gaussian noise
//			LLR[i] = 2*EbN0*(1 + Normal(0, Sigma));
//			// Using Wallace method to generate Gaussian noise
//			//LLR[i] = 2*EbN0*(1 + Wallace(0, Sigma));
//			LLR_fp[i] = int(LLR[i]*(1<<FRAC_WIDTH));
//			//Debug
//			//LLR_fp[i] = i % 47;
//		}
//		Decoder.setState(PCV);
//		//BlkError = Decoder.decode(LLR);
//		BlkError = Decoder.decode_fixpoint(LLR_fp);
//		if(BlkError > 0)
//			PckError++;
//		BitError += BlkError;
//		Counter++;
//	}
//	cout << BitError << " " << PckError << " " << Counter << endl
//		<< " FER: " << PckError/Counter << " BER: " << BitError/Counter/CWD_LENGTH << endl;
//	//Decoder.getInfoBit();
//	return 0;
//}
//
//
