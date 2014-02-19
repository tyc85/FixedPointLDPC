#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rngs.h"
#include "rvgs.h"
#include "Memory.h"
//#include "ArrayLDPC.h"
#include "ArrayLDPCMacro.h"
using namespace std;
//--------- notes -----------
// 1. don't really need the posterior stored. just the hard decision would be enough
// 
//---------------------------
//---- new part
// new version that outputs character stream
int FP_Decoder::decode(const int *LLR, unsigned char out*)
{
	int i, j;
	// !!!! Assuming out is all zero?
	decode_fixpoint(LLR);
	for(i = 0; i < 276; i++)
	{
		
		// NOTE: this should solve the problem where out might not be null
		// since DecodedCodeword is either 0 or 1, erasing all other bits
		j = 0;
		out[i] = DecodedCodeword[i*8 + j];
		for(j = 1; j < 8; j++)
		{
			out[i] = out[i] + (DecodedCodeword[i*8 + j] << j);
		}
	}
	i = 277;
	j = 0;
	out[i] = DecodedCodeword[i*8 + j];
	for(j = 1; j < 2209 % 8; j++)
	{
		out[i] = out[i] + (DecodedCodeword[i*8 + j] << j);
	}
}


//-------------- Decoder part starts
int FP_Decoder::decode_fixpoint(const int *LLR)
{
	int Iteration; 
	int i, j, k;
	//---- 
	int var_add, chk_add, cir_shift;
	int temp;
	//---- 
	static int VarAddr;
	static int ChkAddr;
	static int Shift;
	static int Accum[CIR_SIZE];
	static int BankSelect;
	static int Reg_v2c_pre[CIR_SIZE][CHK_DEG];
	static int Reg_c2v[CIR_SIZE][CHK_DEG];
	static int Reg_c2v_pre[CIR_SIZE][VAR_DEG]; // actually don't need pre?
	//static int Reg_v2c[CHK_DEG];
	static int Forward[CIR_SIZE][CHK_DEG];
	static int Backward[CIR_SIZE][CHK_DEG];


	if(!hardDecision(LLR)) 
	{
		//resetBER();
		//calculateBER();
		//cout << "check sum pass before decoding (continue to decode anyway)" << endl;
		//return BitError;
		return 0;
	}
	//else
	//{
	//	resetBER();
	//	calculateBER();
	//	cout << BitError << " bit errors before decoding start" << endl;
	//}
	
	/*----
	The BRCM is check oriented. The memory address specity which check it is, 
	E.g. EdgeRAM[i].setAddress(j) specity the ith edge of the jth check node. 
	-----*/
	if(FSM.getState() == PCV) // prepare channel value state is on
	{
		// this is a totally rewritted code
		// writing channel value to each edge that connects with the each vnode
		for(i = 0; i < NUM_CGRP; i++) // go through all groups of cnode
		{
			for(k = 0; k < CIR_SIZE; k++) // for each cnode in each cgroup
			{
				//-- Parallel Process: writing one elem to each bank
				for(j = 0; j < CHK_DEG; j++)	
				{
					//-- debug
					Shift = CodeROM.getCirShift(i, j);
					ChkAddr = k + i*CIR_SIZE;
					VarAddr = (k + Shift) % CIR_SIZE + j*CIR_SIZE;// if bug verify this address
					BankSelect = j; //simply which vgroup we are at
					EdgeRAM[BankSelect].setAddress(ChkAddr);	
					EdgeRAM[BankSelect].wrtData(LLR[VarAddr]);
					//-- debug end
				}
			}
		}
		FSM.setState(C2V);	// go to check to variable state
	}

	Iteration = 0;
	while(Iteration < MAX_ITER && FSM.getState() == C2V)
	{
		//-- Process one check group at a time. 
		for(i = 0; i < NUM_CGRP; i++)
		{
			//-- CIR_SIZE banks of RAM and serially retrieve CHK_DEG messages from the memory
			for(j = 0; j < CIR_SIZE; j++)
			{
				//-- Check node parallel process denoted in for loop 
				//-- (loop through all check nodes in one group)
				for(k = 0; k < CHK_DEG; k++)		
				{
					//-- Prepare all the V2C message
					//Shift = CodeROM.getCirShift(k, i);
					ChkAddr = j + i*CIR_SIZE;
					//VarAddr = (k + Shift) % CIR_SIZE + j*CIR_SIZE;
					BankSelect = k;
					EdgeRAM[BankSelect].setAddress(ChkAddr); 
					Reg_v2c_pre[j][k] = EdgeRAM[BankSelect].rdData();	
				}
			}
			//-- Do binary operation of the sxor: forward backward computation
			for(j = 0; j < CIR_SIZE; j++)	//-- Parallel process
			{
				Forward[j][0] = Reg_v2c_pre[j][0];
				Backward[j][CHK_DEG-1] = Reg_v2c_pre[j][CHK_DEG-1];
			}
			//-- keeping this part as it is...
			for(k = 1; k < CHK_DEG-1; k++) //-- Serial process?
			{
				for(j = 0; j < CIR_SIZE; j++)	//-- Parallel process?
				{
					Forward[j][k] = sxor(Forward[j][k-1], Reg_v2c_pre[j][k]);
					Backward[j][CHK_DEG - k - 1] 
					= sxor(Backward[j][CHK_DEG - k], Reg_v2c_pre[j][CHK_DEG - 1 - k]);
				}
			}
			//-- Compute MC2V and accumulate the posteriori:
			//-- NOTE: can ignore the buffer Reg_c2v in software
			//-- 1. Compute the first and last MC2V and write back to RAM
			for(j = 0; j < CIR_SIZE; j++)	//-- Parallel process
			{
				k = 0; // 1st
				Reg_c2v[j][k] = Backward[j][k+1];
				ChkAddr = j + i*CIR_SIZE;
				BankSelect = k;
				EdgeRAM[BankSelect].setAddress(ChkAddr); 
				EdgeRAM[BankSelect].wrtData(Reg_c2v[j][k]);

				//Reg_c2v[CHK_DEG-1] = Forward[CHK_DEG-1];
				k = CHK_DEG-1; // last
				Reg_c2v[j][k] = Forward[j][k-1];
				ChkAddr = j + i*CIR_SIZE;
				BankSelect = k;
				EdgeRAM[BankSelect].setAddress(ChkAddr); 
				EdgeRAM[BankSelect].wrtData(Reg_c2v[j][k]);
			}
			//-- 2. Compute rest of the MC2V message and write back to the RAM
			for(k = 1; k < CHK_DEG-1; k++)
			{							
				for(j = 0; j < CIR_SIZE; j++) //-- Parallel process
				{
					Reg_c2v[j][k] = sxor(Forward[j][k-1], Backward[j][k+1]);
					//VarAddr = k + i*CHK_DEG;
					//EdgeRAM[j].setAddress(VarAddr); 
					//EdgeRAM[j].wrtData(Reg_c2v[j][k]);
					ChkAddr = j + i*CIR_SIZE;
					BankSelect = k;
					EdgeRAM[BankSelect].setAddress(ChkAddr); 
					EdgeRAM[BankSelect].wrtData(Reg_c2v[j][k]);
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
					//Select which bank of RAM is active
					BankSelect = i;  // vgroup index 	
					//  need to deduce ChkAddr here: the right shift is upshift => minus
					temp = j - Shift;
					if(temp < 0)
						temp += CIR_SIZE;
					//The Address within the selected bank of RAM
					ChkAddr = temp + k*CIR_SIZE;	
					EdgeRAM[BankSelect].setAddress(ChkAddr);	
					Reg_c2v_pre[j][k] = EdgeRAM[BankSelect].rdData();
					Accum[j] = Accum[j] + Reg_c2v_pre[j][k];
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
					BankSelect = i; // the vgroup index is the BankSelect
					//  need to deduce ChkAddr here
					temp = j - Shift;
					if(temp < 0)
						temp += CIR_SIZE;
					//The Address within the selected bank of RAM
					ChkAddr = temp + k*CIR_SIZE;	
					EdgeRAM[BankSelect].setAddress(ChkAddr);	
					EdgeRAM[BankSelect].wrtData(Accum[j] - Reg_c2v_pre[j][k]);
				}
			}
		}
		Iteration++;
		// If there is no error then break (cannot detect codeword error)
		if(!checkPost_fp())
		{
			//cout << "check sum of decoder output pass at iteration "
			//	  << Iteration << endl;
			FSM.setState(IDLE);
		}
		else
		{
			FSM.setState(C2V);
			//-- for debugging
			//resetBER();
			//calculateBER();
			//cout << BitError << " bit errors in iteration " << Iteration << endl;
		}
	}// While Iteration loop end
	//checkPost_fp();
	// return total number of iterations
	return Iteration;
}


//---- only for performance test, set the information bits to calculate BER
void FP_Decoder::setInfoBit(char* in, int in_len)
{
	int i, j;
	int char_size = sizeof(char)*8;
	int counter = 0;
	for(i = 0; i < in_len-1; i++)
	{
		for(j = 0; j < char_size; j++)// per byte
		{	
			TrueInfoBit[counter] = (in[i] >> j) & 1;
			counter ++;
		}
	}
	// unload the remaining bits from the last element 
	for(j = 0; j < INFO_LENGTH % char_size; j++)
	{
		TrueInfoBit[counter] = (in[in_len-1] >> j) & 1;
		counter ++;
	}
}

void FP_Decoder::setCodeword(int* in)
{
	int i;
	for(i = 0; i < CWD_LENGTH; i++)
	{
		TrueCodeword[i] = in[i];
	}
}


//---- only perform parity check on a vector
int FP_Decoder::check_fp(int *in)
{
	int i, j, k, shift;
	unsigned int varaddr, chkaddr;
	unsigned int checksum;
	for(i = 0; i < NUM_CGRP; i ++)
	{
		for(j = 0; j < CIR_SIZE; j++)
		{
			checksum = 0;
			chkaddr = i*CIR_SIZE + j;
			for(k = 0; k < NUM_VGRP; k++)
			{
				shift = CodeROM.getCirShift(i, k);
				varaddr = j + k*NUM_VGRP + shift;
				checksum ^= in[varaddr];		
			}
			if(checksum != 0)
				return 1; //check sum no pass
		}
	}
	return 0; //check sum pass
}

int FP_Decoder::check()
{
	int i, j, k, shift;
	unsigned int varaddr, chkaddr;
	unsigned int checksum;
	
	for(i = 0; i < NUM_CGRP; i ++)
	{
		for(j = 0; j < CIR_SIZE; j++)
		{
			checksum = 0;
			chkaddr = i*CIR_SIZE + j;
			for(k = 0; k < NUM_VGRP; k++)
			{
				// key problem: differences between shift back and shift forward: my program used 
				// shift forward, but standard array code use shift back
				// before changing decoder, try changing the generation of H and G
				shift = CodeROM.getCirShift(i, k);
				//if(j - shift < 0)
				//	dummy = j - shift + CIR_SIZE;
				//else
				//	dummy = j - shift;
				//dummy = (shift + j)% CIR_SIZE;
				varaddr = (shift + j)% CIR_SIZE  + k*NUM_VGRP;
				//if(varaddr != clist[chkaddr][k])
				//	cout << "mistmatch!" << endl;
				checksum ^= TrueCodeword[varaddr];		
			}
			if(checksum != 0)
				return 1; //check sum no pass
		}
	}
	
	return 0; //check sum pass
}

int FP_Decoder::hardDecision(const int *in)
{
	int i, j, k, checksum, shift;
	unsigned int varaddr;
	for(i = 0; i < CWD_LENGTH; i++)
	{
		DecodedCodeword[i] = in[i] > 0? 0:1;
	}
	for(i = 0; i < NUM_CGRP; i ++)
	{
		for(j = 0; j < CIR_SIZE; j++)
		{
			checksum = 0;
			for(k = 0; k < NUM_VGRP; k++)
			{
				shift = CodeROM.getCirShift(i, k);
				varaddr = (shift + j)% CIR_SIZE + k*NUM_VGRP;
				checksum ^= (DecodedCodeword[varaddr]);
			}
			if(checksum != 0)
				return 1; // check sum does not pass
		}
	}
	return 0;
}
//---- use posterior to determined decoded codeword and perform check sum
int FP_Decoder::checkPost_fp()
{
	int i, j, k;
	//int temp[NUM_VAR];
	int shift;
	unsigned int varaddr;
	unsigned int checksum = 0, checksumtemp = 0;

	// full check
	for(i = 0; i < CWD_LENGTH; i++)
	{
		DecodedCodeword[i] = getPost_fp(i) > 0 ? 0:1;
	}

	for(i = 0; i < NUM_CGRP; i ++)
	{
		for(j = 0; j < CIR_SIZE; j++)
		{
			checksum = 0;
			checksumtemp = 0;
			for(k = 0; k < NUM_VGRP; k++)
			{
				shift = CodeROM.getCirShift(i, k);
				varaddr = (shift + j)% CIR_SIZE + k*NUM_VGRP;
				checksum ^= (DecodedCodeword[varaddr]);
				checksumtemp ^= (TrueCodeword[varaddr]);	
			}
			if(checksum != 0)
				return 1; // check sum does not pass
		}
	}
	return 0; // check sum pass
}


// for debugging use, just in case we need it. 
void ReadH()
{
	int i, j;
	//-------------- H matrix check
	ifstream filestr("H_array_p47_r5_forward.txt");
	int vnum, cnum, vdeg_max, cdeg_max;
	int dummy, vdeg[2209], cdeg[235], vlist[2209][5], clist[235][47];
	int par_count = 0, info_count = 0;
	filestr >> vnum;
	filestr >> cnum;
	filestr >> vdeg_max;
	filestr >> cdeg_max;
	for(i = 0; i < vnum; i++)
	{
		filestr >> vdeg[i];
	}
	for(i = 0; i < cnum; i++)
	{
		filestr >> cdeg[i];
	}
	for(i = 0; i < vnum; i++)
	{
		for(j = 0; j < vdeg[i]; j++)
			filestr >> vlist[i][j];
	}
	for(i = 0; i < cnum; i++)
	{
		for(j = 0; j < cdeg[i]; j++)
			filestr >> clist[i][j];
	}
	filestr.close();
	//-------------- H matrix check end
}

//---- no mod
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


// specify the info index
void FP_Decoder::setInfoIndex(int* in)
{
	int i = 0;
	for(i = 0; i < INFO_LENGTH; i ++)
	{
		InfoIndex[i] = in[i];
	}
}

int FP_Decoder::calculateBER()
{
	int i = 0; 
	//cout << "bit error position: ";
	for(i = 0; i < INFO_LENGTH; i++)
	{
		if(DecodedCodeword[InfoIndex[i]] != TrueInfoBit[i])
		{
			BitError ++;
			//cout << i << " "; 
		}
		//cout << endl;
	}
	return BitError;
	//cout << endl;
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


/// error version...
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

	// for debugging
	double Received[NUM_VAR];
	for(i = 0; i < NUM_VAR; i++)
	{
		Received[i] = LLR[i];
	}
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
					//Shift = CodeROM.getCirShift(i, k);
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
					if(MV2C[k] > 8)
						cout << MV2C[k];
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
					//-- old code
					//Shift = CodeROM.getCirShift(k, i);
					//ChkAddr = (j + Shift) % CIR_SIZE;
					//EdgeRAM[k].setAddress(ChkAddr + i*CIR_SIZE);	
					//MC2V[k] = EdgeRAM[k].getData();
					//temp = temp + MC2V[k];
					//-- new code
					EdgeRAM[j].setAddress(i + k*CIR_SIZE);	
					MC2V[k] = EdgeRAM[j].getData();
					temp = temp + MC2V[k];
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
	//checkPost();
	return BitError;
}
/* !!!!!!!!!
need to modify to really check the parity check matrix!!!
not urgent yet
*/ 
int FP_Decoder::checkPost()
{
	int i = 0;
	BitError = 0;
	double temp = 0;
	for(i = 0; i < CWD_LENGTH; i++)
	{
		//temp = getPost(i);
		if(getPost(i) < 0)
			BitError++;
	}
	if(BitError != 0)
		return 1;
	else
		return 0;
}
