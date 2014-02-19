#include <iostream>
#include <fstream>
#include <bitset>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rngs.h"
#include "rvgs.h"
#include "MemoryMacro.h"
#include "Memory.h"
//#include "ArrayLDPC.h"
#include "ArrayLDPCMacro.h"
using namespace std;
//--------- notes -----------
// 1. don't really need the posterior stored. just the hard decision would be enough
// 
//---------------------------


//-------------- Encoder part starts
FP_Encoder::~FP_Encoder()
{
	;
	//int i;
	//delete [] ParityIndex;
	//delete [] InfoIndex;
	//for(i = 0; i < Gdim_col; i++ )
	//	delete [] GeneratorMat[i];
	//
	//delete [] InfoBuffer;
	//delete [] GeneratorMat;
}
FP_Encoder::FP_Encoder(char* Filename, int flag)
{
	//flag == 1 : verbose, flag == 0 supress message
	if(flag) 
		cout << "reading file " << Filename << endl;
	try{
		//ifstream FileStr("G_array_zindx.txt");
		ifstream FileStr(Filename);
		FileStr.exceptions ( ifstream::failbit | ifstream::badbit );
		if(1)// alist file format, current test mode
		{
			int vnum, cnum, vdeg_max, cdeg_max;
			int par_count = 0, info_count = 0;
			FileStr >> vnum;
			
			FileStr >> cnum;
			Gdim_col = cnum;
			Gdim_row = vnum - cnum;
			FileStr >> vdeg_max;
			FileStr >> cdeg_max;
			int i, j;
			for(i = 0; i < vnum; i++)
			{
				 FileStr >> ColumnFlag[i];
				 if(ColumnFlag[i] == 1)
				 {
					 ParityIndex[par_count] = i;
					 par_count++;
				 }
				 else
				 {
					 InfoIndex[info_count] = i;
					 info_count++;
				 }
			}
			for(i = 0; i < cnum;i++)
			{
				 FileStr >> ChkDeg[i];
			}
			//for(i = 0; i < vnum; i++)
			//{
			//	for(j = 0; j < VarDeg[i]; j++)
			//		FileStr >> dummy;
			//}
			for(i = 0; i < cnum; i++)
			{
				for(j = 0; j < ChkDeg[i]; j++)
					FileStr >> G_mlist[i][j];
			}
			//FileStr.close();
			//ofstream OutStr("G_array_zindx.txt");
			//OutStr << vnum << " ";
			//OutStr << cnum << " " << endl;
			//OutStr << vdeg_max << " ";
			//OutStr << cdeg_max << " " << endl;
			//for(i = 0; i < vnum;i++)
			//{
			//	 OutStr << ColumnFlag[i] << " ";
			//}
			//OutStr << endl;
			//for(i = 0; i < cnum;i++)
			//{
			//	 OutStr << ChkDeg[i] << " ";
			//}
			//OutStr << endl;
			//for(i = 0; i < cnum; i++)
			//{
			//	for(j = 0; j < ChkDeg[i]; j++)
			//	{
			//		OutStr << G_mlist[i][j]-1 << " ";
			//	}
			//	OutStr << endl;
			//}
			//OutStr.close();
		}
		else // binary file format (the generator matrix in binary)
		{
			//int par_num, info_num;
			int i, j;
			FileStr >> Gdim_row; // row dimension
			//info_num; //
			FileStr >> Gdim_col;
			//Gdim_row = info_num;
			//Gdim_col = par_num;
		
			//InfoIndex = new unsigned int[Gdim_row];
			//ParityIndex = new unsigned int[Gdim_col];
			//GeneratorMat = new unsigned int*[Gdim_col];
			//InfoBuffer = new unsigned int [InfoLen];
			//for(i = 0; i < Gdim_col; i++)
			//{
			//	GeneratorMat[i] = new unsigned int [Gdim_row];
			//}
			for(i = 0; i < Gdim_row; i++)
				FileStr >> InfoIndex[i];
			for(i = 0; i < Gdim_col; i++)
				FileStr >> ParityIndex[i];
			for(i = 0; i < Gdim_row; i++)
			{
				for(j = 0; j < Gdim_col; j++)
				{
					// this is the transposed version, i.e., the true parity part of a systematic generator matrix
					//FileStr >> GeneratorMat[i][j];
					;
				}
			}
		}
		FileStr.close();
	}

	catch (ifstream::failure e) {
		cout << "Exception opening/reading file " << e.what() << endl;
		system("pause");
		exit(0);
	}
	catch (std::bad_alloc& ba)
	{
		std::cerr << "bad_alloc caught: " << ba.what() << '\n';
		exit(0);
	}
	if(flag)
		cout << "encoder initialized" << endl;
}

// an overloaded function that outputs simple integer array for debugging
int FP_Encoder::encode(char *in, int in_len)
{
	int i, j, k, counter = 0;
	int out_len = 2209;  // hard coded for now
	int num_blks;
	int char_size = sizeof(char)*8;
	unsigned int temp;

	// ideally in_len is matching Gdim_col/8
	// if not matching just pad zero
	// output codeword bitwise stored in output
	// calculate the number of codeword blocks to output
	num_blks = in_len*char_size/INFO_LENGTH + (in_len*char_size % INFO_LENGTH > 0? 1:0);
	
	// already allocated in the top block. unsafe to allocate the memory in the function
	//out = new int[out_len*num_blks]; 

	//for(i = 0; i<num_blks; i++)
	//{
		// unload the information bits from the char array till the 
		// second last element
	for(i = 0; i < in_len-1; i++)
	{
		for(j = 0; j < char_size; j++)// per byte
		{	
			InfoBuffer[counter] = (in[i] >> j) & 1;
			counter ++;
		}
	}
	// unload the remaining bits from the last element 
	for(j = 0; j < Gdim_row % char_size; j++)
	{
		InfoBuffer[counter] = (in[in_len-1] >> j) & 1;
		counter ++;
	}
	//}
	// NOTE for improvement: can ignore InfoBuffer and just use systematic codewords
	for(i = 0; i < Gdim_row; i++)
	{
		Codeword[InfoIndex[i]] = InfoBuffer[i];
	}
	//c = uG : each column is a parity bit.
	/*for(i = 0; i < Gdim_col; i++)
	{
		temp = 0;
		for(j = 0; j < Gdim_row; j++)
		{
			temp ^= (GeneratorMat[j][i] & InfoBuffer[j]);
		}
		Codeword[ParityIndex[i]] = temp;
	}*/
	for(i = 0; i < Gdim_col; i++)
	{
		temp = 0;
		for(j = 0; j < ChkDeg[i]; j++)
		{
			if(ColumnFlag[G_mlist[i][j]] == 0)
			{
				temp ^= Codeword[G_mlist[i][j]];
				//cout << G_mlist[i][j] << " ";
			}
		}
		Codeword[ParityIndex[i]] = temp;
	}
	return out_len;
}

// not done yet
int FP_Encoder::encode(char *in, char *out, int in_len)
{
	int i, j, k, counter = 0;
	int num_blks;
	int out_len;
	int char_size = sizeof(char)*8;
	unsigned int temp;

	// ideally in_len is matching Gdim_col/8
	// if not matching just pad zero
	// output codeword bitwise stored in output
	// calculate the number of codeword blocks to output
	num_blks = in_len*char_size/INFO_LENGTH + (in_len*char_size % INFO_LENGTH > 0? 1:0);

	for(i = 0; i < in_len-1; i++)
	{
		for(j = 0; j < char_size; j++)// per byte
		{	
			InfoBuffer[counter] = (in[i] >> j) & 1;
			counter ++;
		}
	}
	// unload the remaining bits from the last element 
	for(j = 0; j < Gdim_row % char_size; j++)
	{
		InfoBuffer[counter] = (in[in_len-1] >> j) & 1;
		counter ++;
	}
	//}
	// NOTE for improvement: can ignore InfoBuffer and just use systematic codewords
	for(i = 0; i < Gdim_row; i++)
	{
		Codeword[InfoIndex[i]] = InfoBuffer[i];
	}
	//c = uG : each column is a parity bit.
	/*for(i = 0; i < Gdim_col; i++)
	{
		temp = 0;
		for(j = 0; j < Gdim_row; j++)
		{
			temp ^= (GeneratorMat[j][i] & InfoBuffer[j]);
		}
		Codeword[ParityIndex[i]] = temp;
	}*/
	for(i = 0; i < Gdim_col; i++)
	{
		temp = 0;
		for(j = 0; j < ChkDeg[i]; j++)
		{
			if(ColumnFlag[G_mlist[i][j]] == 0)
			{
				temp ^= Codeword[G_mlist[i][j]];
				//cout << G_mlist[i][j] << " ";
			}
		}
		Codeword[ParityIndex[i]] = temp;
	}
	out_len = num_blks*NUM_VAR/8 + (num_blks*NUM_VAR % 8 > 0? 1:0);
	//out = new char[out_len];
	// stuff the codeword into an array of char
	cout << out;
	counter = 0;
	for(i = 0; i < 276; i++)
	{
		//cout << "out[i] = " << std::bitset<8>(out[i]) << endl;
		// increase the index by one every 8 bits (one byte)
		for(j = 0 ;j < 8; j++)
		{
			//cout << 0xff;
			out[i] = out[i] ^ (Codeword[i*8 + j] << j);
			//cout << "codeword " << j << "=" << Codeword[j] << endl;
			//cout << std::bitset<8>(out[i]) << endl;
		// shift the result to the left one (stuffing in)
			//out[i] = (out[i] & 0xff) >> 1;
			//cout << std::bitset<8>(out[i]) << endl;
		}
		cout << std::bitset<8>(out[i]) << endl;
		for(j = 7; j >= 0; j--)
			cout << Codeword[i*8 + j];
		cout << endl;
	}
	for(j = 0 ;j < 2209 % 8; j++)
	{
		out[i] = out[i] ^ (Codeword[i*8 + j] << j);
	}
	//for(i)
	//for(i = 0; i < 277; i++)
	//{
	//		cout << std::bitset<8>(out[i]) << ", ";
	//		for(j = 0; j < 8; j++)
	//			cout << Codeword[i*8 + j];
	//		cout << endl;
	//}
	cout <<endl;
	// returning the padding length (or should I return the codeword length?)
	return out_len/num_blks;
}

//---------------- Ecoder part ends

//FP_Encoder::~FP_Encoder()
//{
//	;
//	//int i;
//	//delete [] ParityIndex;
//	//delete [] InfoIndex;
//	//for(i = 0; i < Gdim_col; i++ )
//	//	delete [] GeneratorMat[i];
//	//
//	//delete [] InfoBuffer;
//	//delete [] GeneratorMat;
//}
//FP_Encoder::FP_Encoder(char* Filename, int flag)
//{
//	cout << "reading file " << Filename << endl;
//	int dummy;
//	try{
//		//ifstream FileStr("G_array_zindx.txt");
//		ifstream FileStr(Filename);
//		FileStr.exceptions ( ifstream::failbit | ifstream::badbit );
//		//if(!FileStr)
//		//{
//		//	cout << "fail to read the file" <<endl;
//		//}
//		if(1)// alist file format, current test mode
//		{
//			int vnum, cnum, vdeg_max, cdeg_max;
//			int par_count = 0, info_count = 0;
//			FileStr >> vnum;
//			
//			FileStr >> cnum;
//			Gdim_col = cnum;
//			Gdim_row = vnum - cnum;
//			FileStr >> vdeg_max;
//			FileStr >> cdeg_max;
//			int i, j;
//			for(i = 0; i < vnum; i++)
//			{
//				 FileStr >> ColumnFlag[i];
//				 if(ColumnFlag[i] == 1)
//				 {
//					 ParityIndex[par_count] = i;
//					 par_count++;
//				 }
//				 else
//				 {
//					 InfoIndex[info_count] = i;
//					 info_count++;
//				 }
//			}
//			for(i = 0; i < cnum;i++)
//			{
//				 FileStr >> ChkDeg[i];
//			}
//			//for(i = 0; i < vnum; i++)
//			//{
//			//	for(j = 0; j < VarDeg[i]; j++)
//			//		FileStr >> dummy;
//			//}
//			for(i = 0; i < cnum; i++)
//			{
//				for(j = 0; j < ChkDeg[i]; j++)
//					FileStr >> G_mlist[i][j];
//			}
//			//FileStr.close();
//			//ofstream OutStr("G_array_zindx.txt");
//			//OutStr << vnum << " ";
//			//OutStr << cnum << " " << endl;
//			//OutStr << vdeg_max << " ";
//			//OutStr << cdeg_max << " " << endl;
//			//for(i = 0; i < vnum;i++)
//			//{
//			//	 OutStr << ColumnFlag[i] << " ";
//			//}
//			//OutStr << endl;
//			//for(i = 0; i < cnum;i++)
//			//{
//			//	 OutStr << ChkDeg[i] << " ";
//			//}
//			//OutStr << endl;
//			//for(i = 0; i < cnum; i++)
//			//{
//			//	for(j = 0; j < ChkDeg[i]; j++)
//			//	{
//			//		OutStr << G_mlist[i][j]-1 << " ";
//			//	}
//			//	OutStr << endl;
//			//}
//			//OutStr.close();
//		}
//		else if(flag == 0) // binary file format (the generator matrix in binary)
//		{
//			//int par_num, info_num;
//			int i, j;
//			FileStr >> Gdim_row; // row dimension
//			//info_num; //
//			FileStr >> Gdim_col;
//			//Gdim_row = info_num;
//			//Gdim_col = par_num;
//			InfoLen = Gdim_row;
//		
//			//InfoIndex = new unsigned int[Gdim_row];
//			//ParityIndex = new unsigned int[Gdim_col];
//			//GeneratorMat = new unsigned int*[Gdim_col];
//			//InfoBuffer = new unsigned int [InfoLen];
//			//for(i = 0; i < Gdim_col; i++)
//			//{
//			//	GeneratorMat[i] = new unsigned int [Gdim_row];
//			//}
//			for(i = 0; i < Gdim_row; i++)
//				FileStr >> InfoIndex[i];
//			for(i = 0; i < Gdim_col; i++)
//				FileStr >> ParityIndex[i];
//			for(i = 0; i < Gdim_row; i++)
//			{
//				for(j = 0; j < Gdim_col; j++)
//				{
//					// this is the transposed version, i.e., the true parity part of a systematic generator matrix
//					//FileStr >> GeneratorMat[i][j];
//					;
//				}
//			}
//		}
//		FileStr.close();
//	}
//
//	catch (ifstream::failure e) {
//		cout << "Exception opening/reading file " << e.what() << endl;
//		system("pause");
//		exit(0);
//	}
//	catch (std::bad_alloc& ba)
//	{
//		std::cerr << "bad_alloc caught: " << ba.what() << '\n';
//		exit(0);
//	}
//	
//	cout << "encoder initialized" << endl;
//}
//
//// an overloaded function that outputs simple integer array for debugging
//int FP_Encoder::encode(char *in, int in_len)
//{
//	int i, j, k, counter = 0;
//	int out_len = 2209;  // hard coded for now
//	int num_blks;
//	int char_size = sizeof(char)*8;
//	//unsigned int *info_temp;
//	//unsigned int *codeword;
//	unsigned int temp;
//	// ideally in_len is matching Gdim_col/8
//	// if not matching just pad zero
//	// output codeword bitwise stored in output
//	
//	// calculate the number of codeword blocks to output
//	num_blks = in_len*char_size/InfoLen + (in_len*char_size % InfoLen > 0? 1:0);
//	
//	// already allocated in the top block. unsafe to allocate the memory in the function
//	//out = new int[out_len*num_blks]; 
//
//	//for(i = 0; i<num_blks; i++)
//	//{
//		// unload the information bits from the char array till the 
//		// second last element
//	for(i = 0; i < in_len-1; i++)
//	{
//		for(j = 0; j < char_size; j++)// per byte
//		{	
//			InfoBuffer[counter] = (in[i] >> j) & 1;
//			counter ++;
//		}
//	}
//	// unload the remaining bits from the last element 
//	for(j = 0; j < Gdim_row % char_size; j++)
//	{
//		InfoBuffer[counter] = (in[in_len-1] >> j) & 1;
//		counter ++;
//	}
//	//}
//	// NOTE for improvement: can ignore InfoBuffer and just use systematic codewords
//	for(i = 0; i < Gdim_row; i++)
//	{
//		Codeword[InfoIndex[i]] = InfoBuffer[i];
//	}
//	//c = uG : each column is a parity bit.
//	/*for(i = 0; i < Gdim_col; i++)
//	{
//		temp = 0;
//		for(j = 0; j < Gdim_row; j++)
//		{
//			temp ^= (GeneratorMat[j][i] & InfoBuffer[j]);
//		}
//		Codeword[ParityIndex[i]] = temp;
//	}*/
//	for(i = 0; i < Gdim_col; i++)
//	{
//		temp = 0;
//		for(j = 0; j < ChkDeg[i]; j++)
//		{
//			if(ColumnFlag[G_mlist[i][j]] == 0)
//			{
//				temp ^= Codeword[G_mlist[i][j]];
//				//cout << G_mlist[i][j] << " ";
//			}
//		}
//		Codeword[ParityIndex[i]] = temp;
//		// self check on generator mat
//		//temp = 0;
//		//for(j = 0; j < ChkDeg[i]; j++)
//		//{
//		//	temp ^= Codeword[G_mlist[i][j]];
//		//}
//		//cout << temp;
//	}
//	return out_len;
//}
//
//// not done yet
//int FP_Encoder::encode(char *in, char *out, int in_len)
//{
//	int i, j, k, counter = 0;
//	int num_blks;
//	unsigned int *info_temp;
//	//unsigned int *codeword;
//	unsigned int temp;
//	int out_len;
//	// ideally in_len is matching Gdim_col/8
//	// if not matching just pad zero
//	// output codeword bitwise stored in output
//	info_temp = new unsigned int [InfoLen];
//
//
//	// calculate the number of codeword blocks to output
//	num_blks = in_len*sizeof(char)/InfoLen + (in_len*sizeof(char) % InfoLen > 0? 1:0);
//
//	// unload the information bits from the char array till the 
//	// second last element
//	for(i = 0; i < in_len-1; i ++)
//	{
//		for(j = 0; j < 8; j++)// per byte
//		{	
//			info_temp[counter] = (in[i] >> j) & 1;
//			counter ++;
//		}
//	}
//	// unload the remaining bits from the last element 
//	for(j = 0; j < InfoLen % 8; j++)
//	{
//		info_temp[counter] = (in[in_len-1] >> j) & 1;
//		counter ++;
//	}
//	
//	
//	// temp test: 2209 /8 = 151
//	// temp test residual: 1
//	// zero padding = 7
//	int pad_len = 8 - NUM_VAR %8;
//	// should be improved to do xoring 32 bits at a time
//	for(i = 0; i < Gdim_col; i++)
//	{
//		Codeword[InfoIndex[i]] = info_temp[i];
//	}
//	for(i = 0; i < Gdim_row; i++)
//	{
//		temp = 0;
//		for(j = 0; j < Gdim_col; j++)
//		{
//			temp ^= (GeneratorMat[j][i]&info_temp[j]);
//		}
//		Codeword[ParityIndex[i]] = temp;
//	}
//	// allocate the memory in the function or not?
//	// hard coded output length. Potentiall we should take more than one block
//	// of input streams and coded into multiple codeword blocks
//
//	// calculate how many char is needed to output the codeword
//	out_len = num_blks*NUM_VAR/8 + (num_blks*NUM_VAR % 8 > 0? 1:0);
//	out = new char[out_len];
//	// stuff the codeword into an array of char
//	for(i = 0; i < NUM_VAR; i++)
//	{
//		// increase the index by one every 8 bits (one byte)
//		out[i/8] = out[i/8] + Codeword[i];
//		// shift the result to the left one (stuffing in)
//		out[i/8] = out[i/8]<< 1;
//	}
//	// returning the padding length (or should I return the codeword length?)
//	return out_len;
//}

//---------------- Encoder part ends



