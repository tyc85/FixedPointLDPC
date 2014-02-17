#ifndef ARRAY_MACRO_H
#define ARRAY_MACRO_H


#include <fstream>
#include <iostream>
#include <math.h>
//---------- Shared constants
//#define ALLZERO
//#define INFO_LENGTH 192
//#define XMITTED_BIT 384
//#define CWD_LENGTH 608
//#define CWD_LENGTH 576
//#define CIRCULANT_SIZE 32

using namespace std;
enum Simulation {MAX_ITER = 10, NUM_PEEK = 1000000, SEED = 100};
//enum CodeWifi {
//		NUM_VAR = 1944, NUM_CHK = 972, NUM_CGRP = 12, VAR_DEG = 24,
//		P = 81, CIR_SIZE = 81, INFO_LENGTH = 1978, CWD_LENGTH = 1944};

enum Code {
		NUM_VAR = 2209, NUM_CHK = 235, NUM_CGRP = 5, VAR_DEG = 5, NUM_VGRP = 47, 
		CHK_DEG = 47, P = 47, CIR_SIZE = 47, INFO_LENGTH = 1978, CWD_LENGTH = 2209};
enum RAM_Const {
		RAM_WIDTH = 32, RAM_SLICE = 8, RAM_DEPTH = CHK_DEG*VAR_DEG};
enum Precision {
		WIDTH_MASK = 0x000000ff, 
		SIGN_MASK = 0x00000080, 
		/*INT_WIDTH = 4, 
		FRAC_WIDTH = 3, 
		INT_WIDTH_NOISE = 4, 
		FRAC_WIDTH_NOISE = 12*/
		INT_WIDTH = 2, 
		FRAC_WIDTH = 4, 
		INT_WIDTH_NOISE = 4, 
		FRAC_WIDTH_NOISE = 6
};
enum StateFSM {IDLE, PCV, V2C, SXOR, C2V, SIMEND};

class ROM
{
public:
	ROM()
	{
		RowDeg = NUM_VGRP;
		ColDeg = NUM_CGRP;
		CirSize = P;
		NumVar = NUM_VAR;
		NumChk = NUM_CHK;
		int i, j;
		for(i = 0; i < NUM_CGRP; i++)
		{
			for(j = 0; j < NUM_VGRP; j++)
			{
				CirShift[i][j] = (i*j) % CirSize;
			}
		}
		CodeRate = 1 - double(NUM_CGRP*CirSize -NUM_CGRP+1)/(CirSize*CirSize);
		cout <<"Array code with circulant size "  << P <<", blocklength " 
			  << NUM_VAR << ", code rate " << getRate() << endl;
	};
	//~ROM();
	double getRate(){ return CodeRate; };
	int getCirShift(int Chk, int Var){ return CirShift[Chk][Var]; };
private:
	int RowDeg;	// Set as variable to accommodate irregular code for future
	int ColDeg;
	int CirSize;
	int Width;
	int Slice;
	int Depth;
	int Address;
	int CirShift[NUM_CGRP][NUM_VGRP];
	int DataBus;
	int WrE;
	int RdE;
	int NumVar;
	int NumChk;
	double CodeRate;
};
class Memory
{
public:
	//Memory();
	//{
	//	set_new_handler(noMoreMemory);
	//	BRAM = new int[RAM_DEPTH];//CHK_DEG*CIR_SIZE
	//};
	//~Memory(){ delete [] BRAM; };
	int rdData(){	return BRAM_fp[Address]; };
	double getData(){	return BRAM[Address]; };
	void wrtData(int in){ BRAM_fp[Address] = in; };
	void wrtData(double in){ BRAM[Address] = in; };
	void setAddress(int Addr){ Address = Addr; };
private:
	int Address;
	int BRAM_fp[RAM_DEPTH];
	double BRAM[RAM_DEPTH];
	int DataBus;
	int WrE;
	int RdE;
};


class ControlFSM
{
public:
	ControlFSM(){ CurState = IDLE; };
	// Should change it to comb. circuit afterward
	void setState(int in){ CurState = in; };	
	int getState(){ return CurState;};
private:
	int CurState;
};


class FP_Decoder
{
public:
	int InfoBit[INFO_LENGTH];
	int decode_fixpoint(const int *LLR);
	int decode(const double *LLR);
	int sgn(double);
	int sgn(int);
	int sxor(int, int);
	int fmin(int, int);
	int fmax(int, int);
	int checkPost();
	int checkPost_fp();
	int getPost_fp(int Addr){ return Posteriori_fp[Addr]; };
	int getState(){ return FSM.getState();  };
	void setState(int in){ FSM.setState(in);  };

	double fmin(double, double);
	double getPost(int Addr){ return Posteriori[Addr]; };
	double getRate(){ return CodeROM.getRate(); };
	double sxor(double, double);
	void wrtPost(int Addr, int in){ Posteriori_fp[Addr] = in;};
	void wrtPost(int Addr, double in){ Posteriori[Addr] = in;};

	//const int* getPostBlk() {return Posteriori;};
private:
	class ROM CodeROM;
	class ControlFSM FSM;
	class Memory EdgeRAM[CIR_SIZE];
	int Posteriori_fp[CWD_LENGTH];
	double Posteriori[CWD_LENGTH];
	int BitError;
	// shifting a fixed point fractional numbers to a binary representation
	static const int Constant = (5.0/8.0)*(1 << FRAC_WIDTH);
};




inline int FP_Decoder::sgn(double x){
	return (x > 0)?1:-1;
}

inline int FP_Decoder::sgn(int x){
	return (x > 0)?1:-1;
}

inline double FP_Decoder::fmin(double x, double y){		
	if(x <= y)
		return x;
	else
		return y;
}
inline int FP_Decoder::fmin(int x, int y){		
	if(x <= y)
		return x;
	else
		return y;
}

//--------------------- Encoder class
class FP_Encoder
{
public:
	FP_Encoder();
	~FP_Encoder();
	FP_Encoder(char*m, int);
	int encode(char *, char *, int);
private:
	int Codeword[NUM_VAR];
	//int ChkDeg[NUM_VAR];
	int Gdimx;
	int Gdimy;
	int InfoLen;
	int Codelen;
	class ROM CodeROM;
	unsigned int *ParityIndex;
	unsigned int *InfoIndex;
	unsigned int **GeneratorMat;
};

#endif

