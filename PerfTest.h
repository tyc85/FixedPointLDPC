#ifndef PERF_TEST_H
#define PERF_TEST_H

void noMoreMemory();
int ArrayLDPC_Debug();
int ArrayLDPC_Debug_Wifi();
int ArrayLDPC_PerfTest(double db_start, double db_end, double db_step, char* Filename);
int ArrayLDPC_TimeTrial(double db, int MaxPckNum, char* Filename);
int DecodeTrial(double EbN0_dB, int MaxPacket);
int EncodeTrial(char *info, int MaxPacket);
int ArrayLDPC_Debug_Shorten(int short_len);
#endif