#ifndef PERF_TEST_H
#define PERF_TEST_H

void noMoreMemory();
int ArrayLDPC_PerfTest(double db_start, double db_end, double db_step, char* Filename);
int ArrayLDPC_TimeTrial(double db, int MaxPckNum, char* Filename);

#endif