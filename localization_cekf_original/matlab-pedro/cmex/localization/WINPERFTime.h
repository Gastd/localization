#include <windows.h>

typedef struct {
	LARGE_INTEGER	StartTickCount;
	LARGE_INTEGER	CurrentTickCount;
	LARGE_INTEGER	ElapsedTickCount;
	LARGE_INTEGER	TickFrequency;
} WINPERF_COUNTER,*PWINPERF_COUNTER;


/*********************************************************************************************
**********************************************************************************************
** WINPERF_ResetCounter: 
**		Reset performance counter. 
**********************************************************************************************
*********************************************************************************************/
__inline  void WINPERF_ResetCounter(PWINPERF_COUNTER pWinperfCounter)
{
	if (!QueryPerformanceCounter(&pWinperfCounter->StartTickCount)){
		pWinperfCounter->StartTickCount.QuadPart = -1;
	}
}

/*********************************************************************************************
**********************************************************************************************
** WINPERF_GetElapsedTime: 
**		Get relapsed time since last reset of the performance counter. 
**********************************************************************************************
*********************************************************************************************/
__inline double WINPERF_GetElapsedTime(PWINPERF_COUNTER pWinperfCounter)
{
	if (!QueryPerformanceCounter(&pWinperfCounter->CurrentTickCount)){
		pWinperfCounter->CurrentTickCount.QuadPart = -1;
	}
	else{
		QueryPerformanceFrequency(&pWinperfCounter->TickFrequency);
		pWinperfCounter->ElapsedTickCount.QuadPart = pWinperfCounter->CurrentTickCount.QuadPart - pWinperfCounter->StartTickCount.QuadPart;
		return( ((double)(pWinperfCounter->ElapsedTickCount.LowPart))/((double)(pWinperfCounter->TickFrequency.LowPart)) );
	}
	return -1.0;
}

