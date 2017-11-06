#ifdef _MSC_VER
#include <windows.h>

class PerfChrono
{
    __int64 freq, t0; 

public:
    PerfChrono() {
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    }

    void Start() {
        QueryPerformanceCounter((LARGE_INTEGER*)&t0);
    }

    DWORD GetDiffMs(){
        __int64 t1;
        QueryPerformanceCounter((LARGE_INTEGER*)&t1); 
        return (DWORD)(((t1 - t0) * 1000) / freq); 
    }

    DWORD GetDiffUs() { //micro sec
        __int64 t1; 
        QueryPerformanceCounter((LARGE_INTEGER*)&t1); 
        return (DWORD)(((t1 - t0) * 1000000) / freq); 
    }

    DWORD GetDiffNs(){
        __int64 t1; 
        QueryPerformanceCounter((LARGE_INTEGER*)&t1); 
        return (DWORD)(((t1 - t0) * 1000000000) / freq);
    }

    DWORD GetDiff(UINT unit){
        __int64 t1;
        QueryPerformanceCounter((LARGE_INTEGER*)&t1);
        return (DWORD)(((t1 - t0) * unit) / freq);
    }

    DWORD GetFreq(){
        return (DWORD)freq;
    }
};

#else

#include <time.h>

class PerfChrono
{
    
public:
    double GetDiffMs(){
        return (clock()-start)*1000./CLOCKS_PER_SEC; 
    }
    void Start() {
        start = clock();
    }

clock_t start; 
};

#endif