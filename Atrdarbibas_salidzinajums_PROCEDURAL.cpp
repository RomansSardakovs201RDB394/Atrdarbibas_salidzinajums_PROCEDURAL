#include <iostream>
#include <cstdlib>  
#include <random>
#include <iostream>
#include <windows.h>
#include <psapi.h>
#include <thread>

size_t allocatedBytesALL = 0;

size_t allocatedBytes = 0;
void* operator new(size_t size) {
    allocatedBytes = allocatedBytes + size;
    return malloc(size);
}

size_t AverageRAM = 0;
size_t MaxRAM = 0;
int RAMcounter = 0;
HANDLE process = GetCurrentProcess();
void RAMMONITOR() {
    while (true) {
        PROCESS_MEMORY_COUNTERS_EX pmc;
        GetProcessMemoryInfo(process, (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
        SIZE_T physMemUsedByMe = pmc.PrivateUsage;
        MaxRAM = pmc.PeakWorkingSetSize;
        AverageRAM = AverageRAM + physMemUsedByMe;
        RAMcounter++;
    }
}

int CPUcounter = 0;
double AverageCPU = 0.0;
void CPUMONITOR() {
    FILETIME idleTime, kernelTime, userTime;
    ULARGE_INTEGER lastIdleTime{}, lastKernelTime{}, lastUserTime{};
    double cpuUsage = 0.0;

    while (true) {
        GetSystemTimes(&idleTime, &kernelTime, &userTime);
        ULARGE_INTEGER currentIdleTime = { idleTime.dwLowDateTime, idleTime.dwHighDateTime };
        ULARGE_INTEGER currentKernelTime = { kernelTime.dwLowDateTime, kernelTime.dwHighDateTime };
        ULARGE_INTEGER currentUserTime = { userTime.dwLowDateTime, userTime.dwHighDateTime };

        if (lastIdleTime.QuadPart != 0) {
            ULONGLONG idleDiff = currentIdleTime.QuadPart - lastIdleTime.QuadPart;
            ULONGLONG kernelDiff = currentKernelTime.QuadPart - lastKernelTime.QuadPart;
            ULONGLONG userDiff = currentUserTime.QuadPart - lastUserTime.QuadPart;

            ULONGLONG totalDiff = kernelDiff + userDiff;
            cpuUsage = (1.0 - (double)idleDiff / (double)totalDiff) * 100.0;
        }

        lastIdleTime = currentIdleTime;
        lastKernelTime = currentKernelTime;
        lastUserTime = currentUserTime;

        AverageCPU = AverageCPU + cpuUsage;
        CPUcounter++;
        Sleep(1000);
    }
}

void CalculateAverageCPU() {
    std::cout << "Average CPU usage: " << (AverageCPU / CPUcounter) << '%' << std::endl;
}

void CalculateAverageRAM() {
    std::cout << "Average RAM usage: " << (AverageRAM / RAMcounter) << " bytes" << std::endl;
}


const int ARRAY_SIZE = 100000;
int RandomArray[ARRAY_SIZE];

const int FirstSize = 1000;
const int SecondSize = 10000;
const int ThirdSize = 100000;

void FillRandomArray() {

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uni(1, 100);

    for (int i = 0; i < ARRAY_SIZE; i++) {
        int random_integer = uni(rng);
        RandomArray[i] = random_integer;
    }

}

struct SelectionSort {

    int FirstSelectionSort[FirstSize];

    int SecondSelectionSort[SecondSize];

    int ThirdSelectionSort[ThirdSize];

    void FillFirstArray() {
        for (int i = 0; i < FirstSize; i++) {
            FirstSelectionSort[i] = RandomArray[i];
        }
    }

    void FillSecondArray() {
        for (int i = 0; i < SecondSize; i++) {
            SecondSelectionSort[i] = RandomArray[i];
        }
    }

    void FillThirdArray() {
        for (int i = 0; i < ThirdSize; i++) {
            ThirdSelectionSort[i] = RandomArray[i];
        }
    }

}SelectionSort;
typedef struct SelectionSort SelectiveInstance;

void selectionSort(int firstarray[], int secondarray[], int thirdarray[]) {
    PROCESS_MEMORY_COUNTERS_EX pmc{};
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
    SIZE_T physMemUsedByMe = pmc.WorkingSetSize;

    for (int i = 0; i < FirstSize - 1; i++) {
        int minIndex = i;
        for (int j = i + 1; j < FirstSize; j++) {
            if (firstarray[j] < firstarray[minIndex]) {
                minIndex = j;
            }
        }
        int temp = firstarray[i];
        firstarray[i] = firstarray[minIndex];
        firstarray[minIndex] = temp;
    }

    std::cout << "First SelectionSort done" << std::endl;

    for (int i = 0; i < SecondSize - 1; i++) {
        int minIndex = i;
        for (int j = i + 1; j < SecondSize; j++) {
            if (secondarray[j] < secondarray[minIndex]) {
                minIndex = j;
            }
        }
        int temp = secondarray[i];
        secondarray[i] = secondarray[minIndex];
        secondarray[minIndex] = temp;
    }

    std::cout << "Second SelectionSort done" << std::endl;

    for (int i = 0; i < ThirdSize - 1; i++) {
        int minIndex = i;
        for (int j = i + 1; j < ThirdSize; j++) {
            if (thirdarray[j] < thirdarray[minIndex]) {
                minIndex = j;
            }
        }
        int temp = thirdarray[i];
        thirdarray[i] = thirdarray[minIndex];
        thirdarray[minIndex] = temp;
    }

    std::cout << "Third SelectionSort done" << std::endl;

    std::cout << "SelectionSort RAM USAGE: " << physMemUsedByMe << std::endl;
    std::cout << "SelectionSort ALLOCATED BYTES: " << allocatedBytes << std::endl;

}

struct BubbleSort {

    int FirstBubbleSort[FirstSize];

    int SecondBubbleSort[SecondSize];

    int ThirdBubbleSort[ThirdSize];

    void FillFirstArray() {
        for (int i = 0; i < FirstSize; i++) {
            FirstBubbleSort[i] = RandomArray[i];
        }
    }

    void FillSecondArray() {
        for (int i = 0; i < SecondSize; i++) {
            SecondBubbleSort[i] = RandomArray[i];
        }
    }

    void FillThirdArray() {
        for (int i = 0; i < ThirdSize; i++) {
            ThirdBubbleSort[i] = RandomArray[i];
        }
    }

}BubbleSort;
typedef struct BubbleSort BubbleSortInstance;

void bubbleSort(int firstarray[], int secondarray[], int thirdarray[]) {
    PROCESS_MEMORY_COUNTERS_EX pmc{};
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
    SIZE_T physMemUsedByMe = pmc.WorkingSetSize;

    for (int i = 0; i < FirstSize - 1; i++) {
        for (int j = 0; j < FirstSize - i - 1; j++) {
            if (firstarray[j] > firstarray[j + 1]) {
                int temp = firstarray[j];
                firstarray[j] = firstarray[j + 1];
                firstarray[j + 1] = temp;
            }
        }
    }

    std::cout << "First BubbleSort done" << std::endl;

    for (int i = 0; i < SecondSize - 1; i++) {
        for (int j = 0; j < SecondSize - i - 1; j++) {
            if (secondarray[j] > secondarray[j + 1]) {
                int temp = secondarray[j];
                secondarray[j] = secondarray[j + 1];
                secondarray[j + 1] = temp;
            }
        }
    }

    std::cout << "Second BubbleSort done" << std::endl;

    for (int i = 0; i < ThirdSize - 1; i++) {
        for (int j = 0; j < ThirdSize - i - 1; j++) {
            if (thirdarray[j] > thirdarray[j + 1]) {
                int temp = thirdarray[j];
                thirdarray[j] = thirdarray[j + 1];
                thirdarray[j + 1] = temp;
            }
        }
    }

    std::cout << "Third BubbleSort done" << std::endl;

    std::cout << "BubbleSort RAM USAGE: " << physMemUsedByMe << std::endl;
    std::cout << "BubbleSort ALLOCATED BYTES: " << allocatedBytes << std::endl;
}

struct CountingSort {

    int FirstCountingSort[FirstSize];

    int SecondCountingSort[SecondSize];

    int ThirdCountingSort[ThirdSize];

    void FillFirstArray() {
        for (int i = 0; i < FirstSize; i++) {
            FirstCountingSort[i] = RandomArray[i];
        }
    }

    void FillSecondArray() {
        for (int i = 0; i < SecondSize; i++) {
            SecondCountingSort[i] = RandomArray[i];
        }
    }

    void FillThirdArray() {
        for (int i = 0; i < ThirdSize; i++) {
            ThirdCountingSort[i] = RandomArray[i];
        }
    }

}CountingSort;
typedef struct CountingSort CountingSortInstance;

void countingSort(int firstarray[], int secondarray[], int thirdarray[]) {
    PROCESS_MEMORY_COUNTERS_EX pmc{};
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
    SIZE_T physMemUsedByMe = pmc.WorkingSetSize;

    int firstcount[100 + 1] = { 0 };
    int* firstoutput = new int[FirstSize];

    for (int i = 0; i < FirstSize; i++)
        firstcount[firstarray[i]]++;

    for (int i = 1; i <= 100; i++)
        firstcount[i] += firstcount[i - 1];

    for (int i = FirstSize - 1; i >= 0; i--) {
        firstoutput[firstcount[firstarray[i]] - 1] = firstarray[i];
        firstcount[firstarray[i]]--;
    }

    for (int i = 0; i < FirstSize; i++)
        firstarray[i] = firstoutput[i];

    std::cout << "First CountingSort done" << std::endl;

    int secondcount[100 + 1] = { 0 };
    int* secondoutput = new int[SecondSize];

    for (int i = 0; i < SecondSize; i++)
        secondcount[secondarray[i]]++;

    for (int i = 1; i <= 100; i++)
        secondcount[i] += secondcount[i - 1];

    for (int i = SecondSize - 1; i >= 0; i--) {
        secondoutput[secondcount[secondarray[i]] - 1] = secondarray[i];
        secondcount[secondarray[i]]--;
    }

    for (int i = 0; i < SecondSize; i++)
        secondarray[i] = secondoutput[i];

    std::cout << "Second CountingSort done" << std::endl;

    int thirdcount[100 + 1] = { 0 };
    int* thirdoutput = new int[ThirdSize];

    for (int i = 0; i < ThirdSize; i++)
        thirdcount[thirdarray[i]]++;

    for (int i = 1; i <= 100; i++)
        thirdcount[i] += thirdcount[i - 1];

    for (int i = ThirdSize - 1; i >= 0; i--) {
        thirdoutput[thirdcount[thirdarray[i]] - 1] = thirdarray[i];
        thirdcount[thirdarray[i]]--;
    }

    for (int i = 0; i < ThirdSize; i++)
        thirdarray[i] = thirdoutput[i];

    std::cout << "Third CountingSort done" << std::endl;

    std::cout << "CountingSort RAM USAGE: " << physMemUsedByMe << std::endl;
    std::cout << "CountingSort ALLOCATED BYTES: " << allocatedBytes << std::endl;

}




int main()
{
    std::thread CpuMonitorThread(CPUMONITOR);
    std::thread ramMonitorThread(RAMMONITOR);

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 1; i <= 10; i++) {
        FillRandomArray();

        std::cout << "----------SelectiveSort----------" << std::endl;
        allocatedBytesALL = allocatedBytesALL + allocatedBytes;
        allocatedBytes = 0;
        SelectiveInstance* Selective = new SelectiveInstance();
        Selective->FillFirstArray();
        Selective->FillSecondArray();
        Selective->FillThirdArray();
        selectionSort(Selective->FirstSelectionSort, Selective->SecondSelectionSort, Selective->ThirdSelectionSort);


        std::cout << "----------BubbleSort----------" << std::endl;
        allocatedBytesALL = allocatedBytesALL + allocatedBytes;
        allocatedBytes = 0;
        BubbleSortInstance* BubbleSort = new BubbleSortInstance();
        BubbleSort->FillFirstArray();
        BubbleSort->FillSecondArray();
        BubbleSort->FillThirdArray();
        bubbleSort(BubbleSort->FirstBubbleSort, BubbleSort->SecondBubbleSort, BubbleSort->ThirdBubbleSort);


        std::cout << "----------CountingSort----------" << std::endl;
        allocatedBytesALL = allocatedBytesALL + allocatedBytes;
        allocatedBytes = 0;
        CountingSortInstance* CountingSort = new CountingSortInstance();
        CountingSort->FillFirstArray();
        CountingSort->FillSecondArray();
        CountingSort->FillThirdArray();
        countingSort(CountingSort->FirstCountingSort, CountingSort->SecondCountingSort, CountingSort->ThirdCountingSort);

    }
    allocatedBytesALL = allocatedBytesALL + allocatedBytes;

    std::cout << "----------RESULTS----------" << std::endl;
    CalculateAverageRAM();
    CalculateAverageCPU();
    std::cout << "All Time ALLOCATED: " << allocatedBytesALL << " bytes" << std::endl;
    std::cout << "PEAK RAM USAGE: " << MaxRAM << " bytes" << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> endClock = end - start;
    std::cout << "Execution time: " << endClock.count() << " seconds";
    std::cout << std::endl;

    CloseHandle(process);
    exit(-1);
    return 0;
}
