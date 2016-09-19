
#include "HyJobsImpl.h"

#include "concurrentqueue.h"

#include <windows.h>
#include <thread>
#include <vector>
#include <queue>
#include <mutex>
#include <atomic>

using namespace std;

namespace HyVoxel {

class JobGroup : public JobGroupBase
{
public:
	std::atomic<unsigned int> groupId;
	std::atomic<int> refCount;

	virtual uint32_t GetGroupId() override { return groupId.load(); }
	virtual void SetGroupId(uint32_t val)
	{
		groupId.store(val);
	}

	virtual int GetRefCount() override { return refCount.load();  }
	virtual void IncRefCount() override { ++refCount; }
	virtual void DecRefCount() override { --refCount; }
};

static moodycamel::ConcurrentQueue<GeneralJob*> sWorkerThreadsJobQueue;
static moodycamel::ConcurrentQueue<GeneralJob*>	sMainThreadJobQueue;

static atomic<bool>								sIsQuiting = true;
static atomic<int>								sRunningThreadCount = 0;
static vector<HANDLE>							sThreadHandles;
int												sThreadIds[MAX_THREAD_COUNT];

static unsigned long WINAPI DoProcessJobs(LPVOID lpParameter)
{
	sRunningThreadCount++;
	int threadIdx = *(int*)lpParameter;

	while (!sIsQuiting.load())
	{
		while (true)
		{
			GeneralJob* nextJob = NULL;
			bool found = sWorkerThreadsJobQueue.try_dequeue(nextJob);
			if (!found)
				break;

			if (!nextJob->jobGroupRef.IsValid())
			{
				delete nextJob;
			}
			else if (nextJob->ReadyToStartFunc && (!(*nextJob->ReadyToStartFunc)(nextJob)))
			{
				sWorkerThreadsJobQueue.enqueue(nextJob);
				// avoid busy loop
				break;
			}
			else
			{
				nextJob->threadId = threadIdx;
				(*nextJob->JobFunc)(nextJob);
				delete nextJob;
			}
		}

		// Sleep 1 ms to avoid busy-loop.
		Sleep(1);
	}

	sRunningThreadCount--;
	return 0;
}


void JobMgr::StartWorkerThreads(int theadCount)
{
	if (!sIsQuiting.load())
		return;

	sIsQuiting.store(false);

	for (int i = 0; i < MAX_THREAD_COUNT; ++i)
		sThreadIds[i] = i;

	if (theadCount == 0)
		theadCount = std::thread::hardware_concurrency();

	for (int threadIdx = 0; threadIdx < theadCount; ++threadIdx)
	{
		HANDLE newThread = ::CreateThread(NULL, 0, DoProcessJobs, (LPVOID)&(sThreadIds[threadIdx]), 0, NULL);
		sThreadHandles.push_back(newThread);
	}
}

void JobMgr::ProcessJobsOnMainThread(int jobCount)
{
	for (int i = 0; i < jobCount; ++i)
	{
		GeneralJob* nextJob = NULL;
		if (!sMainThreadJobQueue.try_dequeue(nextJob))
			return;

		if (!nextJob->jobGroupRef.IsValid())
		{
			delete nextJob;
		}
		else
		{
			(*nextJob->JobFunc)(nextJob);
			delete nextJob;
		}
	}
}

void JobMgr::ScheduleJob(GeneralJob* job, bool isMainThreadJob)
{
	if (job->jobGroupRef.jobGroup)
		job->jobGroupRef.jobGroup->IncRefCount();

	if (isMainThreadJob)
		sMainThreadJobQueue.enqueue(job);
	else
		sWorkerThreadsJobQueue.enqueue(job);
}

void JobMgr::ClearPendingJobsAndWait()
{
	sIsQuiting.store(true);

	while (sRunningThreadCount != 0)
	{
		Sleep(1);
	}

	while (true)
	{
		GeneralJob* nextJob = NULL;
		bool found = sWorkerThreadsJobQueue.try_dequeue(nextJob);
		if (!found)
			break;

		if (nextJob->ClearInputResourcesFunc)
			(*nextJob->ClearInputResourcesFunc)(nextJob);
		delete nextJob;
	}

	while (true)
	{
		GeneralJob* nextJob = NULL;
		if (!sMainThreadJobQueue.try_dequeue(nextJob))
			break;

		if (nextJob->ClearInputResourcesFunc)
			(*nextJob->ClearInputResourcesFunc)(nextJob);
		delete nextJob;
	}

}

struct JobGroupPool
{
	mutex				accessLock;
	queue<JobGroup*>	freeJobGroups;
	vector<JobGroup*>	deallocPtrs;

	const static int BLOCK_NUMBER = 100;

	~JobGroupPool()
	{
		for (auto& ptr : deallocPtrs)
		{
			delete[] ptr;
		}
		deallocPtrs.clear();
	}

	JobGroupRef AcquireJobGroup()
	{
		lock_guard<std::mutex> guard(accessLock);

		if (freeJobGroups.empty())
		{
			JobGroup* groups = new JobGroup[BLOCK_NUMBER];
			for (int i = 0; i < BLOCK_NUMBER; ++i)
				freeJobGroups.push(&(groups[i]));

			deallocPtrs.push_back(groups);
		}

		JobGroupRef ref;
		ref.jobGroup = freeJobGroups.front();
		freeJobGroups.pop();
		return ref;
	}

	void RecycleJobGroup(JobGroupRef& ref)
	{
		lock_guard<std::mutex> guard(accessLock);
		freeJobGroups.push(static_cast<JobGroup*>(ref.jobGroup));
	}
};

static JobGroupPool sJobGroupPool;

JobGroupRef JobMgr::GenerateJobGroup()
{
	static unsigned int sGroupId = 1;
	++sGroupId;
	if (sGroupId == 0)
		sGroupId = 1;

	JobGroupRef ref = sJobGroupPool.AcquireJobGroup();
	ref.jobGroup->SetGroupId(sGroupId);
	static_cast<JobGroup*>(ref.jobGroup)->refCount = 0;
	ref.groupId = sGroupId;

	return ref;
}

void JobMgr::ReleaseJobGroup(JobGroupRef& ref)
{
	if (ref.jobGroup == NULL)
		return;

	ref.jobGroup->SetGroupId(0);
	sJobGroupPool.RecycleJobGroup(ref);
}

void JobMgr::DumpDebugInfo()
{
	std::vector<GeneralJob*> jobs;

	// Get all jobs out of the concurrent queue.
	while (true)
	{
		GeneralJob* nextJob = NULL;
		bool found = sWorkerThreadsJobQueue.try_dequeue(nextJob);
		if (found)
			jobs.push_back(nextJob);
		else
			break;
	}

	char buf[1024];
	sprintf(buf, "Worker thread job count: %d\n", (int)jobs.size());
	OutputDebugStringA(buf);

	int jobIndex = 0;
	for (auto& job : jobs)
	{
		std::string jobName = job->GetName();
		sprintf(buf, "\t[#%d] JobName: %s; GroupId1: %d, GroupId2: %d\n",
			jobIndex, jobName.c_str(), job->jobGroupRef.jobGroup->GetGroupId(), job->jobGroupRef.groupId);
		OutputDebugStringA(buf);

		sWorkerThreadsJobQueue.enqueue(job);
		++jobIndex;
	}


}

}