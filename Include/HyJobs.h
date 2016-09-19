#pragma once

#include <cstdint>

namespace HyVoxel {

#pragma pack(push, 8)

struct GeneralJob;

#define MAX_THREAD_COUNT 32

typedef void (*JOBFUNC)(GeneralJob* job);

typedef bool (*READYTOSTARTFUNC)(GeneralJob* job);

typedef void (*CLEARRESOURCESFUNC)(GeneralJob* job);


class JobGroupBase
{
public:
	virtual ~JobGroupBase() {}

	virtual uint32_t GetGroupId() = 0;
	virtual void SetGroupId(uint32_t val) = 0;

	virtual int GetRefCount() = 0;
	virtual void IncRefCount() = 0;
	virtual void DecRefCount() = 0;
};

struct JobGroupRef
{
	JobGroupBase*	jobGroup;
	unsigned int	groupId;

	JobGroupRef()
		: jobGroup(NULL)
		, groupId(0)
	{ }

	JobGroupRef(const JobGroupRef& other)
		: jobGroup(other.jobGroup)
		, groupId(other.groupId)
	{ }

	JobGroupRef& operator=(const JobGroupRef& other)
	{
		this->jobGroup = other.jobGroup;
		this->groupId = other.groupId;
		return *this;
	}

	// So that JobGroupRef can be used in a std::set
	bool operator<(const JobGroupRef& other) const { return groupId < other.groupId; }

	bool IsValid()
	{
		if (jobGroup == NULL)
			return true;
		else
			return jobGroup->GetGroupId() == groupId;
	}
};

struct GeneralJob
{
	JOBFUNC				JobFunc;
	READYTOSTARTFUNC	ReadyToStartFunc;
	CLEARRESOURCESFUNC	ClearInputResourcesFunc;

	JobGroupRef			jobGroupRef;
	unsigned int		threadId;

	GeneralJob()
		: JobFunc(NULL)
		, ReadyToStartFunc(NULL)
		, ClearInputResourcesFunc(NULL)
		, threadId(0)
	{ }

	virtual const char* GetName() { return NULL; }

	virtual ~GeneralJob()
	{
		if (jobGroupRef.jobGroup)
		{
			jobGroupRef.jobGroup->DecRefCount();
		}
	}
};

class IHyJobMgrInterface : public IHyVoxelInterfaceBase
{
public:
	// 0: automatically decide how many threads are suitable.
	virtual void StartWorkerThreads(int theadCount = 0) = 0;

	// Please DONOT schedule new jobs on other threads (worker threads are fine) after you call ClearPendingJobsAndWait().
	virtual void ProcessJobsOnMainThread(int jobCount = 1) = 0;

	virtual void ScheduleJob(GeneralJob* job, bool isMainThreadJob = false) = 0;

	virtual void ClearPendingJobsAndWait() = 0;

	virtual JobGroupRef GenerateJobGroup() = 0;
	virtual void ReleaseJobGroup(JobGroupRef& ref) = 0;

	virtual void DumpDebugInfo() = 0;
};

#pragma pack(pop)

}