#pragma once

#include "HyVoxelLib.h"

namespace HyVoxel {

	class JobMgr : public IHyJobMgrInterface
	{
	public:
		virtual void Release() override { delete this; }

		// 0: automatically decide how many threads are suitable.
		virtual void StartWorkerThreads(int theadCount = 0) override;

		// Please DONOT schedule new jobs on other threads (worker threads are fine) after you call ClearPendingJobsAndWait().
		virtual void ProcessJobsOnMainThread(int jobCount = 1) override;

		virtual void ScheduleJob(GeneralJob* job, bool isMainThreadJob = false) override;

		virtual void ClearPendingJobsAndWait() override;

		virtual JobGroupRef GenerateJobGroup() override;
		virtual void ReleaseJobGroup(JobGroupRef& ref) override;

		virtual void DumpDebugInfo() override;
	};

}