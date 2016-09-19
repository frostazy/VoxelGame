#pragma once
#ifndef EXTERNAL_MUTEX_H
#define EXTERNAL_MUTEX_H

#include <stdlib.h>

// some helpful macros
#define DECLARE_NOCOPY(ClassName) private: ClassName(const ClassName &); ClassName &operator=(const ClassName &)
#define CHECKSTATIC(expr,msg) typedef char ERROR_##msg[(expr) ? 1 : -1]

// some constant definitions
#define DEFAULT_SPIN_COUNT 4000

namespace ExternalMutex
{

	class Mutex
	{
		DECLARE_NOCOPY(Mutex);

	public:
		Mutex(int spinCount = DEFAULT_SPIN_COUNT);    // 0 = don't ever spin (some platforms may not implement spin feature)
		Mutex(const char* name);
		~Mutex();

		void Lock() const;
		void Unlock() const;
		bool TryLock() const;     // same as Lock(), except returns false instead of blocking if not available

	protected:
		// work-space for implementation pre-allocated to avoid allocations
		// temporary work-space needed by windows to hold its CRITICAL_SECTION object
		// normally we would just allocate this on the heap, but we needed to be able to use
		// the critical section object inside of a custom memory manager, which meant that we
		// had to have the critical section before doing any dynamic allocations (chicken and
		// egg thing).  This will also have better performance since it avoids an allocation.

		// note: we make the workspace an array of 6 unsigned 64s's instead of 48 char's, because that will
		// ensure that the start of our m_workSpaceBuffer is at least 8-byte aligned.  It turns out
		// the PS3 required this for the data it was attempting to stick in there
		mutable unsigned long long m_workSpaceBuffer[15];
	};

// this helper-class is used to enter/leave a mutex based on stack scoping
// use this instead of the manual locks and unlocks for higher maintainability
	class MutexGuard
	{
		DECLARE_NOCOPY(MutexGuard);

	public:
		MutexGuard(const Mutex *mutex); // Lock the mutex on construction
		~MutexGuard();                  // Unlock the mutex on destruction
		void UnlockEarly();             // Unlock the mutex earlier than destruction (leave will then not occur during destructor)

	private:
		const Mutex *m_mutex;
	};

	inline MutexGuard::MutexGuard(const Mutex *mutex)
	{
		// assert that mutex != null?
		if (mutex != NULL)
		{
			m_mutex = mutex;
			m_mutex->Lock();
		}
	}

	inline MutexGuard::~MutexGuard()
	{
		if (m_mutex != NULL)
		{
			m_mutex->Unlock();
			m_mutex = NULL;
		}
	}

	inline void MutexGuard::UnlockEarly()
	{
		if (m_mutex != NULL)
		{
			m_mutex->Unlock();
			m_mutex = NULL;
		}
	}

} // end of ExternalMutex namespace

#endif  // EXTERNAL_MUTEX_H