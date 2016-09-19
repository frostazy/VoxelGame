// Windows Implementation of Mutex
#include "HyVoxelPrivatePCH.h"


#include "ExternalMutex.h"

namespace ExternalMutex
{
	enum { cNamedID = 'NAME' };

#if defined(_WINDOWS)
#define SPINLOCK WindowsSpinLock
#define MUTEXLOCK WindowsMutexLock

#elif defined (__APPLE__)
#define SPINLOCK AppleSpinLock
#define MUTEXLOCK PthreadMutexLock

#elif defined(__ANDROID_API__)
#define SPINLOCK AndroidSpinLock
#define MUTEXLOCK PthreadMutexLock

#endif



#if defined(__ANDROID_API__)
	// define a spinlock for andoid platforms
#include <atomic>
#include <linux/spinlock.h>

	/*
	** Spin lock
	*/
	class AndroidSpinLock
	{
	public:
		static inline void init(spinlock_t& lock, int spinCount = 0)
		{
			lock = SPIN_LOCK_UNLOCKED;
			spin_lock_init(&lock);
		}

		static inline void kill(spinlock_t& lock)
		{
		}

		static inline bool tryLock(spinlock_t& lock)
		{
			return spin_trylock(&lock) != 0;
		}

		static inline void lock(spinlock_t& lock)
		{
			spin_lock(&lock);
		}

		static inline void unlock(spinlock_t& lock)
		{
			spin_unlock(&lock);
		}
	};


#endif


#if defined(_WINDOWS)

#define _WIN32_WINNT 0x403      // needed to access spin count functionality of critical section
#include "windows.h"

	/*
	** Spin lock
	*/
	class WindowsSpinLock
	{
	public:
		static inline void init(CRITICAL_SECTION &lock, int spinCount = 0)
		{
			InitializeCriticalSectionAndSpinCount(&lock, spinCount);
		}

		static inline void kill(CRITICAL_SECTION &lock)
		{
			DeleteCriticalSection(&lock);
		}

		static inline bool tryLock(CRITICAL_SECTION &lock)
		{
			return (1 == TryEnterCriticalSection(&lock));
		}

		static inline void lock(CRITICAL_SECTION &lock)
		{
			EnterCriticalSection(&lock);
		}

		static inline void unlock(CRITICAL_SECTION &lock)
		{
			LeaveCriticalSection(&lock);
		}
	};

	/*
	** Mutex lock
	*/
	class WindowsMutexLock
	{
	public:
		static inline void init(HANDLE &handle, const char *name)
		{
			handle = CreateMutexA(0, false, name);
		}

		static inline void kill(HANDLE &handle)
		{
			CloseHandle(handle);
		}

		static inline bool tryLock(HANDLE &handle)
		{
			return (0 == WaitForSingleObject(handle, 0));
		}

		static inline void lock(HANDLE &handle)
		{
			WaitForSingleObject(handle, INFINITE);
		}

		static inline void unlock(HANDLE &handle)
		{
			ReleaseMutex(handle);
		}
	};

#endif


#if defined (__APPLE__)

	#include <libkern/OSAtomic.h>

	/*
	** Spin lock
	*/
	class AppleSpinLock
	{
	public:
		static inline void init(OSSpinLock& lock, int spinCount = 0)
		{
			lock = OS_SPINLOCK_INIT;
		}

		static inline void kill(OSSpinLock& lock)
		{
            OSSpinLockUnlock(&lock);
            lock = OS_SPINLOCK_INIT;
		}

		static inline bool tryLock(OSSpinLock& lock)
		{
			return OSSpinLockTry(&lock);
		}

		static inline void lock(OSSpinLock& lock)
		{
			OSSpinLockLock(&lock);
		}

		static inline void unlock(OSSpinLock& lock)
		{
			OSSpinLockUnlock(&lock);
		}
	};
#endif

#if defined(__GNUC__)

	#include <pthread.h>

	/*
	** Mutex lock
	*/
	class PthreadMutexLock
	{
	public:
		static inline void init(pthread_mutex_t& handle, const char *name)
		{
			pthread_mutex_init(&handle, NULL);
		}

		static inline void kill(pthread_mutex_t& handle)
		{
			pthread_mutex_destroy(&handle);
		}

		static inline bool tryLock(pthread_mutex_t& handle)
		{
			return (pthread_mutex_trylock(&handle) == 0);
		}

		static inline void lock(pthread_mutex_t& handle)
		{
			pthread_mutex_lock(&handle);
		}

		static inline void unlock(pthread_mutex_t& handle)
		{
			pthread_mutex_unlock(&handle);
		}
	};
#endif





	// lock storage
	union MutexUnion
	{
		unsigned int id;

#if defined(_WINDOWS)
		CRITICAL_SECTION spinLock;
		HANDLE mutexLock;
#elif defined(__ANDROID_API__)
		spinlock_t spinLock;
		pthread_mutex_t mutexLock;
#elif defined (__APPLE__)
		OSSpinLock spinLock;
		pthread_mutex_t mutexLock;
#endif
	};


	// implementation of lock manipulation functions

	Mutex::Mutex(int spinCount)
	{
		CHECKSTATIC(sizeof(MutexUnion) <= sizeof(m_workSpaceBuffer), WorkSpaceBufferTooSmall);
		MutexUnion& workspace = (MutexUnion&)m_workSpaceBuffer;
		SPINLOCK::init(workspace.spinLock, spinCount);
	}

	Mutex::Mutex(const char *name)
	{
		if (name != NULL)       // if you use this constructor, you have to give it a name
		{
			MutexUnion& workspace = (MutexUnion&)m_workSpaceBuffer;
			workspace.id = cNamedID;
			MUTEXLOCK::init(workspace.mutexLock, name);
		}
	}

	Mutex::~Mutex()
	{
		MutexUnion& workspace = (MutexUnion&)m_workSpaceBuffer;
		if (workspace.id == cNamedID)
			MUTEXLOCK::kill(workspace.mutexLock);
		else
			SPINLOCK::kill(workspace.spinLock);		
	}

	bool Mutex::TryLock() const
	{
		MutexUnion& workspace = (MutexUnion&)m_workSpaceBuffer;
		if (workspace.id == cNamedID)
			return MUTEXLOCK::tryLock(workspace.mutexLock);
		else
			return SPINLOCK::tryLock(workspace.spinLock);
	}

	void Mutex::Lock() const
	{
		MutexUnion& workspace = (MutexUnion&)m_workSpaceBuffer;
		if (workspace.id == cNamedID)
			MUTEXLOCK::lock(workspace.mutexLock);
		else
			SPINLOCK::lock(workspace.spinLock);
	}

	void Mutex::Unlock() const
	{
		MutexUnion& workspace = (MutexUnion&)m_workSpaceBuffer;
		if (workspace.id == cNamedID)
			MUTEXLOCK::unlock(workspace.mutexLock);
		else
			SPINLOCK::unlock(workspace.spinLock);
	}

} // end of namespace ExternalMutex