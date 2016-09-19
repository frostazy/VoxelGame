/*
 * Configuration properties for HyVoxel library
 */

#pragma once

#include <set>
#include <map>
#include <vector>
#include <list>
#include <string>
#include <queue>
#include <deque>

#include <stdint.h>

namespace HyVoxel
{
	/************************************************************************/
	/* Memory management                                                    */
	/************************************************************************/
#define VF_NEW new
#define VF_DELETE delete

#define HYVOXEL_USE_DEFAULT_ALLOC    1
#define VF_RAWALLOC(bytes) (HyVoxel::MallocDefault(bytes))
#define VF_RAWCALLOC(num, bytes) (HyVoxel::CallocDefault(num, bytes))
#define VF_RAWREALLOC(ptr, bytes) (HyVoxel::ReallocDefault(ptr, bytes))
#define VF_ALLOC(typ, num) ((typ *)VF_RAWALLOC((num) * sizeof(typ)))
#define VF_FREE(ptr) (free((ptr)))

#define VF_RAWALLOC_ALIGNED(bytes) (HyVoxel::MallocAlignedDefault(bytes))
#define VF_ALLOC_ALIGNED(typ, num) ((typ*)VF_RAWALLOC_ALIGNED((num)*sizeof(typ)))
#define VF_FREE_ALIGNED(ptr) (HyVoxel::FreeAlignedDefault(ptr))

#if 0
	/************************************************************************/
	/* Basic types                                                          */
	/************************************************************************/

	typedef signed char int8_t;
	typedef unsigned char uint8_t;

	typedef signed short int16_t;
	typedef unsigned short uint16_t;

	typedef signed int int32_t;
	typedef unsigned int uint32_t;

#ifdef _WIN32
	typedef signed __int64 int64_t;
	typedef unsigned __int64 uint64_t;
#endif
#endif

	typedef float float32_t;
	typedef double float64_t;

	/************************************************************************/
	/* Containers                                                           */
	/************************************************************************/
	template <typename T>
	using TVFVector = std::vector<T>;

	template <typename T>
	using TVFList = std::list<T>;

	template <typename K, typename V>
	using TVFMap = std::map<K, V>;

	template <typename K>
	using TVFSet = std::set<K>;

	template <typename K>
	using TVFQueue = std::queue<K>;

	template <typename K>
	using TVFDeque = std::deque<K>;

	typedef std::string VFString;

	extern std::ofstream logstream;

	#if HYVOXEL_USE_DEFAULT_ALLOC
	static inline void* MallocDefault(size_t nBytes)
	{
		void* result = malloc(nBytes);
		if (result == nullptr)
		{
			#if !TARGET_OS_IPHONE && !defined(__ANDROID_API__)
			throw std::bad_alloc();
			#else
			printf("bad alloc");
			exit(EXIT_FAILURE);
			#endif
		}
		return result;
	}

	static inline void* CallocDefault(size_t nObjects, size_t nSize)
	{
		void* result = calloc(nObjects, nSize);
		if (result == nullptr)
		{
			#if !TARGET_OS_IPHONE && !defined(__ANDROID_API__)
			throw std::bad_alloc();
			#else
			printf("bad alloc");
			exit(EXIT_FAILURE);
			#endif
		}
		return result;
	}

	static inline void* ReallocDefault(void* ptr, size_t nBytes)
	{
		void* result = realloc(ptr, nBytes);
		if (result == nullptr)
		{
			#if !TARGET_OS_IPHONE && !defined(__ANDROID_API__)
			throw std::bad_alloc();
			#else
			printf("bad alloc");
			exit(EXIT_FAILURE);
			#endif
		}
		return result;
	}

	static inline void* MallocAlignedDefault(size_t nBytes)
	{

		#ifdef _WIN32
		void* result = _aligned_malloc(nBytes, 16);
		if (result == nullptr)
			#else
		// int posix_memalign(void **memptr, size_t alignment, size_t size);
		void* result = NULL; // was nullptr
		int memResultCode = posix_memalign(&result, 16, nBytes);
		if (result == NULL || memResultCode != 0)
			#endif
		{
			#if !TARGET_OS_IPHONE && !defined(__ANDROID_API__)
			throw std::bad_alloc();
			#else
			printf("bad alloc");
			exit(EXIT_FAILURE);
			#endif
		}
		return result;
	}

	static inline void FreeAlignedDefault(void* ptr)
	{
		#ifdef _WIN32
		_aligned_free(ptr);
		#else
		free(ptr);
		#endif

	}
	#endif

}


