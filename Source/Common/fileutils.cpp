/************************************************************
* (C) Voxel Farm Inc. 2015
*/
#include "HyVoxelPrivatePCH.h"

#if defined(_WINDOWS)
	#include "windows.h"
#else
	#include <sys/types.h>
	#include <sys/stat.h>
#endif
#include "fileutils.h"
#include <stdio.h>

long HyVoxel::getFileSize(const char *fileName)
{
#if defined(_WINDOWS)
	BOOL                        fOk;
	WIN32_FILE_ATTRIBUTE_DATA   fileInfo;

	if (NULL == fileName)
	{
		return -1;
	}

	fOk = GetFileAttributesExA(fileName, GetFileExInfoStandard, (void*)&fileInfo);
	if (!fOk)
	{
		return -1;
	}

	return (long)fileInfo.nFileSizeLow;
#else // NOTE: This might work for windows but not tested for windows user permissions
	struct stat filestats;
	stat(fileName, &filestats);

	return filestats.st_size;
#endif
}

bool HyVoxel::fileExists(const char* fileName)
{
#if defined(_WINDOWS)
	bool returnValue = false;
	unsigned long attrib = GetFileAttributesA(fileName);
	if (attrib != 0xFFFFFFFF)
	{
		returnValue = true;
	}
	return returnValue;
#else // NOTE: This might work for windows but not tested for windows user permissions
	struct stat filestats;
	return stat(fileName, &filestats) == 0;
#endif
}

