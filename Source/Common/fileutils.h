/************************************************************
* (C) Voxel Farm Inc. 2015
*/

#pragma once

namespace HyVoxel
{
	/// Returns the size in bytes of the specified file
	long getFileSize(const char *fileName);
	/// Returns whether the specified file exists
	bool fileExists(const char* fileName);
}
