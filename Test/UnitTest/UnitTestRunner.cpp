#include <windows.h>

#include <string>
#include <fstream>
#include <iostream>

using namespace std;

typedef int (*RunAllUnitTestsPtr)(int argc, char* argv[]);


bool FileExist(const string& filePath)
{
	ifstream file(filePath);

	return !file.fail();
}

int main(int argc, char* argv[])
{
	HMODULE hModule = GetModuleHandle(NULL);
	char runnerFilePath_c[MAX_PATH];
	GetModuleFileName(hModule, runnerFilePath_c, MAX_PATH);

	string runnerFilePath = runnerFilePath_c;
	string binDir;
	const size_t last_slash_idx = runnerFilePath.rfind('\\');
	if (std::string::npos != last_slash_idx)
	{
		binDir = runnerFilePath.substr(0, last_slash_idx);
	}
	
	string hyVoxelFilePath = binDir + "\\HyVoxel.dll";
	if (!FileExist(hyVoxelFilePath))
	{
		cout << "HyVoxel.dll is not found." << endl;
		return 1;
	}

	HINSTANCE hinstLib = LoadLibrary(hyVoxelFilePath.c_str());
	if (hinstLib == NULL)
	{
		cout << "Failed to load HyVoxel.dll" << endl;
		return 1;
	}

	RunAllUnitTestsPtr runAllUnitTestsPtr = (RunAllUnitTestsPtr)GetProcAddress(hinstLib, "RunAllUnitTests");
	if (runAllUnitTestsPtr == NULL)
	{
		cout << "RunAllUnitTests not found in HyVoxel.dll. Check if HyVoxel.dll is shipping version." << endl;
		return 1;
	}

	// The unit tests will be run during this query.
	return (*runAllUnitTestsPtr)(argc, argv);
}