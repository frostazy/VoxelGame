
#include "HyVoxelPrivatePCH.h"

#include "Common/RTree.h"

bool HyVoxel::RTFileStream::OpenRead(const char* a_fileName)
    {
        if (fopen_s(&m_file, a_fileName, "rb") != 0)
        {
            m_file = NULL;
        }
        if (!m_file)
        {
            return false;
        }
        return true;
    }
    
bool HyVoxel::RTFileStream::OpenWrite(const char* a_fileName)
    {
        if (fopen_s(&m_file, a_fileName, "wb") != 0)
        {
            m_file = NULL;
        }
        if (!m_file)
        {
            return false;
        }
        return true;
    }