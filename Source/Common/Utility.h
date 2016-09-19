#pragma once

#include <assert.h>

namespace HyVoxel
{
    inline void Clamp(int & value, int upper, int lower)
    {
        assert(lower <= upper);

        value = (value < lower) ? lower : value;
        value = (value > upper) ? upper : value;
    }
}