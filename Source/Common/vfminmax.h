// Purpose of this file is to undefine platform specific min and max defines and instead use the std::min and std::max.
// This will improve cross platform compatability.

//#pragma once

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

#include <algorithm>
