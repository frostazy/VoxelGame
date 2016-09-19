#pragma once

#include <vector>

namespace HyVoxel {

	template<typename T>
	class hyvector : public std::vector<T>
	{
	public:
		int Find(const T& val) const {
			const_iterator findRes = std::find(begin(), end(), val);
			if (findRes != end())
				return (int)(findRes - begin());
			else
				return -1;
		}

		hyvector& operator+=(const hyvector& other) {
			insert(end(), other.begin(), other.end());
			return *this;
		}

		void Append(const T* arr, int length) {
			insert(end(), arr, arr + length);
		}

	};

	template<typename T>
	class OutVectorImpl : public OutVector<T>, public hyvector<T>
	{
	public:
		virtual void Release() { delete this; }

		void SetOutVectorValues()
		{
			length = (int)this->size();
			pdata = length > 0 ? this->data() : NULL;
		}
	};

}