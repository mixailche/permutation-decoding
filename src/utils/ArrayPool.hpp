#pragma once

#include <vector>
#include <memory>
#include <stack>

#include "Utils.hpp"

namespace utils {
	template <typename T>
	class ArrayPool {
	public:
		ArrayPool(size_t numLayers, size_t maxPaths);
		~ArrayPool() = default;

		T* PopArray(size_t layer, size_t pathIdx);
		T* GetArray(size_t layer, size_t pathIdx) const;

		void Clear(size_t layer);
		void ForkPath(size_t srcIdx, size_t dstIdx);

	private:
		using Pointer = T*;

		size_t mNumLayers;
		size_t mMaxPaths;

		std::unique_ptr<char[]> mBuffer;
		std::vector<Pointer> mLayerPointers;
		
		std::vector<Pointer> mNewArrays;
		std::vector<Pointer> mArrayPointers;

		//std::vector<std::vector<size_t>> mFreeIndices;
		std::vector<size_t> mStackSizes;

		Pointer PopFreePointer(size_t layer);
	};
}

template<typename T>
utils::ArrayPool<T>::ArrayPool(size_t numLayers, size_t maxPaths)
	: mNumLayers(numLayers)
	, mMaxPaths(maxPaths)
	, mLayerPointers(numLayers + 1)
	, mBuffer(std::make_unique<char[]>(8 * (1ull << numLayers) * maxPaths * sizeof(T)))
	, mArrayPointers((numLayers + 1) * maxPaths, nullptr)
	, mNewArrays((numLayers + 1) * maxPaths, nullptr)
	//, mFreeIndices(numLayers + 1, utils::Iota(maxPaths))
	, mStackSizes(numLayers + 1, 4 * maxPaths)
{
	mLayerPointers[0] = (T*) mBuffer.get();

	for (size_t layer = 1; layer <= mNumLayers; layer++) {
		mLayerPointers[layer] = mLayerPointers[layer - 1] + 4 * (1ull << (layer - 1)) * maxPaths;
	}
}

template<typename T>
T* utils::ArrayPool<T>::PopFreePointer(size_t layer)
{
	return mLayerPointers[layer] + (--mStackSizes[layer]) * (1ull << layer);
}

template<typename T>
T* utils::ArrayPool<T>::PopArray(size_t layer, size_t pathIdx)
{
	return mArrayPointers[pathIdx * (mNumLayers + 1) + layer] = PopFreePointer(layer);
}

template<typename T>
T* utils::ArrayPool<T>::GetArray(size_t layer, size_t pathIdx) const
{
	return mArrayPointers[pathIdx * (mNumLayers + 1) + layer];
}

template<typename T>
void utils::ArrayPool<T>::Clear(size_t layer)
{
	mStackSizes[layer] = 4 * mMaxPaths;
}

template<typename T>
void utils::ArrayPool<T>::ForkPath(size_t srcIdx, size_t dstIdx)
{
	Pointer* curPointers = mArrayPointers.data() + srcIdx * (mNumLayers + 1);
	Pointer* newPointers = mArrayPointers.data() + dstIdx * (mNumLayers + 1);
	std::memcpy(newPointers, curPointers, (mNumLayers + 1) * sizeof(Pointer));
}
