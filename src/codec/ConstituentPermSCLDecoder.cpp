#include "Construct/PermSet.h"
#include "utils/Utils.hpp"
#include "math/MathUtils.h"
#include "ConstituentPermSCLDecoder.h"
#include "PermSCDecoder.h"
#include "SCDecoder.h"

#include <iostream>

using codec::ConstituentPermSCLDecoder;

ConstituentPermSCLDecoder::ConstituentPermSCLDecoder(const PolarSpecification* spec, size_t numPerms, size_t minNumECs)
	: mSpec(spec)
	, mNumLayers(utils::IntLog2(mSpec->Length))
	, mNumOperations(0)
	, mOuterDecoders(mNumLayers + 1)
	, mLLRs(mNumLayers + 1)
	, mSymbols(mNumLayers + 1)
	, mDynFrozenSymbols(mSpec->Length, 0)
{
	for (size_t layer = 0; layer <= mNumLayers; layer++) {
		size_t length = 1ull << layer;

		mOuterDecoders[layer].assign(1ull << (mNumLayers - layer), nullptr);
		mLLRs[layer] = new double[length];
		mSymbols[layer] = new bool[length];
	}

	InitOuterDecoders(numPerms, minNumECs, mNumLayers, 0);
}

ConstituentPermSCLDecoder::~ConstituentPermSCLDecoder()
{
	for (size_t layer = 0; layer <= mNumLayers; layer++) {
		delete[] mLLRs[layer];
		delete[] mSymbols[layer];

		for (auto* decoder : mOuterDecoders[layer]) {
			delete decoder;
		}
	}

	for (auto* outerSpec : mOuterSpecs) {
		delete outerSpec;
	}
}

namespace details {
	static class HardDecisionDecoder : public codec::Decoder {
	public:
		HardDecisionDecoder(bool frozen)
			: mFrozen(frozen)
		{}

		std::vector<bool> Decode(const std::vector<double>& inputLLRs) const override
		{
			if (inputLLRs.size() != 1) {
				throw std::invalid_argument("Only 1-element vectors allowed for HardDecisionDecoder");
			}
			return  std::vector<bool>({ mFrozen ? false : inputLLRs[0] < 0 });
		}

		virtual size_t NumOperations() const override
		{
			return 0;
		}

	private:
		bool mFrozen;
	};
} // namespace details

void ConstituentPermSCLDecoder::InitOuterDecoders(size_t numPerms, size_t minNumECs, size_t layer, size_t phase)
{
	if (layer == 0) {
		mOuterDecoders[layer][phase] = new details::HardDecisionDecoder(mSpec->StaticFrozen[phase] || mSpec->Dynamic.Frozen[phase]);
		std::cout << layer << " " << phase << "\n";
		return;
	}

	size_t length = 1ull << layer;
	size_t left  = length * phase;
	size_t right = length * (phase + 1);
	auto* outerSpec = new PolarSpecification(length);

	for (size_t i = left; i < right; i++) {
		outerSpec->StaticFrozen[i - left] = mSpec->StaticFrozen[i] || mSpec->Dynamic.Frozen[i];
	}
	outerSpec->Dimension = std::count(outerSpec->StaticFrozen.begin(), outerSpec->StaticFrozen.end(), false);

	auto blockSizes = construct::GetStabBlocksStructure(outerSpec->StaticFrozen);
	size_t numClasses = 1;
	for (auto b : blockSizes) {
		for (size_t j = 2; j <= b; j++) {
			numClasses *= (1ull << j) - 1;
		}
	}
	numClasses /= 3;

	if (numClasses >= minNumECs) {
		mOuterSpecs.push_back(outerSpec);
		auto perms = construct::BuildBLTAPermSet(std::min(numPerms, numClasses), layer, 0, 0, std::move(blockSizes));
		mOuterDecoders[layer][phase] = new PermSCDecoder(*outerSpec, std::move(perms));
		//mOuterDecoders[layer][phase] = new SCDecoder(outerSpec);

		std::cout << layer << " " << phase << " [ ";
		for (auto b : blockSizes) {
			std::cout << b << " ";
		}
		std::cout << "] (" << numClasses << " ECs)\n";
		return;
	}
	else {
		InitOuterDecoders(numPerms, minNumECs, layer - 1, 2 * phase);
		InitOuterDecoders(numPerms, minNumECs, layer - 1, 2 * phase + 1);
	}
}

size_t ConstituentPermSCLDecoder::NumOperations() const
{
	return mNumOperations;
}

std::vector<bool> ConstituentPermSCLDecoder::Decode(const std::vector<double>& inputLLRs) const
{
	mNumOperations = 0;
	std::fill(mDynFrozenSymbols.begin(), mDynFrozenSymbols.end(), 0);
	std::copy(inputLLRs.begin(), inputLLRs.end(), mLLRs[mNumLayers]);
	ProcessNode(mNumLayers, 0);
	return std::vector(mSymbols[mNumLayers], mSymbols[mNumLayers] + mSpec->Length);
}

void ConstituentPermSCLDecoder::ProcessNode(size_t layer, size_t phase) const
{
	if (mOuterDecoders[layer][phase] != nullptr) {
		ApplyOuterDecoder(layer, phase);
		return;
	}

	auto half = 1ull << (layer - 1);

	// Negative branch
	for (size_t i = 0; i < half; i++) {
		auto a = mLLRs[layer][i];
		auto b = mLLRs[layer][i + half];
		mLLRs[layer - 1][i] = utils::NegativeMerge(a, b);
		mNumOperations++;
	}
	ProcessNode(layer - 1, 2 * phase);

	// Positive branch
	for (size_t i = 0; i < half; i++) {
		auto a = mLLRs[layer][i];
		auto b = mLLRs[layer][i + half];
		auto u = mSymbols[layer - 1][i];
		mSymbols[layer][i] = u;
		mLLRs[layer - 1][i] = utils::PositiveMerge(a, b, u);
		mNumOperations++;
	}
	ProcessNode(layer - 1, 2 * phase + 1);

	for (size_t i = 0; i < half; i++) {
		auto u = mSymbols[layer - 1][i];
		mSymbols[layer][i] ^= u;
		mSymbols[layer][i + half] = u;
	}
}

void ConstituentPermSCLDecoder::ApplyOuterDecoder(size_t layer, size_t phase) const
{
	// TODO: handle dynamic frozen symbols

	size_t length = 1ull << layer;
	//size_t left  = length * phase;
	//size_t right = length * (phase + 1);
	
	Decoder* decoder = mOuterDecoders[layer][phase];
	std::vector<double> inputLLRs(mLLRs[layer], mLLRs[layer] + length);
	std::vector<bool> result = decoder->Decode(inputLLRs);

	mNumOperations += decoder->NumOperations();
	std::copy(result.begin(), result.end(), mSymbols[layer]);
	
	// math::VecGF2 messageVec(result);
}
