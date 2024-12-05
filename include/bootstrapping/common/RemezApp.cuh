#pragma once

#include <NTL/RR.h>

#include "MinicompFunc.cuh"
#include "MinicompRemez.cuh"
#include "PolyUpdate.cuh"

using namespace minicomp;

NTL::RR GetError(int d, NTL::RR t, bool is_first_function, int type, NTL::RR scale);
NTL::RR GetErrorCoeff(int d, NTL::RR t, std::vector<NTL::RR> &coeff, bool is_first_function, int type, NTL::RR scale);
