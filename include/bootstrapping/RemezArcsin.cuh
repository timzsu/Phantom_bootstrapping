#pragma once

#include "common/Remez.cuh"

class RemezArcsin : public boot::Remez {
 public:
  RemezArcsin(RemezParam _param, double _log_width, long _deg)
      : boot::Remez(_param, 1, _log_width, _deg) {}

  NTL::RR function_value(NTL::RR x) {
    return 1 / (2 * NTL::ComputePi_RR()) * arcsin(x);
  }
};
