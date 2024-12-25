#pragma once

#include "ckks_evaluator.cuh"
#include "phantom.h"

namespace nexus {
using namespace std;
using namespace phantom;

class GELUEvaluator {
 private:
  int d_g = 2;
  int d_f = 2;

  std::shared_ptr<CKKSEvaluator> ckks;

 public:
  GELUEvaluator(std::shared_ptr<CKKSEvaluator> ckks) : ckks(ckks) {}

  void gelu(PhantomCiphertext &x, PhantomCiphertext &res);
};
}  // namespace nexus
