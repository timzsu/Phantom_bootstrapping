#pragma once

#include "ckks_evaluator.cuh"
#include "phantom.h"

namespace nexus {
using namespace std;
using namespace phantom;

class LNEvaluator {
 private:
  std::shared_ptr<CKKSEvaluator> ckks;

 public:
  LNEvaluator(std::shared_ptr<CKKSEvaluator> ckks) : ckks(ckks) {}
  void layer_norm(PhantomCiphertext &x, PhantomCiphertext &res, int len);
};
}  // namespace nexus
