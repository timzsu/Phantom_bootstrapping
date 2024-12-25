#pragma once

#include "ckks_evaluator.cuh"
#include "phantom.h"

namespace nexus {
using namespace std;
using namespace phantom;

class SoftmaxEvaluator {
 private:
  std::shared_ptr<CKKSEvaluator> ckks;

 public:
  SoftmaxEvaluator(std::shared_ptr<CKKSEvaluator> ckks) : ckks(ckks) {}

  void softmax(PhantomCiphertext &x, PhantomCiphertext &res, int len);
};
}  // namespace nexus
