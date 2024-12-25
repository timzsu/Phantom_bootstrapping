#pragma once

#include "bootstrapping/Bootstrapper.cuh"
#include "ckks_evaluator.cuh"
#include "phantom.h"

namespace nexus {
using namespace std;
using namespace phantom;

class ArgmaxEvaluator {
 private:
  std::shared_ptr<CKKSEvaluator> ckks;
  std::shared_ptr<Bootstrapper> bootstrapper;

 public:
  ArgmaxEvaluator(std::shared_ptr<CKKSEvaluator> ckks, std::shared_ptr<Bootstrapper> bootstrapper, int main_mod_count) {
    this->ckks = ckks;
    this->bootstrapper = bootstrapper;
  }

  void argmax(PhantomCiphertext &x, PhantomCiphertext &res, int len);

  void bootstrap(PhantomCiphertext &x);
};
}  // namespace nexus
