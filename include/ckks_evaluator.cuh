#pragma once

#include <complex>
#include <cuComplex.h>

#include "phantom.h"

namespace nexus {

inline std::vector<cuDoubleComplex> pack_complex(std::vector<std::complex<double>>& values) {
  std::vector<cuDoubleComplex> cu_complex_values(values.size());
  for (size_t i = 0; i < values.size(); i++) {
    cu_complex_values[i] = make_cuDoubleComplex(values[i].real(), values[i].imag());
  }
  return cu_complex_values;
}
class Encoder {
 private:
  std::shared_ptr<PhantomContext> context;
  std::shared_ptr<PhantomCKKSEncoder> encoder;

 public:
  Encoder() = default;

  Encoder(std::shared_ptr<PhantomContext> context, std::shared_ptr<PhantomCKKSEncoder> encoder) {
    this->context = context;
    this->encoder = encoder;
  }

  inline size_t slot_count() { return encoder->slot_count(); }


  // Vector (of doubles or complexes) inputs
  inline void encode(std::vector<double> values, size_t chain_index, double scale, PhantomPlaintext &plain) {
    if (values.size() == 1) {
      encode(values[0], chain_index, scale, plain);
      return;
    }
    values.resize(encoder->slot_count(), 0.0);
    encoder->encode(*context, values, scale, plain, chain_index);
  }

  inline void encode(std::vector<double> values, double scale, PhantomPlaintext &plain) {
    if (values.size() == 1) {
      encode(values[0], scale, plain);
      return;
    }
    values.resize(encoder->slot_count(), 0.0);
    encoder->encode(*context, values, scale, plain);
  }

  inline void encode(std::vector<std::complex<double>> complex_values, double scale, PhantomPlaintext &plain) {
    if (complex_values.size() == 1) {
      encode(complex_values[0], scale, plain);
      return;
    }
    complex_values.resize(encoder->slot_count(), 0.0);
    encoder->encode(*context, pack_complex(complex_values), scale, plain);
  }

  // Value inputs (fill all slots with that value)
  inline void encode(double value, size_t chain_index, double scale, PhantomPlaintext &plain) {
    std::vector<double> values(encoder->slot_count(), value);
    encoder->encode(*context, values, scale, plain, chain_index);
  }

  inline void encode(double value, double scale, PhantomPlaintext &plain) {
    std::vector<double> values(encoder->slot_count(), value);
    encoder->encode(*context, values, scale, plain);
  }

  inline void encode(std::complex<double> complex_value, double scale, PhantomPlaintext &plain) {
    std::vector<cuDoubleComplex> complex_values(encoder->slot_count(), make_cuDoubleComplex(complex_value.real(), complex_value.imag()));
    encoder->encode(*context, complex_values, scale, plain);
  }

  template <typename T, typename = std::enable_if_t<std::is_same<std::remove_cv_t<T>, double>::value || std::is_same<std::remove_cv_t<T>, cuDoubleComplex>::value>>
  inline void decode(PhantomPlaintext &plain, std::vector<T> &values) {
    encoder->decode(*context, plain, values);
  }

  inline void decode(PhantomPlaintext &plain, std::vector<std::complex<double>> &values) {
    auto wrapper = std::vector<cuDoubleComplex>(values.size());
    decode(plain, wrapper);
    for (size_t i = 0; i < values.size(); i++) {
      values[i] = std::complex<double>(cuCreal(wrapper[i]), cuCimag(wrapper[i]));
    }
  }
};

class Encryptor {
 private:
  std::shared_ptr<PhantomContext> context;
  std::shared_ptr<PhantomPublicKey> encryptor;

 public:
  Encryptor() = default;

  Encryptor(std::shared_ptr<PhantomContext> context, std::shared_ptr<PhantomPublicKey> encryptor) {
    this->context = context;
    this->encryptor = encryptor;
  }

  inline void encrypt(PhantomPlaintext &plain, PhantomCiphertext &ct) {
    encryptor->encrypt_asymmetric(*context, plain, ct);
  }
};

// Symmetric Encryptor
// class Encryptor {
//  private:
//   std::shared_ptr<PhantomContext> context;
//   std::shared_ptr<PhantomSecretKey> encryptor;

//  public:
//   Encryptor() = default;

//   Encryptor(std::shared_ptr<PhantomContext> context, std::shared_ptr<PhantomSecretKey> encryptor) {
//     this->context = context;
//     this->encryptor = encryptor;
//   }

//   inline void encrypt(PhantomPlaintext &plain, PhantomCiphertext &ct) {
//     encryptor->encrypt_symmetric(*context, plain, ct);
//   }
// };

class Evaluator {
 private:
  std::shared_ptr<PhantomContext> context;
  std::shared_ptr<PhantomCKKSEncoder> encoder;

 public:
  Evaluator() = default;
  Evaluator(std::shared_ptr<PhantomContext> context, std::shared_ptr<PhantomCKKSEncoder> encoder) {
    this->context = context;
    this->encoder = encoder;
  }

  // Mod switch
  inline void mod_switch_to_next_inplace(PhantomCiphertext &ct) {
    phantom::mod_switch_to_next_inplace(*context, ct);
  }

  inline void mod_switch_to_inplace(PhantomCiphertext &ct, size_t chain_index) {
    phantom::mod_switch_to_inplace(*context, ct, chain_index);
  }

  inline void mod_switch_to_next_inplace(PhantomPlaintext &pt) {
    phantom::mod_switch_to_next_inplace(*context, pt);
  }

  inline void mod_switch_to_inplace(PhantomPlaintext &pt, size_t chain_index) {
    phantom::mod_switch_to_inplace(*context, pt, chain_index);
  }

  inline void rescale_to_next_inplace(PhantomCiphertext &ct) {
    phantom::rescale_to_next_inplace(*context, ct);
  }

  // Relinearization
  inline void relinearize_inplace(PhantomCiphertext &ct, const PhantomRelinKey &relin_keys) {
    phantom::relinearize_inplace(*context, ct, relin_keys);
  }

  // Multiplication
  inline void square(PhantomCiphertext &ct, PhantomCiphertext &dest) {
    multiply(ct, ct, dest);
  }

  inline void square_inplace(PhantomCiphertext &ct) {
    multiply_inplace(ct, ct);
  }

  inline void multiply(PhantomCiphertext &ct1, const PhantomCiphertext &ct2, PhantomCiphertext &dest) {
    if (&ct2 == &dest) {
      multiply_inplace(dest, ct1);
    } else {
      dest = ct1;
      multiply_inplace(dest, ct2);
    }
  }

  inline void multiply_inplace(PhantomCiphertext &ct1, const PhantomCiphertext &ct2) {
    phantom::multiply_inplace(*context, ct1, ct2);
  }

  inline void multiply_plain(PhantomCiphertext &ct, PhantomPlaintext &plain, PhantomCiphertext &dest) {
    dest = phantom::multiply_plain(*context, ct, plain);
  }

  inline void multiply_plain_inplace(PhantomCiphertext &ct, PhantomPlaintext &plain) {
    phantom::multiply_plain_inplace(*context, ct, plain);
  }

  // Addition
  inline void add_plain(const PhantomCiphertext &ct, PhantomPlaintext &plain, PhantomCiphertext &dest) {
    dest = phantom::add_plain(*context, ct, plain);
  }

  inline void add_plain_inplace(PhantomCiphertext &ct, PhantomPlaintext &plain) {
    phantom::add_plain_inplace(*context, ct, plain);
  }

  inline void add(PhantomCiphertext &ct1, const PhantomCiphertext &ct2, PhantomCiphertext &dest) {
    dest = phantom::add(*context, ct1, ct2);
  }

  inline void add_inplace(PhantomCiphertext &ct1, const PhantomCiphertext &ct2) {
    phantom::add_inplace(*context, ct1, ct2);
  }

  inline void add_many(std::vector<PhantomCiphertext> &cts, PhantomCiphertext &dest) {
    size_t size = cts.size();
    if (size < 2) throw std::invalid_argument("add_many requires at least 2 ciphertexts");

    add(cts[0], cts[1], dest);
    for (size_t i = 2; i < size; i++) {
      add_inplace(dest, cts[i]);
    }
  }

  // Subtraction
  inline void sub_plain(PhantomCiphertext &ct, PhantomPlaintext &plain, PhantomCiphertext &dest) {
    dest = ct;
    sub_plain_inplace(dest, plain);
  }

  inline void sub_plain_inplace(PhantomCiphertext &ct, PhantomPlaintext &plain) {
    phantom::sub_plain_inplace(*context, ct, plain);
  }

  inline void sub(PhantomCiphertext &ct1, const PhantomCiphertext &ct2, PhantomCiphertext &dest) {
    if (&ct2 == &dest) {
      sub_inplace(dest, ct1);
      negate_inplace(dest);
    } else {
      dest = ct1;
      sub_inplace(dest, ct2);
    }
  }

  inline void sub_inplace(PhantomCiphertext &ct1, const PhantomCiphertext &ct2) {
    phantom::sub_inplace(*context, ct1, ct2);
  }

  // Rotation
  inline void rotate_vector(PhantomCiphertext &ct, int steps, PhantomGaloisKey &galois_keys, PhantomCiphertext &dest) {
    dest = phantom::rotate(*context, ct, steps, galois_keys);
    cudaStreamSynchronize(ct.data_ptr().get_stream());  // this is currently required, rotation is unstable
  }

  inline void rotate_vector_inplace(PhantomCiphertext &ct, int steps, PhantomGaloisKey &galois_keys) {
    phantom::rotate_inplace(*context, ct, steps, galois_keys);
    cudaStreamSynchronize(ct.data_ptr().get_stream());  // this is currently required, rotation is unstable
  }

  // Negation
  inline void negate(PhantomCiphertext &ct, PhantomCiphertext &dest) {
    dest = ct;
    negate_inplace(dest);
  }

  inline void negate_inplace(PhantomCiphertext &ct) {
    phantom::negate_inplace(*context, ct);
  }

  // Galois
  inline void apply_galois(PhantomCiphertext &ct, uint32_t elt, PhantomGaloisKey &galois_keys, PhantomCiphertext &dest) {
    dest = phantom::apply_galois(*context, ct, elt, galois_keys);
  }

  inline void apply_galois_inplace(PhantomCiphertext &ct, int step, PhantomGaloisKey &galois_keys) {
    auto elt = context->key_galois_tool_->get_elts_from_steps({step})[0];
    phantom::apply_galois_inplace(*context, ct, elt, galois_keys);
  }

  // Complex Conjugate
  inline void complex_conjugate(PhantomCiphertext &ct, const PhantomGaloisKey &galois_keys, PhantomCiphertext &dest) {
    dest = ct;
    complex_conjugate_inplace(dest, galois_keys);
  }

  inline void complex_conjugate_inplace(PhantomCiphertext &ct, const PhantomGaloisKey &galois_keys) {
    phantom::rotate_inplace(*context, ct, 0, galois_keys);
  }

  // Matrix Multiplication
  inline void transform_from_ntt(const PhantomCiphertext &ct, PhantomCiphertext &dest) {
    dest = ct;
    transform_from_ntt_inplace(dest);
  }

  inline void transform_from_ntt_inplace(PhantomCiphertext &ct) {
    auto rns_coeff_count = ct.poly_modulus_degree() * ct.coeff_modulus_size();

    const auto stream = ct.data_ptr().get_stream();

    for (size_t i = 0; i < ct.size(); i++) {
      uint64_t *ci = ct.data() + i * rns_coeff_count;
      nwt_2d_radix8_backward_inplace(ci, context->gpu_rns_tables(), ct.coeff_modulus_size(), 0, stream);
    }

    ct.set_ntt_form(false);
    // cudaStreamSynchronize(stream);
  }

  inline void transform_to_ntt(const PhantomCiphertext &ct, PhantomCiphertext &dest) {
    dest = ct;
    transform_to_ntt_inplace(dest);
  }

  inline void transform_to_ntt_inplace(PhantomCiphertext &ct) {
    auto rns_coeff_count = ct.poly_modulus_degree() * ct.coeff_modulus_size();
    const auto stream = ct.data_ptr().get_stream();

    for (size_t i = 0; i < ct.size(); i++) {
      uint64_t *ci = ct.data() + i * rns_coeff_count;
      nwt_2d_radix8_forward_inplace(ci, context->gpu_rns_tables(), ct.coeff_modulus_size(), 0, stream);
    }

    ct.set_ntt_form(true);
    // cudaStreamSynchronize(stream);
  }

  // Bootstrapping
  inline void multiply_const(const PhantomCiphertext &ct, double value, PhantomCiphertext &dest) {
    dest = ct;
    multiply_const_inplace(dest, value);
  }

  inline void multiply_const_inplace(PhantomCiphertext &ct, double value) {
    PhantomPlaintext const_plain;

    std::vector<double> values(encoder->slot_count(), value);
    encoder->encode(*context, values, ct.scale(), const_plain);
    mod_switch_to_inplace(const_plain, ct.chain_index());
    multiply_plain_inplace(ct, const_plain);
  }

  inline void add_const(PhantomCiphertext &ct, double value, PhantomCiphertext &dest) {
    dest = ct;
    add_const_inplace(dest, value);
  }

  inline void add_const_inplace(PhantomCiphertext &ct, double value) {
    PhantomPlaintext const_plain;

    std::vector<double> values(encoder->slot_count(), value);
    encoder->encode(*context, values, ct.scale(), const_plain);
    mod_switch_to_inplace(const_plain, ct.chain_index());
    add_plain_inplace(ct, const_plain);
  }

  inline void add_reduced_error(const PhantomCiphertext &ct1, const PhantomCiphertext &ct2, PhantomCiphertext &dest) {
    if (&ct2 == &dest) {
      add_inplace_reduced_error(dest, ct1);
    } else {
      dest = ct1;
      add_inplace_reduced_error(dest, ct2);
    }
  }

  void add_inplace_reduced_error(PhantomCiphertext &ct1, const PhantomCiphertext &ct2);

  inline void sub_reduced_error(const PhantomCiphertext &ct1, const PhantomCiphertext &ct2, PhantomCiphertext &dest) {
    dest = ct1;
    sub_inplace_reduced_error(dest, ct2);
  }

  void sub_inplace_reduced_error(PhantomCiphertext &ct1, const PhantomCiphertext &ct2);

  inline void multiply_reduced_error(const PhantomCiphertext &ct1, const PhantomCiphertext &ct2, const PhantomRelinKey &relin_keys, PhantomCiphertext &dest) {
    if (&ct2 == &dest) {
      multiply_inplace_reduced_error(dest, ct1, relin_keys);
    } else {
      dest = ct1;
      multiply_inplace_reduced_error(dest, ct2, relin_keys);
    }
  }

  void multiply_inplace_reduced_error(PhantomCiphertext &ct1, const PhantomCiphertext &ct2, const PhantomRelinKey &relin_keys);

  inline void double_inplace(PhantomCiphertext &ct) const {
    phantom::add_inplace(*context, ct, ct);
  }

  template <typename T, typename = std::enable_if_t<std::is_same<std::remove_cv_t<T>, double>::value || std::is_same<std::remove_cv_t<T>, std::complex<double>>::value>>
  inline void multiply_vector_reduced_error(PhantomCiphertext &ct, std::vector<T> &values, PhantomCiphertext &dest) {
    dest = ct;
    multiply_vector_inplace_reduced_error(dest, values);
  }

  inline void multiply_vector_inplace_reduced_error(PhantomCiphertext &ct, std::vector<double> &values) {
    PhantomPlaintext plain;

    values.resize(encoder->slot_count(), 0.0);
    encoder->encode(*context, values, ct.scale(), plain);
    mod_switch_to_inplace(plain, ct.chain_index());
    multiply_plain_inplace(ct, plain);
  }

  inline void multiply_vector_inplace_reduced_error(PhantomCiphertext &ct, std::vector<std::complex<double>> &values) {
    PhantomPlaintext plain;

    values.resize(encoder->slot_count(), 0.0);
    encoder->encode(*context, pack_complex(values), ct.scale(), plain);
    mod_switch_to_inplace(plain, ct.chain_index());
    multiply_plain_inplace(ct, plain);
  }
};

class Decryptor {
 private:
  std::shared_ptr<PhantomContext> context;
  std::shared_ptr<PhantomSecretKey> decryptor;

 public:
  Decryptor() = default;
  Decryptor(std::shared_ptr<PhantomContext> context, std::shared_ptr<PhantomSecretKey> decryptor) {
    this->context = context;
    this->decryptor = decryptor;
  }

  inline void decrypt(PhantomCiphertext &ct, PhantomPlaintext &plain) {
    decryptor->decrypt(*context, ct, plain);
  }

  inline void create_galois_keys_from_elts(std::vector<uint32_t> &elts, PhantomGaloisKey &galois_keys) {
    const auto &s = cudaStreamPerThread;

    int log_n = phantom::arith::get_power_of_two(context->poly_degree_);
    bool is_bfv = (context->first_context_data().parms().scheme() == phantom::scheme_type::bfv);
    
    context->key_galois_tool_.reset();
    context->key_galois_tool_ = std::make_unique<phantom::util::PhantomGaloisTool>(elts, log_n, s, is_bfv);

    galois_keys = decryptor->create_galois_keys(*context);
  }

  inline void create_galois_keys_from_steps(std::vector<int> &steps, PhantomGaloisKey &galois_keys) {
    auto elts = context->key_galois_tool_->get_elts_from_steps(steps);
    create_galois_keys_from_elts(elts, galois_keys);
  }
};

class CKKSEvaluator {
 private:
  // Sign function coefficients
  std::vector<double> F4_COEFFS = {0, 315, 0, -420, 0, 378, 0, -180, 0, 35};
  int F4_SCALE = (1 << 7);
  std::vector<double> G4_COEFFS = {0, 5850, 0, -34974, 0, 97015, 0, -113492, 0, 46623};
  int G4_SCALE = (1 << 10);

  // Helper functions
  uint64_t get_modulus(PhantomCiphertext &x, int k);

  PhantomCiphertext init_guess(PhantomCiphertext x);
  PhantomCiphertext eval_line(PhantomCiphertext x, PhantomPlaintext m, PhantomPlaintext c);

  // Evaluation functions
  PhantomCiphertext newton_iter(PhantomCiphertext x, PhantomCiphertext res, int iter);
  std::pair<PhantomCiphertext, PhantomCiphertext> goldschmidt_iter(PhantomCiphertext v, PhantomCiphertext y, int d = 1);
  void eval_odd_deg9_poly(std::vector<double> &a, PhantomCiphertext &x, PhantomCiphertext &dest);

 public:
  // Memory managed outside of the evaluator
  std::shared_ptr<PhantomContext> context;
  std::shared_ptr<PhantomRelinKey> relin_keys;
  std::shared_ptr<PhantomGaloisKey> galois_keys;
  std::vector<std::uint32_t> galois_elts;

  // Component classes
  Encoder encoder;
  Encryptor encryptor;
  Evaluator evaluator;
  Decryptor decryptor;

  size_t degree;
  double scale;
  size_t slot_count;

  CKKSEvaluator(
    std::shared_ptr<PhantomContext> context, 
    std::shared_ptr<PhantomPublicKey> encryptor, 
    std::shared_ptr<PhantomSecretKey> decryptor,
    std::shared_ptr<PhantomCKKSEncoder> encoder, 
    std::shared_ptr<PhantomRelinKey> relin_keys, 
    std::shared_ptr<PhantomGaloisKey> galois_keys,
    double scale, 
    std::vector<uint32_t> galois_elts = {}
  ) {
    this->context = context;
    this->relin_keys = relin_keys;
    this->galois_keys = galois_keys;
    this->galois_elts = galois_elts;

    this->scale = scale;
    this->slot_count = encoder->slot_count();
    this->degree = this->slot_count * 2;

    // Instantiate the component classes
    Encoder ckks_encoder(context, encoder);
    this->encoder = ckks_encoder;

    Encryptor ckks_encryptor(context, encryptor);
    this->encryptor = ckks_encryptor;

    Evaluator ckks_evaluator(context, encoder);
    this->evaluator = ckks_evaluator;

    Decryptor ckks_decryptor(context, decryptor);
    this->decryptor = ckks_decryptor;
  }

  // Helper functions
  std::vector<double> init_vec_with_value(double value);
  PhantomPlaintext init_plain_power_of_x(size_t exponent);

  void re_encrypt(PhantomCiphertext &ct);
  void print_decrypted_ct(PhantomCiphertext &ct, int num);
  void print_decoded_pt(PhantomPlaintext &pt, int num);

  // Evaluation functions
  PhantomCiphertext sgn_eval(PhantomCiphertext x, int d_g, int d_f, double sgn_factor = 0.5);
  PhantomCiphertext invert_sqrt(PhantomCiphertext x, int d_newt = 20, int d_gold = 1);
  PhantomCiphertext exp(PhantomCiphertext x);
  PhantomCiphertext inverse(PhantomCiphertext x, int iter = 4);

  // Metrics calcuation functions
  double calculate_MAE(std::vector<double> &y_true, PhantomCiphertext &ct, int N);
};

}  // namespace nexus
