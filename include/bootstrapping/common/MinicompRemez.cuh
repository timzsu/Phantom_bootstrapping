#pragma once

#include <NTL/mat_RR.h>

#include <string>

#include "Choosemax.cuh"
#include "MinicompFunc.cuh"
#include "Point.cuh"
#define _USE_MATH_DEFINES

namespace minicomp {
class Remez {
 private:
  // input variables
  long deg, iter, prec;
  int type;
  NTL::RR scale, sc;
  //	NTL::RR sc, *coeff, *temp_ext, maxerr;
  std::vector<NTL::RR> coeff;
  NTL::RR(*func)
  (NTL::RR);
  size_t inter_num;
  std::vector<NTL::RR> inter_start, inter_end;
  //	bool is_first_function;
  bool is_opt_sampling;

  // variables for algorithm
  Point *sam, *ext;
  NTL::RR maxerr;
  long ext_count;
  int *ext_maxsum_index;
  NTL::RR *temp_ext;
  //	string filename;
  NTL::vec_RR v, w, v_0;
  NTL::mat_RR m;

 public:
  Remez(NTL::RR (*_func)(NTL::RR), size_t _inter_num, std::vector<NTL::RR> _inter_start, std::vector<NTL::RR> _inter_end, NTL::RR _sc, long _prec, long _deg, long _iter, int _type, NTL::RR _scale, bool _is_opt_sampling);
  ~Remez();
  void initialize();
  void getcoeffwitherr();
  int getextreme();
  void choosemaxs();
  void printgraph(std::ofstream &out, NTL::RR start, NTL::RR end, NTL::RR sc);
  bool test();
  NTL::RR getmaxerr();
  void getcoeff(NTL::RR _coeff[]);
  void getcoeff(std::vector<NTL::RR> &_coeff);
  void getext_xpos(std::vector<NTL::RR> &_ext_xpos);
};
}  // namespace minicomp
