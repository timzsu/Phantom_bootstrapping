#pragma once

#include <NTL/RR.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "Point.cuh"

namespace minicomp {
NTL::ZZ mod(NTL::ZZ a, NTL::ZZ q);
long pmod(long a, long b);
long pow2(long n);
size_t ceil_divide(size_t x, size_t y);
size_t ceil_to_int(double x);
int floor_to_int(double x);
long log2_long(long n);
long num_one(long n);
NTL::RR GetApproxError(int d, NTL::RR t, bool is_first_function);
NTL::RR GetInvApproxError(int d, NTL::RR t, NTL::RR *X, NTL::RR *Y, long num);
NTL::RR GetInvApproxError(int d, NTL::RR t, std::vector<NTL::RR> &X, std::vector<NTL::RR> &Y, long num);
NTL::RR real(NTL::RR p);
int dep(std::list<int> &lt);
int dep(int deg);
int mult(int deg);
int enough_mult(int a);
int enough_dep(int a);
NTL::RR exptoreal(NTL::RR x);
NTL::RR realtoexp(NTL::RR x);
NTL::RR sgn(NTL::RR x);
NTL::RR comp(NTL::RR a, NTL::RR b);
NTL::RR fracpart(NTL::RR x);
NTL::RR ReLU(NTL::RR x);
double comp(double a, double b);
double ReLU(double x);
bool yabscompare(Point a, Point b);
bool ycompare(Point a, Point b);
bool xcompare(Point a, Point b);
bool isin(NTL::RR x, long K, NTL::RR width);
NTL::RR eval(long deg, NTL::RR *coeff, NTL::RR val, int type, NTL::RR scale);
NTL::RR eval(long deg, std::vector<NTL::RR> &coeff, NTL::RR val, int type, NTL::RR scale);
void showgraph(std::ofstream &out, NTL::RR *coeff, long deg, NTL::RR start, NTL::RR end, NTL::RR sc, int type, NTL::RR scale);
void showgraph(std::ofstream &out, std::vector<NTL::RR> coeff, long deg, NTL::RR start, NTL::RR end, NTL::RR sc, int type, NTL::RR scale);
bool oddevennextcombi(long *arr, long arrlen, long len);
NTL::RR expmaxerr(long deg, NTL::RR expx);
NTL::RR invexpmaxerr(long deg, NTL::RR expy, NTL::RR *X, NTL::RR *Y, long num);
NTL::RR invexpmaxerr(long deg, NTL::RR expy, std::vector<NTL::RR> &X, std::vector<NTL::RR> &Y, long num);
NTL::RR getmaxerr(NTL::RR (*func)(NTL::RR), std::vector<NTL::RR> &coeff, long deg, NTL::RR start, NTL::RR end, int type, NTL::RR scale, long prec, long num, bool is_opt_sampling);
NTL::RR find_extreme(NTL::RR (*func)(NTL::RR), Point *&ext, int &ext_count, std::vector<NTL::RR> coeff, long deg, NTL::RR start, NTL::RR end, long prec, NTL::RR scan, int type, NTL::RR scale, bool is_opt_sampling);
}  // namespace minicomp
