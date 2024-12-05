#pragma once

#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

#include "Point.cuh"

NTL::RR fracpart(NTL::RR x);
NTL::RR sqrtfracpart(NTL::RR x, NTL::RR a);
NTL::RR fraccos(NTL::RR x, long scale);
NTL::RR arcsin(NTL::RR x);
bool yabscompare(Point a, Point b);
bool ycompare(Point a, Point b);
bool xcompare(Point a, Point b);
bool isin(NTL::RR x, long K, NTL::RR width);
NTL::RR chebeval(long deg, NTL::RR *coeff, NTL::RR val);
void showgraph(std::ofstream &out, NTL::RR *coeff, long deg, long K, NTL::RR sc);
bool oddevennextcombi(long *arr, long arrlen, long len);
void oddbabycount(long &k, long &m, long deg);
void babycount(long &k, long &m, long deg);

void add(std::complex<double> *&rtn, std::complex<double> *&vec1, std::complex<double> *&vec2, long n);
void addinplace(std::complex<double> *&vec, std::complex<double> *&addvec, long n);
void subt(std::complex<double> *&rtn, std::complex<double> *&vec1, std::complex<double> *&vec2, long n);
void subtinplace(std::complex<double> *&vec, std::complex<double> *&subtvec, long n);
void mul(std::complex<double> *&rtn, std::complex<double> *&vec1, std::complex<double> *&vec2, long n);
void mulinplace(std::complex<double> *&vec, std::complex<double> *&mulvec, long n);
void constmul(std::complex<double> *&rtn, std::complex<double> *&vec, std::complex<double> constant, long n);
void constmulinplace(std::complex<double> *&vec, std::complex<double> constant, long n);
void text_to_array(std::ifstream &in, NTL::RR *&array, long n);
int giantstep(int M);
void rotation(int logslot, int Nh, int shiftcount, const std::vector<std::complex<double>> &vec, std::vector<std::complex<double>> &rtnvec);

int max_index(double *array, int length);
