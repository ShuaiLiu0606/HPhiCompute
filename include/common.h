#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <complex>
#include <vector>
#include <tuple>
#include <iomanip>
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>
#include <cstdint>
#include <map>
#include <array>
#include <chrono>
#include <cassert>
#include <omp.h>

#define __MKL
#ifdef __MKL
#define MKL_Complex16 std::complex<double>
#include <mkl.h>
#include <mkl_types.h>
#include <mkl_blas.h>
#include <mkl_lapack.h>
#define EIGEN_USE_MKL_ALL
#endif

typedef std::chrono::high_resolution_clock::time_point TimeVar;
#define duration(a) std::chrono::duration_cast<std::chrono::seconds>(a).count()
#define durationms(a) std::chrono::duration_cast<std::chrono::milliseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

using std::abs;
using std::array;
using std::cerr;
using std::cin;
using std::complex;
using std::conj;
using std::cos;
using std::cout;
using std::endl;
using std::fill;
using std::fixed;
using std::flush;
using std::get;
using std::ifstream;
using std::ios;
using std::istringstream;
using std::left;
using std::log;
using std::log2;
using std::map;
using std::max;
using std::min;
using std::ofstream;
using std::ostringstream;
using std::pair;
using std::pow;
using std::right;
using std::setfill;
using std::setprecision;
using std::setw;
using std::sin;
using std::string;
using std::stringstream;
using std::to_string;
using std::tuple;
using std::vector;

/*
using std::cin; using std::cout; using std::cerr; using std::endl; using std::ifstream; using std::ofstream;
using std::complex; using std::vector; using std::array; using std::string; using std::tuple; using std::map; using std::pair;
using std::abs; using std::cos; using std::conj; using std::log; using std::log2; using std::max; using std::min; using std::pow; using std::sin;
using std::istringstream; using std::stringstream; using std::ostringstream;
using std::ios; using std::fill; using std::fixed; using std::flush; using std::setfill; using std::setprecision; using std::setw; using std::left; using std::right; using std::to_string;
using std::get; */

#define PI 3.14159265358979323846264338327950288419
const complex<double> I = complex<double>(0.0, 1.0);

using my_int = MKL_INT;
using my_cplx = complex<double>;
using my_vec = vector<double>;
using my_cplx_vec = vector<complex<double>>;
using my_mat = vector<vector<double>>;
using my_cplx_mat = vector<vector<complex<double>>>;

// -----------------------------------------------------------------------------
// 类型定义
//
// Term: 一个算符项 (key, 系数)
// key 用字符串表示算符 c†_a c_b
//
// 例如:
// key = "21" 表示算符  c†_{2} c_{1}
//
// Operator: 一个算符由多个项组成
// -----------------------------------------------------------------------------

using Term = pair<string, complex<double>>;
using Operator = vector<Term>;

//************************************************************************

void GetAllFilesWithPrefix(const string &folder_path, vector<vector<string>> &All_files, my_vec &extractedValues1, my_mat &extractedValues2);

Operator multiply_terms(const Operator &terms, const Operator &factor);

vector<pair<string, int>> compute_repeated_operator(const vector<pair<string, int>> &operator_, int steps);

Operator build_Sz(double Spin);
Operator build_Sx(double Spin);
Operator build_Sy(double Spin);
Operator build_Splus(double Spin);
Operator build_Sminus(double Spin);
void ComS_Corr(int Sdim, const Operator &result, const my_mat &OneBody, int site, double &realPart, double &imagPart);
void ComSS_Corr(int Nsites, int Sdim, const Operator &result, const my_mat &TwoBody, int sitei, int sitej, double &realPart, double &imagPart);
void ComSSS_Corr(int Nsites, int Sdim, const Operator &result, const my_mat &ThreeBody, int sitei, int sitej, int sitek, double &realPart, double &imagPart);
void ComSSSS_Corr(int Nsites, int Sdim, const Operator &result, const my_mat &FourBody, int sitei, int sitej, int sitem, int siten, double &realPart, double &imagPart);
void ComSSSSSS_Corr(int Nsites, int Sdim, const Operator &result, const my_mat &SixBody, int plaq_index, double &realPart, double &imagPart);

//***********************随机向量归一化******************************
void rd_wf(const int &dim, double *wf);

void rd_wf(const int &dim, my_cplx *wf);

void rd_wf_product(const int &dim, double *wf);

void uni_wf(const int &dim, double *wf);

//**********************bitwise operations*************************
int Bits_Count1(const int &b);

int Bits_CycleRight(const int &s, const int &n);

int Bits_CycleLeft(const int &s, const int &n);

int Bits_Reflect(const int &s, const int &n);

int Bits_Invert(const int &s, const int &n);

int Bits_HammingDis(const int &a, const int &b);

void Bits_Decomposition(const uint32_t &s, const uint32_t &n, const uint32_t &size_A, uint32_t *sites_A, uint32_t &a, uint32_t &b);

void Bits_Print(const int &a, const int &n);

int nchoosek(const int &n, const int &_k);

void Debug_Check_Conj_dense(double *mat, int dim);

//************************读取文件等号后面的数据****************************

void GetPara(const string &fname, const string &string_match, int &para);

void GetPara(const string &fname, const string &string_match, double &para);

void GetPara(const string &fname, const string &string_match, char &para);

//**********************向量,矩阵的写入，读取，输出*************************

void Vec_Read(const string &fname, vector<double> &data);

void Vec_Read(const string &fname, vector<complex<double>> &data);

void Vec_Print(const double *data, const int dsize);

void Vec_Print(const complex<double> *data, const int dsize);

void Mat_Read(const string &filename, vector<vector<int>> &matrix);

void Mat_Read(const string &filename, vector<vector<double>> &matrix);

void Mat_Read(const string &filename, vector<vector<complex<double>>> &matrix);

void Mat_Print(const complex<double> *data, const int rows, const int cols);

void Mat_Print(const double *data, const int rows, const int cols);

void Mat_Pointer(const vector<vector<double>> &matrix, double *flat_array); // 矩阵转化为一维指针数组

void Mat_Pointer(const vector<vector<complex<double>>> &matrix, complex<double> *flat_array);

//***************************二进制向量写入，读取**************************
void Vec_Fwrite(const string &fname, const double *data, const int dsize);

void Vec_Fread(const string &fname, double *data, const int dsize);

void Vec_Fwrite(const string &fname, const complex<double> *data, const int dsize);

void Vec_Fread(const string &fname, complex<double> *data, const int dsize);

//****************************向量，矩阵运算*******************************
double Vec_Dot(const int &dim, double *x, const int incx, double *y, const int incy);

my_cplx Vec_Dot(const int &dim, my_cplx *x, const int incx, my_cplx *y, const int incy);

void Out_Dot(int row, const double *vec1, int col, const double *vec2, double *result);

void Out_Dot(int row, const complex<double> *vec1, int col, const complex<double> *vec2, complex<double> *result);

void Mat_Vec_Dot(const double *matrix, const double *vec, int rows, int cols, double *result);

void Mat_Vec_Dot(const my_cplx *matrix, const my_cplx *vec, int rows, int cols, my_cplx *result);

void Mat_Dot(CBLAS_LAYOUT matrix_layout, CBLAS_TRANSPOSE A_transpose, CBLAS_TRANSPOSE B_transpose, int M, int N, int K, my_cplx *A, int lda, my_cplx *B, int ldb, my_cplx *C, int ldc);

void Mat_Dot(CBLAS_LAYOUT matrix_layout, CBLAS_TRANSPOSE A_transpose, CBLAS_TRANSPOSE B_transpose, int M, int N, int K, double *A, int lda, double *B, int ldb, double *C, int ldc);

void Vec_Cross(const double *vec1, const double *vec2, double *result);

void Vec_Cross(const my_cplx *vec1, const my_cplx *vec2, my_cplx *result);

double *Kronecker(const double *matA, const double *matB, int A_rows, int A_cols, int B_rows, int B_cols);

my_cplx *Kronecker(const my_cplx *matA, const my_cplx *matB, int A_rows, int A_cols, int B_rows, int B_cols);
//************************************有限温cTPQ*********************************************

my_vec cla_mean(const my_mat &data, int axis);

my_vec cal_standard_deviation(const my_mat &data, int axis, bool ddof);

my_mat cla_mean(const vector<my_mat> &data, int axis);

my_mat cal_standard_deviation(const vector<my_mat> &data, int axis, bool ddof);

//*******************************************others*******************************************
double VN_Entropy(const int &dim, double *p);

pair<double, double> Rotate_Coord(double x, double y, double angle);

void CopyFile(const string &sourcePath, const string &destinationPath);

//******************************Matrix evd (use mkl lapacke function)*************************i
void Mat_Svd(int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, double *a, lapack_int lda, double *s, double *u, lapack_int ldu, double *vt, lapack_int ldvt);

void Mat_Svd(int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, my_cplx *a, lapack_int lda, double *s, my_cplx *u, lapack_int ldu, my_cplx *vt, lapack_int ldvt);

void DenseMatrixEigenSolver(my_int matrix_layout, char jobz, char uplo, lapack_int n, double *a, lapack_int lda, double *w);

void DenseMatrixEigenSolver(my_int matrix_layout, char jobz, char uplo, lapack_int n, my_cplx *a, lapack_int lda, double *w);

void Lanczos(const int &dim, double *vals, int *cols, int *PointerBE, const int &M, char job_vec, double *alpha, double *l_vecs, double *l_eigvecs);

void Lanczos(const int &dim, my_cplx *vals, int *cols, int *PointerBE, const int &M, char job_vec, double *alpha, my_cplx *l_vecs, double *l_eigvecs);

// y = a*x
void Vec_ax(const int &dim, const my_cplx &a, my_cplx *x, my_cplx *y);

void Vec_ax(const int &dim, const double &a, double *x, double *y);

// a = real(x)
void Vec_real(const int &dim, my_cplx *x, double *a);

// a = imag(x)
void Vec_imag(const int &dim, my_cplx *x, double *a);

// a = x + i*y, x and y are real vectors
void Vec_xpiy(const int &dim, double *x, double *y, my_cplx *a);

#endif