#pragma once
#include <vector>
#include <complex>
#include <iostream>
#include <random>
#include <algorithm>
#include <iomanip>
#include <locale>
#include <chrono>

using namespace std;

// 0: I, 1: X, 2: Y, 3: Z
enum PauliType { I = 0, X = 1, Y = 2, Z = 3 };

using Complex = std::complex<double>;
using Matrix2x2 = std::vector<Complex>;
using Structure = std::vector<std::vector<Matrix2x2>>;

// --- Матрицы Паули ---
Matrix2x2 pauli_I();
Matrix2x2 pauli_X();
Matrix2x2 pauli_Y();
Matrix2x2 pauli_Z();
Matrix2x2 pauli_matrix(int code);

// --- Построение структур ---
std::vector<std::vector<int>> createH0_structure(int N);
std::vector<std::vector<int>> createM0_structure(int N);
std::vector<double>  factor(int N);
Structure build_structure(const std::vector<std::vector<int>>& vec_full, 
	const std::vector<double>& factors);

// --- Печать ---
void print_matrix(const Matrix2x2& M);
void print_structure_codes(const std::vector<std::vector<int>>& vec_full);
void print_structure_matrices(const Structure& S);
void print_structure_term(const vector< Matrix2x2>& term);
void printVector(const std::vector<Complex>& v);

// --- Сумма ---
Matrix2x2 add_vectors(const std::vector<Complex>& a, const std::vector<Complex>& b);

// --- Произведение Кронекера --- 

void kron_mult_recursive(   ///////////// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! AAAAAAAAAAAAAAAAAAAAAAAA
    const std::vector<Matrix2x2>& factors,
    size_t N,
    const Matrix2x2& vec_in,
    Matrix2x2& vec_out,
    size_t depth,
    size_t offset_in,
    size_t offset_out,
    size_t stride
);

void kron_mult_recursive_old(   ///////////// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! AAAAAAAAAAAAAAAAAAAAAAAA
    const std::vector<Matrix2x2>& factors,
    size_t N,
    const Matrix2x2& vec_in,
    Matrix2x2& vec_out,
    size_t depth,
    size_t offset_in,
    size_t offset_out,
    size_t stride
);

void kron_mult_recursive_double(
    const std::vector<Matrix2x2>& factors,
    size_t depth,
    const Complex* src,
    Complex* dst,
    size_t stride,
    size_t dim
);

/// --- скалярное умножение ---  
Complex dot_product(const std::vector<Complex>& a, const std::vector<Complex>& b);

vector<Complex> random_normal_vector(size_t dim, double mean, double stddev);
vector<Complex> generateRademacherComplexVector(size_t dimension);
vector<Complex> Hutchinson_vector(size_t n);

void printVector_double(const std::vector<double>& v);
