

#include "Matrix.h"

using namespace std::chrono;

//int main() {
//
//    setlocale(LC_ALL, "Russian");
//    size_t N = 2;
//    int dim = 1 << N;
//
//    // генерация матрицы 
//    vector<Matrix2x2> term;
//    random_device rd;
//    mt19937 gen(rd());
//    normal_distribution<double> dist(0.0, 1.0);
//
//    for (int i = 0; i < N; i++) {
//        Matrix2x2 M(4);
//        for (int j = 0; j < 4; j++) {
//            M[j] = Complex(round(dist(gen)), 0.0);
//        }
//        term.push_back(M);
//        print_matrix(M);
//    }
//
//    auto vec = random_normal_vector(dim, 0.0, 1.0);
//    printVector(vec);
//    std::vector<Complex> tmp(dim, 0.0);
//
//    auto start_total = high_resolution_clock::now();
//    //kron_mult_recursive(term, N, vec, result_1, 0, 0, 0, dim);
//    //kron_mult_recursive_old(term, N, vec, result_2, 0, 0, 0, dim);
//    kron_mult_recursive_double(term, 0, vec.data(), tmp.data(), dim, dim);
//    auto end_total = high_resolution_clock::now();
//
//    Complex Sp_H0 = 0.0;
//
//    auto dot_product_H0 = dot_product(vec, tmp); // Умножеие на транспонированный вектор
//    Sp_H0 += dot_product_H0;
//
//    printVector(tmp);
//
//    cout << "\n Sp:" << Sp_H0 << "\n";
//
//    return 0;
//}