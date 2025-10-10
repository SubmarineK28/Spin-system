

#include "Matrix.h"

using namespace std::chrono;

    //kron_mult_recursive_old(term, N, vec, result_2, 0, 0, 0, dim);
//    printVector(result_2);
//int main() {
//
//    setlocale(LC_ALL, "Russian");
//    size_t N = 3;
//    int dim = 1 << N;
//
// // генерация матрицы 
//    vector<Matrix2x2> term;
//    random_device rd;
//    mt19937 gen(rd()); 
//    normal_distribution<double> dist(0.0, 1.0);
//
//    for (int i = 0; i < N; i++) {
//        Matrix2x2 M(4);
//        for (int j = 0; j < 4; j++) {
//            M[j] = Complex(dist(gen), 0.0);
//        }
//        term.push_back(M);
//    }
//    
//
//    auto vec = random_normal_vector(dim, 0.0, 5.0);
//    std::vector<Complex> tmp(dim, 0.0);
//
//    //printVector(vec);
//    //Matrix2x2 result_1(dim, 0.0);
//    Matrix2x2 result_2(dim, 0.0);
//
//    //vector<Complex> result(dim, 0);
//
//    auto start_total = high_resolution_clock::now();
//    //kron_mult_recursive(term, N, vec, result_1, 0, 0, 0, dim);
//    kron_mult_recursive_double(term, 0, vec.data(), tmp.data(), dim, dim);
//    auto end_total = high_resolution_clock::now();
//    auto duration_1 = duration_cast<milliseconds>(end_total - start_total);
//    cout << "Execution total time: " << duration_1.count() << " s" << endl;
//    //printVector(result_1);
//    //cout << "\n";
//    printVector(tmp);
//
//
//
//    return 0;
//}