                       // H0 - [ (A0, A1, ...., An), (A0, A1, ...., An), ... ]
                       // term = (A0, A1, ...., An)
                       // M = An
#include <iomanip>
#include "Matrix.h"  // подключаем функции из Matrix.cpp
#include <chrono>

using namespace std;
using namespace std::chrono;


// --- Обёртка ---
static vector<Complex> kron_mult(const std::vector<std::vector<Matrix2x2>>& factors, const vector<Complex>& vec) {
    size_t N = factors[0].size();
    size_t dim = size_t(1) << N;
    vector<Complex> Sum_0(dim, 0.0);
    for (const auto& term : factors) {
        vector<Complex> result(dim, 0.0);
        kron_mult_recursive(term, N, vec, result, 0, 0, 0, dim);
        Sum_0 = add_vectors(Sum_0, result);
    }

    return Sum_0;
}

// --- Печать вектора ---
static void printVector(const vector<Complex>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        for (size_t i = 0; i < v.size(); ++i) {
            cout << setw(8) << fixed << setprecision(2) << v[i] << " ";
        }
        cout << "\n";
    }
}


// --- main ---
int main() {
    auto start = high_resolution_clock::now();

    size_t k = 1; // power
    double M = 50.0; // samples
    double Betta = 0.0;
    size_t N = 20;  // кол-во частиц (кубитов)
    int dim = 1 << N;
    cout << "Building H0 and M0 structures...\n";

    auto H0_struct = createH0_structure(N);
    auto M0_struct = createM0_structure(N);

    auto konst = factor(N);

    auto H0 = build_structure(H0_struct, konst);
    auto M0 = build_structure(M0_struct, konst);

    //cout << "\nH0 terms (codes):\n";
    //print_structure_codes(H0_struct);
    //cout << "\nH0 terms (matrices):\n";
    //print_structure_matrices(H0);
    //
    //cout << "\n num of factors :\n";
    //printVector_double(konst);
    //
    //cout << "\nM0 terms (codes):\n";
    //print_structure_codes(M0_struct);
    //cout << "\nM0 terms (matrices):\n";
    //print_structure_matrices(M0);

    // Выборка из M
    Complex Sp_H0 = 0.0;
    Complex Sp_M0 = 0.0;
    Complex Sp = 0.0;

    for (size_t i = 0; i < static_cast<size_t>(M); i++) {
        auto vec = random_normal_vector(dim, 0.0, 1.0); // генерация случайного вектора
        //printVector(vec); // print 

        Matrix2x2 Sum_H0;// Возведение в степень k
        Matrix2x2 Sum_M0;

        for (size_t j = 0; j < k; j++) {
            if (j == 0 && Betta == 0) {
                Sum_H0 = kron_mult(H0, vec);
                //printVector(Sum_H0); // print 
            }
            else if (j == 0) {
                Sum_H0 = kron_mult(H0, vec);
                Sum_M0 = kron_mult(M0, vec);
                printVector(Sum_H0); // print 
                printVector(Sum_M0);
            }
            else if (Betta == 0){
                Sum_H0 = kron_mult(H0, Sum_H0); // [factor*( A0 ? ... ? An) + ... + factor*( A0 ? ... ? An)] * vec = vec_1 !!!
            }                                   // Sum_M0 = vec_k = [factor*( A0 ? ... ? An) + ... + factor*( A0 ? ... ? An)] * vec_k-1
            else {
                Sum_H0 = kron_mult(H0, Sum_H0); // [factor*( A0 ? ... ? An) + ... + factor*( A0 ? ... ? An)] * vec = vec_1 !!!
                Sum_M0 = kron_mult(M0, Sum_M0); // ....
            }
        } 
        if (Betta == 0) {
            auto dot_product_H0 = dot_product(vec, Sum_H0); // Умножеие на транспонированный вектор
            Sp_H0 += dot_product_H0;
        }
        else {
            auto dot_product_H0 = dot_product(vec, Sum_H0); // Умножеие на транспонированный вектор
            auto dot_product_M0 = dot_product(vec, Sum_M0);
            Sp_H0 += dot_product_H0;
            Sp_M0 += dot_product_M0;
        }
    }

    Complex coeff = Complex(1.0, 0.0) / static_cast<double>(M);
    Sp = coeff * Sp_H0 + Betta * coeff * Sp_M0;

    cout << "\n Sp:" << Sp << "\n";

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(end - start);
    cout << "Execution time: " << duration.count() << " s" << endl;
    return 0;

    return 0;
}
