                       // H0 - [ (A0, A1, ...., An), (A0, A1, ...., An), ... ]
                       // term = (A0, A1, ...., An)
                       // M = An
#include <iomanip>
#include "Matrix.h"  // ���������� ������� �� Matrix.cpp
#include <chrono>

using namespace std;
using namespace std::chrono;


 //--- ������ ---
static vector<Complex> kron_mult(const std::vector<std::vector<Matrix2x2>>& factors, const vector<Complex>& vec) {
    size_t N = factors[0].size();
    size_t dim = size_t(1) << N;
    vector<Complex> Sum_0(dim, 0.0);

    for (const auto& term : factors) {
        vector<Complex> result(dim, 0.0);
        cout << "\n";
        print_structure_term(term);
        cout << "\n";
        kron_mult_recursive_double(term, 0, vec.data(), result.data(), dim, dim);
        Sum_0 = add_vectors(Sum_0, result);
    }

    return Sum_0;
}

 //--- ������ ������� ---
static void printVector(const vector<Complex>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        for (size_t i = 0; i < v.size(); ++i) {
            cout << setw(8) << fixed << setprecision(2) << v[i] << " ";
        }
        cout << "\n";
    }
}


//--- main ---

int main() {
    auto start_total = high_resolution_clock::now();

    time_point<high_resolution_clock> start_one_kron;
    time_point<high_resolution_clock> end_one_kron;

    size_t k = 2; // power
    double M = 200.0; // samples
    double Betta = 0.0;
    size_t N = 2;  // ���-�� ������ (�������)
    int dim = 1 << N;
    cout << "Building H0 and M0 structures...\n";

    auto H0_struct = createH0_structure(N);
    auto M0_struct = createM0_structure(N);

    auto konst = factor(N);

    auto H0 = build_structure(H0_struct, konst);
    auto M0 = build_structure(M0_struct, konst); // ������ ������� � term ���������� �� factor

    cout << "\nH0 terms (codes):\n";
    print_structure_codes(H0_struct);
    //cout << "\nH0 terms (matrices):\n";
    //print_structure_matrices(H0);
    //
    cout << "\n num of factors :\n";
    printVector_double(konst);
    //
    //cout << "\nM0 terms (codes):\n";
    //print_structure_codes(M0_struct);
    //cout << "\nM0 terms (matrices):\n";
    //print_structure_matrices(M0);

    // ������� �� M
    Complex Sp_H0 = 0.0;
    Complex Sp_M0 = 0.0;
    Complex Sp = 0.0;

    for (size_t i = 0; i < static_cast<size_t>(M); i++) {
        //auto vec = random_normal_vector(dim, 0.0, 1.0); // ��������� ���������� �������
        //auto vec = Hutchinson_vector(dim);
        auto vec = generateRademacherComplexVector(dim);
        cout << "\n";
        cout << i << " vector = ";
        printVector(vec); // print 
        cout << "\n";

        Matrix2x2 Sum_H0(dim, 0.0);
        Matrix2x2 Sum_M0(dim, 0.0);
        for (size_t j = 0; j < k; j++) {
            if (j == 0) {
                auto start_one_kron = high_resolution_clock::now();
                Sum_H0 = kron_mult(H0, vec);    // ���������!
                if (Betta != 0) {
                    Sum_M0 = kron_mult(M0, vec); 
                }
                auto end_one_kron = high_resolution_clock::now();
                auto duration_2 = duration_cast<milliseconds>(end_one_kron-start_one_kron);
                cout << "Execution time of one kron_mult : " << duration_2.count() << " ms" << endl;
            }
            else {
                Sum_H0 = kron_mult(H0, Sum_H0);
                if (Betta != 0) {
                    Sum_M0 = kron_mult(M0, Sum_M0);
                }
            }
        }

        if (Betta == 0) {
            auto dot_product_H0 = dot_product(vec, Sum_H0); // �������� �� ����������������� ������ // ���������! 
            Sp_H0 += dot_product_H0;
        }
        else {
            auto dot_product_H0 = dot_product(vec, Sum_H0); // �������� �� ����������������� ������
            auto dot_product_M0 = dot_product(vec, Sum_M0);
            Sp_H0 += dot_product_H0;
            Sp_M0 += dot_product_M0;
        }
    }

    //auto coeff = Complex(0.5, 0.0) / static_cast<double>(M);
    double coeff = 1.0 / 200.0;
    Sp = coeff * Sp_H0 + Betta * coeff * Sp_M0;

    cout << std::fixed << std::setprecision(4) << coeff; // �������: 0.0005
    cout << "\n Sp:" << Sp << "\n";
    auto end_total = high_resolution_clock::now();
    auto duration_1 = duration_cast<seconds>(end_total - start_total);
    //auto duration_2 = duration_cast<seconds>(end_one_kron - end_one_kron);
    cout << "Execution total time: " << duration_1.count() << " s" << endl;
    //cout << "Execution time of one kron_mult : " << duration_2.count() << " s" << endl;

    return 0;
 }
//Complex(1.0, 0.0) / static_cast<double>(M);
//<< std::cout.precision(3)