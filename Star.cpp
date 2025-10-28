                       // H0 - [ (A0, A1, ...., An), (A0, A1, ...., An), ... ]
                       // term = (A0, A1, ...., An)
                       // M = An
#include <iomanip>
#include "Matrix.h"  // подключаем функции из Matrix.cpp
#include <chrono>

using namespace std;
using namespace std::chrono;


 //--- Обёртка ---
std::vector<Complex> kron_mult_wrapper(const std::vector<std::vector<Matrix2x2>>& factors, const std::vector<Complex>& vec) {
    size_t N = factors[0].size();
    size_t dim = 1ULL << N;

    std::vector<Complex> Sum_0(dim, Complex(0.0, 0.0));
    std::vector<Complex> src(dim), dst(dim), tmp(dim);

    // для каждого терма копируем исходный vec в src
    for (const auto& term : factors) {
        std::copy(vec.begin(), vec.end(), src.begin());
        kron_mult_iterative(term, src.data(), dst.data(), tmp.data(), dim);
        for (size_t i = 0; i < dim; ++i) Sum_0[i] += dst[i];
    }
    return Sum_0;
}



//--- main ---

int main() {

    unsigned int cores = thread::hardware_concurrency();
    if (cores == 0) cores = 1;
    omp_set_num_threads(cores);

    std::cout << "Using " << cores << " CPU threads\n";
    auto start_total = high_resolution_clock::now();

    time_point<high_resolution_clock> start_one_kron;
    time_point<high_resolution_clock> end_one_kron;

    size_t k = 2; // power
    double M = 1000.0; // samples
    double Betta = 1.0;
    size_t N = 2;  // кол-во частиц (кубитов)
    int dim = 1 << N;
    cout << "Building H0 and M0 structures...\n";

    auto H0_struct = createH0_structure(N);
    auto M0_struct = createM0_structure(N);

    auto konst = factor(N);

    auto H0 = build_structure(H0_struct, konst);
    auto M0 = build_structure(M0_struct); 

    auto H = combine(H0, M0, Betta);

    //cout << "\nH0 terms (codes):\n";
    //print_structure_codes(H0_struct);
    //cout << "\nH0 terms (matrices):\n";
    //print_structure_matrices(H0);
    //
    //cout << "\n num of factors :\n";
    //printVector_double(konst);
    //
    cout << "\nM0 terms (codes):\n";
    print_structure_codes(M0_struct);
    cout << "\nM0 terms (matrices):\n";
    print_structure_matrices(M0);

    // Выборка из M
    Complex Sp_H0 = 0.0;
    Complex Sp_M0 = 0.0;
    Complex Sp = 0.0;

    for (size_t i = 0; i < static_cast<size_t>(M); i++) {
        auto vec = random_normal_vector(dim, 0.0, 1.0); // генерация случайного вектора
       //auto vec = Hutchinson_vector(dim);
        //auto vec = generateRademacherComplexVector(dim); - фигня полная 
        //vector<Complex> vec(dim, Complex(1.0, 0.0));
        //cout << "\n";
        //cout << i << " vector = ";
        //printVector(vec); // print  
        //cout << "\n";   

        Matrix2x2 Sum_H0(dim, 0.0);
        Matrix2x2 Sum_M0(dim, 0.0);
        for (size_t j = 0; j < k; j++) {
            if (j == 0) {
                auto start_one_kron = high_resolution_clock::now();
                Sum_H0 = kron_mult_wrapper(H, vec);
                //cout << "\n Sum_H0:\n";
                //printVector(Sum_H0);// корректно?
                auto end_one_kron = high_resolution_clock::now();
                //auto duration_2 = duration_cast<milliseconds>(end_one_kron-start_one_kron);
                //cout << "Execution time of one kron_mult : " << duration_2.count() << " ms" << endl;
            }
            else {
                Sum_H0 = kron_mult_wrapper(H, Sum_H0);

            }
        }

        auto dot_product_H0 = dot_product(vec, Sum_H0); // Умножеие на транспонированный вектор // корректно! 
        Sp_H0 += dot_product_H0;

    }

    //auto coeff = Complex(0.5, 0.0) / static_cast<double>(M);
    Complex coeff = Complex(1.0, 0.0) / M;
    Sp = coeff * Sp_H0;

    cout << std::fixed << std::setprecision(4) << coeff; // Выведет: 0.0005
    cout << "\n Sp:" << Sp << "\n";
    auto end_total = high_resolution_clock::now();
    auto duration_1 = duration_cast<seconds>(end_total - start_total);
    //auto duration_2 = duration_cast<seconds>(end_one_kron - end_one_kron);
    cout << "Execution total time: " << duration_1.count() << " s" << endl;
    //cout << "Execution time of one kron_mult : " << duration_2.count() << " s" << endl;


    // ------------------------- тест ---------------------------- //

    {
        cout << "\n ------------- TEST --------------- \n";
        Complex Sp_H0 = 0.0;
        Complex Sp = 0.0;
    
        auto H_st = build_H0_matrix(H, N);
        auto H2 = power_H0_matrix(H_st, 2);
    
        cout << "\nH0 matrix (" << (1 << N) << "x" << (1 << N) << "):\n";
        //print_full_matrix(H);
    
            
        for (size_t j = 0; j < static_cast<size_t>(M); j++) {
    
            auto vec = random_normal_vector(dim, 0.0, 1.0);
            //vector<Complex> vec(dim, Complex(1.0, 0.0));
    
            auto dot_product_H0 = dot_product_1(H2, vec); // - вектор 
            //cout << "\n dot_product_H:\n";
            //printVector(dot_product_H0);
            auto dot_product_T_H0 = dot_T_product(vec, dot_product_H0); // - число 
    
            Sp_H0 += dot_product_T_H0;
        }
    
        Complex coeff = Complex(1.0, 0.0) / M;
        Sp = coeff * Sp_H0;
    
        cout << std::fixed << std::setprecision(4) << coeff; // Выведет: 0.0005
        cout << "\n Sp:" << Sp << "\n";
    
        Complex tr = trace_matrix(H2, dim);
        std::cout << "Trace of H: " << tr << "\n";
    
    }
    return 0;
 }

 //if (it == factors.begin()) {
//    cout << "\n ------------------------ \n";
//}
//cout << "\n";
//print_structure_term(term);
//cout << "\n";
 //kron_mult_recursive_double(term, 0,
 //    src.data(),      // вход
 //    dst.data(),      // выход
 //    tmp.data(),      // временный буфер
 //    dim, dim);
 //kron_mult_recursive_double(term, 0, vec.data(), result.data(), dim, dim);
 //Sum_0 = add_vectors(Sum_0, dst);
 //cout << "\n Sum_0:\n";
 //printVector(Sum_0);

 //auto H = build_term_matrix(term, N);
 //auto dot_product_H0 = dot_product_1(H, vec);
 //cout << "\n dot_product_H:\n";
 //printVector(dot_product_H0);

 // Для последнего элемента
 //if (it == std::prev(factors.end())) {
 //    cout << "\n ------------------------ \n";
 //}