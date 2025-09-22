#include <iostream>
#include <vector>
#include <complex>
#include <random>

using namespace std;


// 0: I, 1: X, 2: Y, 3: Z
enum PauliType { I = 0, X = 1, Y = 2, Z = 3 };

using Complex = std::complex<double>;
using Matrix2x2 = std::vector<Complex>;
using Structure = std::vector<std::vector<Matrix2x2>>;

// --- Матрицы Паули ---
Matrix2x2 pauli_I() { return { 1.0, 0.0, 0.0, 1.0 }; }
Matrix2x2 pauli_X() { return { 0.0, 1.0, 1.0, 0.0 }; }
Matrix2x2 pauli_Y() {
    using namespace std::complex_literals;
    return { 0.0, -1.0i, 1.0i, 0.0 };
}
Matrix2x2 pauli_Z() { return { 1.0, 0.0, 0.0, -1.0 }; }

Matrix2x2 pauli_matrix(int code) {
    switch (code) {
    case I: return pauli_I();
    case X: return pauli_X();
    case Y: return pauli_Y();
    case Z: return pauli_Z();
    default: return pauli_I();
    }
}

std::vector<double>  factor(int N) {
    std::vector<double> konst;
    for (int j = 0; j < N; ++j) {
        for (int k = j + 1; k < N; ++k) {
            double factor = 1.0 / pow(k - j, 3);
            konst.push_back(factor);
            konst.push_back(factor);
            konst.push_back(factor);
        }
    }
    return konst;
}

// --- Вспомогательные функции ---
std::vector<std::vector<int>> createH0_structure(int N) {
    std::vector<std::vector<int>> vec_full;
    for (int j = 0; j < N; ++j) {
        for (int k = j + 1; k < N; ++k) {
            for (int sig = 1; sig <= 2; ++sig) {
                std::vector<int> vec_state(N, 0);
                vec_state[j] = sig;
                vec_state[k] = sig;
                vec_full.push_back(vec_state);
            }
            std::vector<int> vec_state(N, 0);
            vec_state[j] = Z;
            vec_state[k] = Z;
            vec_full.push_back(vec_state);
        }
    }
    return vec_full;
}

std::vector<std::vector<int>> createM0_structure(int N) {
    std::vector<std::vector<int>> vec_full;
    for (int j = 0; j < N; ++j) {
        std::vector<int> vec_state(N, 0);
        vec_state[j] = Z;
        vec_full.push_back(vec_state);
    }
    return vec_full;
}


Structure build_structure(const std::vector<std::vector<int>>& vec_full,
    const std::vector<double>& factors) {
    Structure H0;
    for (size_t term_idx = 0; term_idx < vec_full.size(); ++term_idx) {
        const auto& v = vec_full[term_idx]; // vec_state
        std::vector<Matrix2x2> term;
        for (size_t i = 0; i < v.size(); ++i) {
            Matrix2x2 M = pauli_matrix(v[i]);
            if (i == 0) { // умножаем только первую матрицу
                for (auto& x : M) x *= factors[term_idx];
            }
            term.push_back(M);
        }
        H0.push_back(term);
    }
    return H0; //  -  вектор матриц
}


// --- Печать матрицы ---
void print_matrix(const Matrix2x2& M) {
    std::cout << "[[" << M[0] << ", " << M[1] << "], "
        << "[" << M[2] << ", " << M[3] << "]]";
}

// --- Печать структуры кодов (I, X, Y, Z) ---
void print_structure_codes(const std::vector<std::vector<int>>& vec_full) {
    for (const auto& v : vec_full) {
        for (int p : v) {
            switch (p) {
            case I: std::cout << "I "; break;
            case X: std::cout << "X "; break;
            case Y: std::cout << "Y "; break;
            case Z: std::cout << "Z "; break;
            }
        }
        std::cout << "\n";
    }
}

// --- Печать структуры матриц ---
void print_structure_matrices(const Structure& S) {
    for (const auto& term : S) {
        for (const auto& M : term) {
            print_matrix(M);
            std::cout << "   ";
        }
        std::cout << "\n";
    }
}

// --- ( A0 ? ... ? An ) * vec --- произведение Кронеккера 

void kron_mult_recursive(   ///////////// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! AAAAAAAAAAAAAAAAAAAAAAAA
    const std::vector<Matrix2x2>& factors,
    size_t N,
    const vector<Complex>& vec_in,
    vector<Complex>& vec_out,
    size_t depth,
    size_t offset_in,
    size_t offset_out,
    size_t stride
) {
    if (depth == N) {
        vec_out[offset_out] = vec_in[offset_in];
        return;
    }

    size_t next_stride = stride / 2;
    const auto& M = factors[depth]; ///// factors - вектор матриц, M - это уже ссылка на адрес одной матрицы 

    vector<Complex> temp0(next_stride, 0.0);
    vector<Complex> temp1(next_stride, 0.0);

    kron_mult_recursive(factors, N, vec_in, temp0, depth + 1,
        offset_in, 0, next_stride);
    kron_mult_recursive(factors, N, vec_in, temp1, depth + 1,
        offset_in + next_stride, 0, next_stride);

    for (size_t i = 0; i < 2; ++i) {
        for (size_t idx = 0; idx < next_stride; ++idx) {
            Complex val = M[i * 2 + 0] * temp0[idx] + M[i * 2 + 1] * temp1[idx];
            vec_out[offset_out + i * next_stride + idx] = val;
        }
    }
}

std::vector<Complex> add_vectors(const std::vector<Complex>& a, const std::vector<Complex>& b) {
    std::vector<Complex> res(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        res[i] = a[i] + b[i];
    }
    return res;
}

Complex dot_product(const std::vector<Complex>& a, const std::vector<Complex>& b) {
    if (a.size() != b.size()) {
        throw std::runtime_error("Размеры векторов не совпадают");
    }
    Complex sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        sum += std::conj(a[i]) * b[i];  
    }
    return sum;
}

Matrix2x2 random_normal_vector(size_t dim, double mean, double stddev) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dist(mean, stddev);

    Matrix2x2 vec(dim);
    for (size_t i = 0; i < dim; ++i) {
        double re = dist(gen);  // действительная часть
        double im = dist(gen);  // мнимая часть
        vec[i] = Complex(re, im);
    }
    return vec;
}

void printVector(const std::vector<Complex>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << "v[" << i << "] = "
            << v[i].real() << " + " << v[i].imag() << "i"
            << std::endl;
    }
}

void printVector_double(const std::vector<double>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << "v[" << i << "] = "
            << v[i]
            << std::endl;
    }
}