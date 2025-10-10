#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <algorithm>


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
            konst.push_back(-2.0*factor);
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
    std::cout << "\n[" << M[0] << ", " << M[1] << "], \n "
        << "[" << M[2] << ", " << M[3] << "]\n";
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

void print_structure_term(const vector< Matrix2x2>& term) {
        for (const auto& M : term) {
            print_matrix(M);
            std::cout << "   ";
        }
        std::cout << "\n";
}


// --- ( A0 ? ... ? An ) * vec --- произведение Кронеккера 

void kron_mult_recursive_old(   ///////////// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! AAAAAAAAAAAAAAAAAAAAAAAA
    const std::vector<Matrix2x2>& factors,
    size_t N,
    const Matrix2x2& vec_in, // ссылка на константу 
    Matrix2x2& vec_out,
    size_t depth,
    size_t offset_in,
    size_t offset_out,
    size_t stride)
{
    if (depth == N) {
        //cout << string(depth * 2, ' ') << "БАЗА: depth=" << depth
        //    << " offset_in=" << offset_in
        //    << " offset_out=" << offset_out << endl;
        //
        vec_out[offset_out] = vec_in[offset_in];
        return;
    }

    size_t next_stride = stride / 2;
    const auto& M = factors[depth]; ///// factors - вектор матриц, M - это уже ссылка на адрес одной матрицы 

    // 

    Matrix2x2 temp0(next_stride, 0.0);
    Matrix2x2 temp1(next_stride, 0.0);
   

    //cout << string(depth * 2, ' ')
    //    << "? Спуск рекурсии depth=" << depth
    //    << " stride=" << stride << "\n";
      
    kron_mult_recursive_old(factors, N, vec_in, temp0, depth + 1, // при вызове kron я копирую вектор temp 
        offset_in, 0, next_stride);
    kron_mult_recursive_old(factors, N, vec_in, temp1, depth + 1,
        offset_in + next_stride, 0, next_stride);

    //cout << string(depth * 2, ' ')         
    //    << "? Выполняется цикл for на depth=" << depth << "\n";


    for (size_t i = 0; i < 2; ++i) {
        for (size_t idx = 0; idx < next_stride; ++idx) {
            Complex  val = M[i * 2 + 0] * temp0[idx] + M[i * 2 + 1] * temp1[idx];
            vec_out[offset_out + i * next_stride + idx] = val;  // vec_out - это вектор предыдущего вызова: vec_out = {  (temp0) , (temp1) }
         
            //cout << string(depth * 2, ' ') << "(" << M[i * 2 + 0] << " * " << temp0[idx] << ")"
            //    << " + "
            //    << "(" << M[i * 2 + 1] << " * " << temp1[idx] << ")"
            //    << " = " << val << "\n";
            //
            //cout << string(depth * 2, ' ') << "vec_out[" << offset_out << "+" << i << "*" << next_stride << "+ " << idx << "]" << " = " << val << "\n";
        }
    }
}


//
void kron_mult_recursive(
    const std::vector<Matrix2x2>& factors,
    size_t N,
    const Matrix2x2& vec_in,
    Matrix2x2& vec_out,
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
    const auto& M = factors[depth];

    size_t temp0_offset = offset_out;
    size_t temp1_offset = offset_out + next_stride;

    kron_mult_recursive(factors, N, vec_in, vec_out, depth + 1,
        offset_in, temp0_offset, next_stride);
    kron_mult_recursive(factors, N, vec_in, vec_out, depth + 1,
        offset_in + next_stride, temp1_offset, next_stride);

    for (size_t i = 0; i < 2; ++i) {
        for (size_t idx = 0; idx < next_stride; ++idx) {
            Complex  val = M[i * 2 + 0] * vec_out[temp0_offset + idx]
                + M[i * 2 + 1] * vec_out[temp1_offset + idx];
            vec_out[offset_out + i * next_stride + idx] = val;
        }
    }
}

void kron_mult_recursive_double(
    const std::vector<Matrix2x2>& factors,
    size_t depth,
    const Complex* src,
    Complex* dst,
    size_t stride,
    size_t dim
) {
    if (depth == factors.size()) {
        // базовый случай: просто копируем вектор
        std::copy(src, src + dim, dst);
        return;
    }

    const auto& M = factors[depth];
    size_t half = stride / 2;

    // применяем матрицу M блок за блоком
    for (size_t block = 0; block < dim; block += stride) {
        for (size_t idx = 0; idx < half; ++idx) {
            Complex x0 = src[block + idx];
            Complex x1 = src[block + half + idx];

            dst[block + idx] = M[0] * x0 + M[1] * x1;
            dst[block + half + idx] = M[2] * x0 + M[3] * x1;
        }
    }

    // рекурсивный вызов — двигаемся глубже
    kron_mult_recursive_double(factors, depth + 1, dst, const_cast<Complex*>(src), half, dim);
};



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


vector<Complex> random_normal_vector(size_t dim, double mean, double stddev) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dist(mean, stddev);

    vector<Complex> vec(dim);
    for (size_t i = 0; i < dim; ++i) {
        double re = (dist(gen));  // действительная часть
        double im = 0.0; // мнимая часть
        vec[i] = Complex(re, im);
    }
    return vec;
}

vector<Complex> Hutchinson_vector(size_t n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 1); // Use int!

    std::vector<Complex> result(n);
    for (size_t i = 0; i < n; ++i) {
        double rademacher_1 = dist(gen);
        double rademacher_2 = 0;
        result[i] = Complex(rademacher_1, rademacher_2); // Real part is ±1, imaginary part is 0
    }
    return result;
}

void printVector(const vector<Complex>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << "v[" << i << "] = "
            << v[i].real() << " + " << v[i].imag() << std::endl;
    }
}

//<< v[i].real() << " + " << v[i].imag() << "i"

void printVector_double(const std::vector<double>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << "v[" << i << "] = "
            << v[i]
            << std::endl;
    }
}


vector<Complex> generateRademacherComplexVector(size_t dimension) {
    std::vector<Complex> result;
    result.reserve(dimension);

    std::random_device rd{};
    std::mt19937 gen{ rd() };
    std::uniform_int_distribution<int> dist(0, 1);

    auto rademacher = [&]() -> double {
        return dist(gen) == 0 ? -1.0 : 1.0;
        };

  
    for (size_t i = 0; i < dimension; ++i) {
        double real_part = rademacher();  
        double imag_part = 0.0;  
        result.emplace_back(real_part, imag_part);
    }

    return result;
}