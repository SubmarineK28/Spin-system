#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <algorithm>
#include <thread>
#include <omp.h>
#include <algorithm>

using namespace std;


// 0: I, 1: X, 2: Y, 3: Z
enum PauliType { I = 0, X = 1, Y = 2, Z = 3 };

using Complex = complex<double>;
using Matrix2x2 = vector<Complex>;
using Structure = vector<vector<Matrix2x2>>;

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
        vector<Matrix2x2> term;
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

Structure build_structure(const std::vector<std::vector<int>>& vec_full) {
    Structure H0;
    for (size_t term_idx = 0; term_idx < vec_full.size(); ++term_idx) {
        const auto& v = vec_full[term_idx]; // vec_state
        vector<Matrix2x2> term;
        for (size_t i = 0; i < v.size(); ++i) {
            Matrix2x2 M = pauli_matrix(v[i]);
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

void print_full_matrix(const Matrix2x2& M) {
    size_t dim = static_cast<size_t>(std::sqrt(M.size()));

    std::cout << "\nMatrix (" << dim << "x" << dim << "):\n";

    for (size_t i = 0; i < dim; ++i) {
        std::cout << "[ ";
        for (size_t j = 0; j < dim; ++j) {
            const Complex& val = M[i * dim + j];
            std::cout << "(" << val.real() << "," << val.imag() << ")";
            if (j + 1 < dim) std::cout << ", ";
        }
        std::cout << " ]\n";
    }
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

void kron_mult_iterative(
    const std::vector<Matrix2x2>& factors,
    const Complex* src_in,
    Complex* dst,
    Complex* tmp,
    size_t dim
) {
    // Рабочие указатели одного типа
    Complex* src = const_cast<Complex*>(src_in);
    Complex* buf = dst; // будем писать в buf
    size_t N = factors.size();

    // копируем входной src_in в src (src и src_in могут совпадать — caller делает копию)
    // caller гарантирует, что src points to modifiable buffer (we used copy before call)

    for (size_t depth = 0; depth < N; ++depth) {
        const auto& M = factors[depth];
        size_t stride = size_t(1) << (N - depth - 1); // блок длины stride, total block 2*stride

    #pragma omp parallel for
        for (long long block = 0; block < (long long)dim; block += 2LL * (long long)stride) {
            for (size_t i = 0; i < stride; ++i) {
                Complex a = src[block + i];
                Complex b = src[block + i + stride];
                buf[block + i] = M[0] * a + M[1] * b;
                buf[block + i + stride] = M[2] * a + M[3] * b;
            }
        }

        // swap roles: next iteration read from buf, write into other buffer
        // we will use tmp to alternate: src <-> buf <-> tmp
        std::swap(src, buf);
        // now buf points to the other buffer (either dst or src_in's buffer)
        // ensure we have somewhere to write on next iteration: if buf == dst then next write uses src buffer, etc.
        // We'll maintain that caller provided two buffers (src_in as modifiable copy and dst), or we used tmp externally.
        if (src == dst) {
            // we want buf to point to tmp for next write
            buf = tmp;
        }
        else {
            buf = dst; // write back to dst if possible
        }
    }

    // after loop, src points to the buffer that contains final result
    if (src != dst) {
        std::copy(src, src + dim, dst);
    }
}

// multiply complex scalar m (mr,mi) by two complex numbers packed as [ar0,ai0,ar1,ai1] in __m256d
static inline __m256d cmplx_scalar_mul_2(const __m256d vec, const double mr, const double mi) {
    // vec = [ar0, ai0, ar1, ai1]
    __m256d mr_vec = _mm256_set1_pd(mr);    // [mr,mr,mr,mr]
    __m256d mi_vec = _mm256_set1_pd(mi);    // [mi,mi,mi,mi]

    __m256d perm = _mm256_permute_pd(vec, 0x5); // permute pairs -> [ai0, ar0, ai1, ar1]
    __m256d t1 = _mm256_mul_pd(mr_vec, vec);    // mr * [ar,ai,...]
    __m256d t2 = _mm256_mul_pd(mi_vec, perm);   // mi * [ai,ar,...]

    // sign vector s = [-1, +1, -1, +1]
    const __m256d sign = _mm256_set_pd(1.0, -1.0, 1.0, -1.0); // note order in set_pd is high->low; we'll multiply perm accordingly
    // But easier: we want multiply t2 by [-1, +1, -1, +1] to add/subtract appropriately.
    // Create s2 such that after multiply we have desired signs at positions [0..3].
    const __m256d s2 = _mm256_setr_pd(-1.0, 1.0, -1.0, 1.0); // setr gives order low->high
    __m256d t2s = _mm256_mul_pd(t2, s2);

    return _mm256_add_pd(t1, t2s); // mr*vec + s2*(mi*perm)
}

// AVX2-optimized iterative Kron mult: factors size N (number of qubits), dim = 2^N
void kron_mult_iterative_avx2(
    const std::vector<Matrix2x2>& factors,
    const Complex* src_in, // caller should provide modifiable src buffer (copy of input vec)
    Complex* dst,
    Complex* tmp,
    size_t dim
) {
    // make modifiable pointer
    Complex* src = const_cast<Complex*>(src_in);
    Complex* buf = dst; // write here
    size_t N = factors.size();

    for (size_t depth = 0; depth < N; ++depth) {
        const auto& M = factors[depth];
        // expect M.size() == 4, matrix elements are Complex
        double m00r = M[0].real(), m00i = M[0].imag();
        double m01r = M[1].real(), m01i = M[1].imag();
        double m10r = M[2].real(), m10i = M[2].imag();
        double m11r = M[3].real(), m11i = M[3].imag();

        size_t stride = size_t(1) << (N - depth - 1); // half-block length (called stride earlier)
        size_t fullBlock = 2 * stride; // the block size we step over

        // Parallelize over blocks
    #pragma omp parallel for schedule(static)
        for (long long block = 0; block < (long long)dim; block += (long long)fullBlock) {
            size_t b = (size_t)block;

            // process in chunks of 2 complex numbers (4 doubles) using AVX2
            size_t i = 0;
            for (; i + 1 < stride; i += 2) { // i indexes complex slots within half-block
                size_t idx0 = b + i;
                // addresses for a and b vectors (two complex numbers each)
                const double* a_ptr = reinterpret_cast<const double*>(&src[idx0]);                // points to re0,im0,...
                const double* b_ptr = reinterpret_cast<const double*>(&src[idx0 + stride]);     // second half longs

                // load 2 complex numbers for a: [ar0, ai0, ar1, ai1]
                __m256d va = _mm256_loadu_pd(a_ptr);
                // load 2 complex numbers for b
                __m256d vb = _mm256_loadu_pd(b_ptr);

                // compute m00 * a
                __m256d r_m00_a = cmplx_scalar_mul_2(va, m00r, m00i);
                // compute m01 * b
                __m256d r_m01_b = cmplx_scalar_mul_2(vb, m01r, m01i);
                // sum => upper part
                __m256d upper = _mm256_add_pd(r_m00_a, r_m01_b);

                // compute m10 * a
                __m256d r_m10_a = cmplx_scalar_mul_2(va, m10r, m10i);
                // compute m11 * b
                __m256d r_m11_b = cmplx_scalar_mul_2(vb, m11r, m11i);
                __m256d lower = _mm256_add_pd(r_m10_a, r_m11_b);

                // store upper into buf[idx0 .. idx0+1] (2 complex doubles)
                _mm256_storeu_pd(reinterpret_cast<double*>(&buf[idx0]), upper);
                // store lower into buf[idx0 + stride .. idx0 + stride +1]
                _mm256_storeu_pd(reinterpret_cast<double*>(&buf[idx0 + stride]), lower);
            }

            // Tail: if stride is odd, or leftover one complex (i < stride)
            for (; i < stride; ++i) {
                size_t idx = b + i;
                Complex a = src[idx];
                Complex bb = src[idx + stride];
                Complex up = M[0] * a + M[1] * bb;
                Complex low = M[2] * a + M[3] * bb;
                buf[idx] = up;
                buf[idx + stride] = low;
            }
        } // end parallel blocks

        // swap src and buf for next depth; ensure buf becomes a valid place to write next
        std::swap(src, buf);

        // set buf to point to an available buffer (either tmp or dst)
        if (buf == dst) {
            buf = tmp;
        }
        else {
            buf = dst;
        }
    } // end depth loop

    // result is in src pointer; ensure it's in dst
    if (src != dst) {
        std::copy(src, src + dim, dst);
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
        double rademacher_2 = 0.0;
        result[i] = Complex(rademacher_1, rademacher_2); // Real part is ±1, imaginary part is 0
    }
    return result;
}
//
//std::vector<Complex> generate_rademacher_vector(size_t dim) {
//    std::random_device rd;
//    std::mt19937 gen(rd());
//    std::uniform_int_distribution<int> dist(0, 1); // 0 или 1
//
//    std::vector<Complex> zi(dim);
//    for (size_t i = 0; i < dim; ++i) {
//        int r = dist(gen);
//        zi[i] = Complex(r == 0 ? -1.0 : 1.0, 0.0); // ±1
//    }
//    return zi;
//}

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
        double imag_part = rademacher();
        result.emplace_back(real_part, imag_part);
    }

    return result;
}


// ------------------------- тест ---------------------------- //

Matrix2x2 kron(const Matrix2x2& A, const Matrix2x2& B) {
    size_t nA = static_cast<size_t>(std::sqrt(A.size()));
    size_t nB = static_cast<size_t>(std::sqrt(B.size()));
    size_t nC = nA * nB;

    Matrix2x2 C(nC * nC);

    for (size_t i = 0; i < nA; ++i) {
        for (size_t j = 0; j < nA; ++j) {
            Complex a = A[i * nA + j];
            for (size_t k = 0; k < nB; ++k) {
                for (size_t l = 0; l < nB; ++l) {
                    C[(i * nB + k) * nC + (j * nB + l)] = a * B[k * nB + l];
                }
            }
        }
    }
    return C;
}

// --- Построение полной матрицы H0 ---
Matrix2x2 build_H0_matrix(const Structure& H0_terms, int N) {
    size_t dim = 1ULL << N;                  // размерность итоговой матрицы
    Matrix2x2 H(dim * dim, 0.0);             // плоская матрица dim x dim

    cout << "Number of terms: " << H0_terms.size() << endl;
    for (const auto& term : H0_terms) {
        cout << "Term has " << term.size() << " matrices" << endl;// начинаем с первой матрицы терма
        Matrix2x2 prod = term[0];

        // перемножаем все матрицы терма
        for (size_t i = 1; i < term.size(); ++i) {
            prod = kron(prod, term[i]);
        }

        // суммируем терм в общую матрицу H
        for (size_t i = 0; i < H.size(); ++i) {
            H[i] += prod[i];
        }
    }

    return H;
}

Matrix2x2 build_term_matrix(const vector<vector<Complex>>& term, int N) {
    size_t dim = 1ULL << N;                  // размерность итоговой матрицы
    Matrix2x2 H(dim * dim, 0.0);             // плоская матрица dim x dim
 
    Matrix2x2 prod = term[0];

    // перемножаем все матрицы терма
    for (size_t i = 1; i < term.size(); ++i) {
        prod = kron(prod, term[i]);
    }

    // суммируем терм в общую матрицу H
    for (size_t i = 0; i < H.size(); ++i) {
        H[i] += prod[i];
    }

    return H;
}

// -- - Матрично - векторное произведение : A * vec-- -
vector<Complex> dot_product_1(const Matrix2x2& A, const std::vector<Complex>& vec) {
    size_t n = static_cast<size_t>(std::sqrt(A.size()));
    std::vector<Complex> result(n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        Complex sum = 0.0;
        for (size_t j = 0; j < n; ++j) {
            sum += A[i * n + j] * vec[j];
        }
        result[i] = sum;
    }

    return result;
}

// --- Скалярное произведение: (vec*)^T * Av ---
Complex dot_T_product(const std::vector<Complex>& vec, const std::vector<Complex>& Av) {
    if (vec.size() != Av.size()) {
        throw std::runtime_error("dot_T_product: vector size mismatch");
    }

    Complex sum = 0.0;
    for (size_t i = 0; i < vec.size(); ++i) {
        sum += std::conj(vec[i]) * Av[i];
    }

    return sum;
}

// --- Перемножение двух квадратных матриц ---
Matrix2x2 multiply_matrix(const Matrix2x2& A, const Matrix2x2& B) {
    size_t n = static_cast<size_t>(std::sqrt(A.size()));
    Matrix2x2 C(n * n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            Complex sum = 0.0;
            for (size_t k = 0; k < n; ++k) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
    return C;
}

// --- Единичная матрица ---
Matrix2x2 identity_matrix(size_t n) {
    Matrix2x2 I(n * n, 0.0);
    for (size_t i = 0; i < n; ++i)
        I[i * n + i] = 1.0;
    return I;
}

// --- Возведение в степень H^k ---
Matrix2x2 power_H0_matrix(Matrix2x2 H, int k) {
    if (k < 0) {
        throw std::runtime_error("power_H0_matrix: negative powers not supported");
    }

    size_t n = static_cast<size_t>(std::sqrt(H.size()));
    Matrix2x2 result = identity_matrix(n);

    while (k > 0) {
        if (k % 2 == 1) {
            result = multiply_matrix(result, H);
        }
        H = multiply_matrix(H, H);
        k /= 2;
    }

    return result;
}

Complex trace_matrix(const std::vector<Complex>& H, size_t n) {
    Complex trace = 0.0;

    for (size_t i = 0; i < n; ++i) {
        trace += H[i * n + i]; // диагональный элемент
    }

    return trace;
}


Structure combine(
    const Structure& H0,
    const Structure& M0,
    double beta
) {
    if (std::abs(beta) < 1e-15) {
        return H0; // Beta == 0 ? Ничего не меняем
    }

    Structure H;
    H.reserve(H0.size() + M0.size());

    // Копируем все из H0
    for (const auto& term : H0) {
        H.push_back(term);
    }

    // Добавляем -beta * (только первую матрицу каждого term в M0)
    for (const auto& term : M0) {
        std::vector<Matrix2x2> processed_term = term; // копия

        //if (!processed_term.empty()) {
        //    for (auto& c : processed_term[0]) {
        //        c *= Complex(-beta, 0.0); // умножаем только первую матрицу
        //    }
        //}
        
        H.push_back(std::move(processed_term));
    }

    return H;
}

