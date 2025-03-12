#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <bitset>
#include <iomanip>
#include <algorithm> 
#include <map>
#include <climits>
#include <limits>
#include <set>
#include <time.h>
#include <omp.h>

using namespace std;

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void Parser(ifstream& input1, ifstream& input2, ifstream& input3, int& blocks, int& iter, double& SNR, long& SEED, int& algo);
double ran1(long* idum);
void normal(double& n1, double& n2, double sigma, long* idum);
void U_generator(int*& u, const int prev_u[6]);
void U_enc_mod(const int* u, int*& c, int*& x);
void AWGN(double*& y, const int* x, const double sigma, long* idum);
double CHK(int algo, double L1, double L2);
bool codeword(const int* x_est, const int* c);
void initialization(const double* y, double sigma, double*& Ll, double (*qml)[1023]);
void BottomUp(int algo, const double (*qml)[1023], double (*rml)[1023]);
void TopDown(const double* Ll, double (*qml)[1023], const double (*rml)[1023]);
void Termination(int* x_est, const double* Ll, const double (*rml)[1023], double*& ql);
void getH(void);
void BifFlipping(double*& y, int*& u, int*& x, int*& c, double sigma, long* idum, int& terminate_b, int& err_bit, int& err_block, double& ber, double& bler);

int blocks, iter, algo;
long SEED;
double SNR;
double threshold = 0.99;    // for BF algo, simulation for 0.99, 1.99, 3.99
int (*G)[1024] = new int[781][1024];    // G[781][1024]
int (*H_s)[32] = new int[1023][32];    // H_s[1023][32]
int (*H_x)[32] = new int[1023][32];    // H_x[1023][32]
int (*H_s_real)[1023] = new int[1023][1023];    // H_x[1023][32]
int (*H_x_real)[1023] = new int[1023][1023];    // H_x[1023][32]

/* AWGN part */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

int main()
{
    // Start the timer to measure the execution time of the simulation
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    ifstream inputFile1("Sim.txt");
    ifstream inputFile2("ldpc_G_1023.txt");
    ifstream inputFile3("ldpc_H_1023.txt");
    if (!inputFile1.is_open() || !inputFile2.is_open() || !inputFile3.is_open()) {
        cerr << "Error opening the input file..." << endl;
        return 1;
    }
    Parser(inputFile1, inputFile2, inputFile3, blocks, iter, SNR, SEED, algo);

    long* idum;
    idum = (long*)malloc(sizeof(long));
    *idum = SEED;
    double sigma = sqrt(1 / ((2.0 * 781.0 / 1023.0) * pow(10.0, SNR / 10.0)));       // R=781/1023 (#input/#output)

    /* decoding algo. part */
    int err_block = 0;
    int err_bit = 0;
    double ber = 0.0;
    double bler = 0.0;
    int x_est[1023];
    int terminate_b = 0;
    bool flag;
    int prev_u[6] = { 1, 0, 0, 0, 0, 0 };

    int* u, * c, * x;   //x=(-1)^c
    double* y;
    u = new int[781];
    c = new int[1023];
    x = new int[1023];
    y = new double[1023];

    /* SPA & MSA */
    if (algo != 2) {
        double* Ll = new double[1023];
        double (*qml)[1023] = new double[1023][1023];
        double (*rml)[1023] = new double[1023][1023];
        double* ql = new double[1023];
        for (int b = 0;b < blocks;b++) {
            U_generator(u, prev_u);
            U_enc_mod(u, c, x);
            AWGN(y, x, sigma, idum);

            initialization(y, sigma, Ll, qml);
            flag = false;  // true means decoding successfully
            for (int it = 0;it < iter;it++) {
                if (flag == false) {
                    BottomUp(algo, qml, rml);
                    TopDown(Ll, qml, rml);
                    Termination(x_est, Ll, rml, ql);
                    flag = codeword(x_est, c);
                }
                else
                    break;
            }
            for (int i = 0;i < 781;i++) {
                if (x_est[i] != c[i])
                    err_bit++;
            }
            if (flag == false)  err_block++;

            cout << "block @ " << b + 1 << endl;
            cout << "#error blocks = " << err_block << endl;
            cout << "#error bits = " << err_bit << endl;

            if (err_block == 100) {
                ber = err_bit / ((b + 1) * 781.0);
                bler = err_block / ((b + 1) * 1.0);
                terminate_b = b + 1;
                break;
            }
            // for next u sequence
            for (int i = 0; i < 6; i++) {
                prev_u[i] = u[781 - 6 + i];
            }
        }
        delete[] Ll;
        delete[] ql;
        delete[] qml;
        delete[] rml;
    }
    /* BF */
    else {
        getH();
        BifFlipping(y, u, x, c, sigma, idum, terminate_b, err_bit, err_block, ber, bler);
    }

    // result
    if (err_block < 100) {
        ber = err_bit / (blocks * 781.0);
        bler = err_block / (blocks * 1.0);
    }
    cout << "----------------------------------------------------" << endl;
    cout << "Number of decoded blocks = " << terminate_b << endl;
    cout << "Number of iterations = " << iter << endl;
    cout << "SNR = " << SNR << endl;
    cout << "Algorithm is " << algo << endl;
    if (algo == 2) cout << "Threshold = " << threshold << endl;
    cout << "# total error bits = " << err_bit << endl;
    cout << "# total error blocks = " << err_block << endl;
    cout << fixed << setprecision(8) << "BER = " << ber << endl;
    cout << fixed << setprecision(5) << "BLER = " << bler << endl;

    /* close the files */
    inputFile1.close();
    inputFile2.close();
    inputFile3.close();
    delete[] u;
    delete[] x;
    delete[] c;
    delete[] y;
    delete[] G;
    delete[] H_s;
    delete[] H_x;

    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("running time: %f sec\n", cpu_time_used);
    //system("pause");
    return 0;
}

void Parser(ifstream& input1, ifstream& input2, ifstream& input3, int& blocks, int& iter, double& SNR, long& SEED, int& algo) {
    // Sim.txt
    string line1;
    int line1Count = 0; // 用來追蹤讀到第幾行
    while (getline(input1, line1)) {
        size_t percentPos = line1.find('%'); // 找到 '%' 的位置
        if (percentPos != string::npos) {
            string value = line1.substr(0, percentPos); // 擷取 '%' 前的部分
            value.erase(value.find_last_not_of(" \t\r\n") + 1); // 移除多餘空白

            // 將數值存入對應的變數
            switch (line1Count) {
            case 0:
                blocks = stoi(value);
                break;
            case 1:
                iter = stoi(value);
                break;
            case 2:
                SNR = stod(value);
                break;
            case 3:
                SEED = stol(value);
                break;
            case 4:
                algo = stoi(value);
                break;
            default:
                cerr << "Unexpected line in Sim.txt" << endl;
                exit(EXIT_FAILURE);
            }
            line1Count++;
        }
    }
    // debug
    // cout << "Parser results: blocks=" << blocks << ", iter=" << iter << ", SNR=" << SNR
    //  << ", SEED=" << SEED << ", algo=" << algo << endl;

    // ldpc_G_1023.txt
    string line2;
    for (int i = 0; i < 781; i++) { // 遍歷每行
        if (getline(input2, line2)) {
            for (int j = 0; j < 1023; j++) {
                G[i][j] = line2[j] - '0'; // 將每個字元轉換為數字
            }
        }
    }

    // ldpc_H_1023.txt
    string line3;
    int row_index = 0;

    while (getline(input3, line3)) {
        // 將行解析為整數
        stringstream ss(line3);

        for (int col = 0; col < 32; col++) {
            int value;
            ss >> value;

            if (row_index < 1023) {
                H_s[row_index][col] = value; // 存入 H_s
            }
            else if (row_index < 2 * 1023) {
                H_x[row_index - 1023][col] = value; // 存入 H_x
            }
            else {
                cerr << "Unexpected extra rows in the file!" << endl;
                exit(EXIT_FAILURE);
            }
        }
        row_index++;
    }
}

double ran1(long* idum) {
    int j;
    long k;
    static long iy = 0;
    static long iv[NTAB];
    double temp;

    if (idum[0] <= 0 || iy == 0) {
        if (-(*idum) < 1)     *idum = 1;
        else     *idum = -(*idum);

        for (j = NTAB + 7; j >= 0; j--) {
            k = *idum / IQ;
            *idum = IA * (*idum - k * IQ) - IR * k;
            if (*idum < 0)
                *idum += IM;
            if (j < NTAB)
                iv[j] = *idum;
        }
        iy = iv[0];
    }

    k = *idum / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;
    if (*idum < 0)
        *idum += IM;

    j = iy / NDIV;
    iy = iv[j];
    iv[j] = *idum;

    if ((temp = AM * iy) > RNMX)
        return RNMX;
    else
        return temp;
}

void normal(double& n1, double& n2, double sigma, long* idum) {
    double temp1, temp2, s;
    do {
        temp1 = ran1(idum);
        temp2 = ran1(idum);
        temp1 = 2 * temp1 - 1;
        temp2 = 2 * temp2 - 1;
        s = temp1 * temp1 + temp2 * temp2;
    } while (s >= 1.0);
    n1 = sigma * temp1 * sqrt((-2.0) * log(s) / s);
    n2 = sigma * temp2 * sqrt((-2.0) * log(s) / s); // log ≡ ln in c++
}

void U_generator(int*& u, const int prev_u[6]) {
    u[0] = prev_u[0];
    u[1] = prev_u[1];
    u[2] = prev_u[2];
    u[3] = prev_u[3];
    u[4] = prev_u[4];
    u[5] = prev_u[5];
    for (int i = 6; i < 781; i++) {
        u[i] = u[i - 5] ^ u[i - 6];
    }
}

void U_enc_mod(const int* u, int*& c, int*& x) {
    /* encoder & modulator */
    for (int j = 0; j < 1023; j++) {
        c[j] = 0;    // initialize
        for (int i = 0; i < 781; i++) {
            c[j] ^= u[i] * G[i][j];
        }
        /* modulator */
        x[j] = (c[j] == 0) ? 1 : -1;
    }
    // debug
    // for (int i = 0;i < 1023 * blocks;i++) {
    //     cout << x[i] << " ";
    // }
}

void AWGN(double*& y, const int* x, const double sigma, long* idum) {
    double n1, n2;
    for (int i = 0;i < 1023;i += 2) {
        normal(n1, n2, sigma, idum);
        y[i] = x[i] + n1;
        if (i != 1022)  y[i + 1] = x[i + 1] + n2;
    }
    // debug
    // for (int i = 0;i < 1023 * blocks;i++) {
    //     cout << y[i] << " ";
    // }
}

double CHK(int algo, double L1, double L2) {
    switch (algo) {
    case 0:     //SPA(alter
        return sgn(L1) * sgn(L2) * min(abs(L1), abs(L2)) + log((1 + exp(-abs(L1 + L2))) / (1 + exp(-abs(L1 - L2))));
        break;
    case 1:     //MSA
        return sgn(L1) * sgn(L2) * min(abs(L1), abs(L2));
        break;
    }
}

bool codeword(const int* x_est, const int* c) {
    for (int length = 0;length < 1023;length++) {
        if (x_est[length] != c[length]) {
            return false;
        }
    }
    return true;
}

void initialization(const double* y, double sigma, double*& Ll, double (*qml)[1023]) {
    for (int i = 0;i < 1023;i++) {
        Ll[i] = 2.0 * y[i] / sigma / sigma;
    }

    // for H_x
    for (int m = 0;m < 1023;m++) {
        for (int l = 0;l < 32;l++) {
            qml[H_x[m][l] - 1][m] = Ll[m];
        }
    }
}

void BottomUp(int algo, const double (*qml)[1023], double (*rml)[1023]) {
#pragma omp parallel for collapse(2)
    for (int m = 0;m < 1023;m++) {
        for (int l = 0;l < 32;l++) {
            rml[m][H_s[m][l] - 1] = 0;
            for (int k = 0;k < 32;k++) {
                if (H_s[m][l] != H_s[m][k]) {
                    rml[m][H_s[m][l] - 1] = (rml[m][H_s[m][l] - 1] == 0) ? qml[m][H_s[m][k] - 1] : CHK(algo, rml[m][H_s[m][l] - 1], qml[m][H_s[m][k] - 1]);
                }
            }
        }
    }
}

void TopDown(const double* Ll, double (*qml)[1023], const double (*rml)[1023]) {
#pragma omp parallel for collapse(2)
    for (int m = 0;m < 1023;m++) {
        for (int l = 0;l < 32;l++) {
            qml[H_x[m][l] - 1][m] = Ll[m];
            for (int k = 0;k < 32;k++) {
                if (H_x[m][l] != H_x[m][k]) {
                    qml[H_x[m][l] - 1][m] += rml[H_x[m][k] - 1][m];
                }
            }
        }
    }
}

void Termination(int* x_est, const double* Ll, const double (*rml)[1023], double*& ql) {
    for (int m = 0;m < 1023;m++) {
        ql[m] = Ll[m];
        for (int l = 0;l < 32;l++) {
            ql[m] += rml[H_x[m][l] - 1][m];
        }
    }
    for (int i = 0;i < 1023;i++) {
        x_est[i] = (ql[i] > 0) ? 0 : 1;
    }
}

void getH(void) {
    for (int i = 0;i < 1023;i++) {
        for (int j = 0;j < 32;j++) {
            H_s_real[i][H_s[i][j] - 1] = 1;
            H_x_real[i][H_x[i][j] - 1] = 1;
        }
    }
}

void BifFlipping(double*& y, int*& u, int*& x, int*& c, double sigma, long* idum, int& terminate_b, int& err_bit, int& err_block, double& ber, double& bler) {
    int prev_u[6] = { 1,0,0,0,0,0 };
    for (int b = 0;b < blocks;b++) {
        U_generator(u, prev_u);
        U_enc_mod(u, c, x);
        AWGN(y, x, sigma, idum);

        int x_est[1023];
        for (int i = 0;i < 1023;i++) {
            x_est[i] = (y[i] > 0) ? 0 : 1;
        }

        bool flag = false;
        for (int it = 0;it < iter;it++) {
            int syndrome[1023] = { 0 };
            for (int i = 0; i < 1023; i++) {
                for (int j = 0; j < 32; j++) { // 每個s連接 32 個x
                    syndrome[i] ^= x_est[H_s[i][j] - 1] & (H_s_real[i][H_s[i][j] - 1]);   // H_s_real[i][H_s[i][j]-1] == 1
                }
            }

            bool all_zero = true;
            for (int i = 0; i < 1023; i++) {
                if (syndrome[i] != 0) {
                    all_zero = false;
                    break;
                }
            }
            if (all_zero) {
                flag = true;
                break;
            }
            // fi
            int fi[1023] = { 0 };
            for (int i = 0; i < 1023; i++) {
                if (syndrome[i]) { // 只計算 Syndrome 不為零的部分
                    for (int j = 0; j < 32; j++) {
                        fi[H_s[i][j] - 1]++;    //每個s連接到的x都要++
                    }
                }
            }

            // filp
            int max_fi = 0;
            for (int j = 0; j < 1023; j++) {
                if (fi[j] > max_fi) {
                    max_fi = fi[j];
                }
            }
            if (max_fi > threshold) {
                for (int i = 0;i < 1023;i++) {
                    if (fi[i] == max_fi) {
                        x_est[i] ^= 1; // flip
                    }
                }
            }
        }

        for (int i = 0;i < 781;i++) {
            if (x_est[i] != c[i])
                err_bit++;
        }
        if (flag == false)  err_block++;

        cout << "block @ " << b + 1 << endl;
        cout << "#error blocks = " << err_block << endl;
        cout << "#error bits = " << err_bit << endl;

        if (err_block == 100) {
            ber = err_bit / ((b + 1) * 781.0);
            bler = err_block / ((b + 1) * 1.0);
            terminate_b = b + 1;
            break;
        }

        // for next u sequence
        for (int i = 0; i < 6; i++) {
            prev_u[i] = u[781 - 6 + i];
        }
    }
}