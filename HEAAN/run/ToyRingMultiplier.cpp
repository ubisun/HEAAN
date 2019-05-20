#include "ToyRingMultiplier.h"

ToyRingMultiplier::ToyRingMultiplier() {

    printParams();

    uint64_t primetest = (1ULL << pbnd) + 1;
    for (long i = 0; i < nprimes; ++i) {
        while (true) {
            primetest += M;
            if (primeTest(primetest)) {
                pVec[i] = primetest;
                break;
            }
        }
    }

    cout << "===== Pre-computation =====" << endl;
    for (long i = 0; i < nprimes; ++i) {
        //red_ss_array[i] = _ntl_general_rem_one_struct_build(pVec[i]);
        pInvVec[i] = inv(pVec[i]);
        prVec[i] = (static_cast<unsigned __int128>(1) << kbar2) / pVec[i];
        uint64_t root = findMthRootOfUnity(M, pVec[i]);
        uint64_t rootinv = invMod(root, pVec[i]);
        uint64_t NInv = invMod(N, pVec[i]);
        mulMod(scaledNInv[i], NInv, (1ULL << 32), pVec[i]);
        mulMod(scaledNInv[i], scaledNInv[i], (1ULL << 32), pVec[i]);
        scaledRootPows[i] = new uint64_t[N]();
        scaledRootInvPows[i] = new uint64_t[N]();
        uint64_t power = 1;
        uint64_t powerInv = 1;

        cout << endl << "pVec[" << i << "] " << pVec[i];
        cout << "\tpInvVec[" << i << "] " << pInvVec[i];
        cout << "\t\troot[" << i << "] " << root;
        cout << "\trootinv[" << i << "] " << rootinv;
        cout << "\tNInv[" << i << "] " << NInv << endl;

        for (long j = 0; j < N; ++j) {
            uint32_t jprime = bitReverse(static_cast<uint32_t>(j)) >> (32 - logN);
            uint64_t rootpow = power;
            mulMod(scaledRootPows[i][jprime], rootpow, (1ULL << 32), pVec[i]);
            mulMod(scaledRootPows[i][jprime], scaledRootPows[i][jprime], (1ULL << 32), pVec[i]);
            uint64_t rootpowInv = powerInv;
            mulMod(scaledRootInvPows[i][jprime], rootpowInv, (1ULL << 32), pVec[i]);
            mulMod(scaledRootInvPows[i][jprime], scaledRootInvPows[i][jprime], (1ULL << 32),
                   pVec[i]);
            mulMod(power, power, root, pVec[i]);
            mulMod(powerInv, powerInv, rootinv, pVec[i]);
            cout << "  => ";
            cout << "j=" << j;
            cout << "\tjprime=" << jprime;
            cout << "\trootpow=" << rootpow;
            cout << "\trootpowInv=" << rootpowInv;
            cout << "\tscaledRootPows[" << i << "][" << jprime << "]=" << scaledRootPows[i][jprime];
            cout << "   \tscaledRootInvPows[" << i << "][" << jprime << "]="
                 << scaledRootInvPows[i][jprime] << endl;
        }
    }
    cout << endl << endl;


    // for (long i = 0; i < nprimes; ++i) {
    //     coeffpinv_array[i] = new mulmod_precon_t[i + 1];
    //     pProd[i] = (i == 0) ? to_ZZ((long)pVec[i]) : pProd[i - 1] * (long)pVec[i];
    //     pProdh[i] = pProd[i] / 2;
    //     pHat[i] = new ZZ[i + 1];
    //     pHatInvModp[i] = new uint64_t[i + 1];
    //     for (long j = 0; j < i + 1; ++j) {
    //         pHat[i][j] = ZZ(1);
    //         for (long k = 0; k < j; ++k) {
    //             pHat[i][j] *= (long)pVec[k];
    //         }
    //         for (long k = j + 1; k < i + 1; ++k) {
    //             pHat[i][j] *= (long)pVec[k];
    //         }
    //         pHatInvModp[i][j] = to_long(pHat[i][j] % (long)pVec[j]);
    //         pHatInvModp[i][j] = invMod(pHatInvModp[i][j], pVec[j]);
    //         coeffpinv_array[i][j] = PrepMulModPrecon(pHatInvModp[i][j], pVec[j]);
    //     }
    // }
}

bool ToyRingMultiplier::primeTest(uint64_t p) {
    if (p < 2)
        return false;
    if (p != 2 && p % 2 == 0)
        return false;
    uint64_t s = p - 1;
    while (s % 2 == 0) {
        s /= 2;
    }
    for (long i = 0; i < 200; i++) {
        uint64_t temp1 = rand();
        temp1 = (temp1 << 32) | rand();
        temp1 = temp1 % (p - 1) + 1;
        uint64_t temp2 = s;
        uint64_t mod = powMod(temp1, temp2, p);
        while (temp2 != p - 1 && mod != 1 && mod != p - 1) {
            mulMod(mod, mod, mod, p);
            temp2 *= 2;
        }
        if (mod != p - 1 && temp2 % 2 == 0)
            return false;
    }
    return true;
}

void ToyRingMultiplier::NTT(uint64_t *a, long index) {
    long t = N;
    long logt1 = logN + 1;
    uint64_t p = pVec[index];
    uint64_t pInv = pInvVec[index];

    for (long m = 1; m < N; m <<= 1) {
        t >>= 1;
        logt1 -= 1;
        for (long i = 0; i < m; i++) {
            long j1 = i << logt1;
            long j2 = j1 + t - 1;
            uint64_t W = scaledRootPows[index][m + i];
            for (long j = j1; j <= j2; j++) {
                butt(a[j], a[j + t], W, p, pInv);
            }
        }
    }
}

void ToyRingMultiplier::INTT(uint64_t *a, long index) {
    uint64_t p = pVec[index];
    uint64_t pInv = pInvVec[index];
    long t = 1;
    for (long m = N; m > 1; m >>= 1) {
        long j1 = 0;
        long h = m >> 1;
        for (long i = 0; i < h; i++) {
            long j2 = j1 + t - 1;
            uint64_t W = scaledRootInvPows[index][h + i];
            for (long j = j1; j <= j2; j++) {
                ibutt(a[j], a[j + t], W, p, pInv);
            }
            j1 += (t << 1);
        }
        t <<= 1;
    }

    uint64_t NScale = scaledNInv[index];
    for (long i = 0; i < N; i++) {
        idivN(a[i], NScale, p, pInv);
    }
}

void ToyRingMultiplier::mulMod(uint64_t &r, uint64_t a, uint64_t b, uint64_t m) {
    unsigned __int128 mul = static_cast<unsigned __int128>(a) * b;
    mul %= static_cast<unsigned __int128>(m);
    r = static_cast<uint64_t>(mul);
}

void ToyRingMultiplier::butt(uint64_t &a, uint64_t &b, uint64_t W, uint64_t p, uint64_t pInv) {
    unsigned __int128 U = static_cast<unsigned __int128>(b) * W;
    uint64_t U0 = static_cast<uint64_t>(U);
    uint64_t U1 = U >> 64;
    uint64_t Q = U0 * pInv;
    unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * p;
    uint64_t H = Hx >> 64;
    uint64_t V = U1 < H ? U1 + p - H : U1 - H;
    b = a < V ? a + p - V : a - V;
    a += V;
    if (a > p)
        a -= p;
}

void ToyRingMultiplier::ibutt(uint64_t &a, uint64_t &b, uint64_t W, uint64_t p, uint64_t pInv) {
    uint64_t T = a < b ? a + p - b : a - b;
    a += b;
    if (a > p)
        a -= p;
    unsigned __int128 UU = static_cast<unsigned __int128>(T) * W;
    uint64_t U0 = static_cast<uint64_t>(UU);
    uint64_t U1 = UU >> 64;
    uint64_t Q = U0 * pInv;
    unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * p;
    uint64_t H = Hx >> 64;
    b = (U1 < H) ? U1 + p - H : U1 - H;
}

void ToyRingMultiplier::idivN(uint64_t &a, uint64_t NScale, uint64_t p, uint64_t pInv) {
    unsigned __int128 U = static_cast<unsigned __int128>(a) * NScale;
    uint64_t U0 = static_cast<uint64_t>(U);
    uint64_t U1 = U >> 64;
    uint64_t Q = U0 * pInv;
    unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * p;
    uint64_t H = Hx >> 64;
    a = (U1 < H) ? U1 + p - H : U1 - H;
}

uint64_t ToyRingMultiplier::invMod(uint64_t x, uint64_t m) { return powMod(x, m - 2, m); }

uint64_t ToyRingMultiplier::powMod(uint64_t x, uint64_t y, uint64_t modulus) {
    uint64_t res = 1;
    while (y > 0) {
        if (y & 1) {
            mulMod(res, res, x, modulus);
        }
        y = y >> 1;
        mulMod(x, x, x, modulus);
    }
    return res;
}

uint64_t ToyRingMultiplier::inv(uint64_t x) { return pow(x, static_cast<uint64_t>(-1)); }

uint64_t ToyRingMultiplier::pow(uint64_t x, uint64_t y) {
    uint64_t res = 1;
    while (y > 0) {
        if (y & 1) {
            res *= x;
        }
        y = y >> 1;
        x *= x;
    }
    return res;
}

uint32_t ToyRingMultiplier::bitReverse(uint32_t x) {
    x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
    x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
    x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
    x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
    return ((x >> 16) | (x << 16));
}

void ToyRingMultiplier::findPrimeFactors(vector<uint64_t> &s, uint64_t number) {
    while (number % 2 == 0) {
        s.push_back(2);
        number /= 2;
    }
    for (uint64_t i = 3; i < sqrt(number); i++) {
        while (number % i == 0) {
            s.push_back(i);
            number /= i;
        }
    }
    if (number > 2) {
        s.push_back(number);
    }
}

uint64_t ToyRingMultiplier::findPrimitiveRoot(uint64_t modulus) {
    vector<uint64_t> s;
    uint64_t phi = modulus - 1;
    findPrimeFactors(s, phi);
    for (uint64_t r = 2; r <= phi; r++) {
        bool flag = false;
        for (auto it = s.begin(); it != s.end(); it++) {
            if (powMod(r, phi / (*it), modulus) == 1) {
                flag = true;
                break;
            }
        }
        if (flag == false) {
            return r;
        }
    }
    return -1;
}

uint64_t ToyRingMultiplier::findMthRootOfUnity(uint64_t M, uint64_t mod) {
    uint64_t res;
    res = findPrimitiveRoot(mod);
    if ((mod - 1) % M == 0) {
        uint64_t factor = (mod - 1) / M;
        res = powMod(res, factor, mod);
        return res;
    } else {
        return -1;
    }
}

void ToyRingMultiplier::printParams() {

    cout << "===== Params Info =====" << endl;
    cout << "[" << logN << "\t] long logN -----> (Origin 16)" << endl;
    cout << "[" << logQ << "\t] long logQ -----> (Origin 1200)" << endl;
    cout << "[" << sigma << "\t] double sigma = 3.2" << endl;
    cout << "[" << h << "\t] long h = 64" << endl;
    cout << "[" << pbnd << "\t] long pbnd = round(log(NTL_SP_BOUND)/log(2.)) - 1 -----> (Origin 59)"
         << endl;
    cout << "[" << kbar << "\t] long kbar = pbnd + 1" << endl;
    cout << "[" << kbar2 << "\t] long kbar2 = 2 * kbar" << endl;
    cout << "[" << logNh << "\t] long logNh = (logN - 1)" << endl;
    cout << "[" << logQQ << "\t] long logQQ = (2 * logQ)" << endl;
    cout << "[" << N << "\t] long N = (1 << logN)" << endl;
    cout << "[" << Nh << "\t] long Nh = (1 << logNh)" << endl;
    cout << "[" << M << "\t] long M = (N << 1)" << endl;
    cout << "[" << nprimes << "\t] long nprimes = (2 + logN + 4 * logQ + pbnd - 1) / pbnd" << endl;
    cout << "[" << Nnprimes << "\t] long Nnprimes = (nprimes << logN)" << endl;
    cout << "[" << cbnd << "\t] long cbnd = (logQQ + NTL_ZZ_NBITS - 1) / NTL_ZZ_NBITS" << endl;
    cout << "[" << Q << "\t] ZZ Q = power2_ZZ(logQ)" << endl;
    cout << "[" << QQ << "\t] ZZ QQ = power2_ZZ(logQQ)" << endl;
    cout << "[" << bignum << "\t] long bignum = 0xfffffff" << endl;
    cout << endl << endl;
}
