#ifndef TOYRINGMULTIPLIER_H
#define TOYRINGMULTIPLIER_H

#include <iostream>
#include <NTL/ZZ.h>
#include <vector>

using namespace std;
using namespace NTL;

//static const long logN = 16;
//static const long logQ = 1200;
static const long logN = 3;
static const long logQ = 4;

static const double sigma = 3.2;
static const long h = 64;
//static const long pbnd = round(log(NTL_SP_BOUND)/log(2.)) - 1;
static const long pbnd = 6;
static const long kbar = pbnd + 1;
static const long kbar2 = 2 * kbar;
static const long logNh = (logN - 1);
static const long logQQ = (2 * logQ);
static const long N = (1 << logN);
static const long Nh = (1 << logNh);
static const long M = (N << 1);
static const long nprimes = (2 + logN + 4 * logQ + pbnd - 1) / pbnd;
static const long Nnprimes = (nprimes << logN);
static const long cbnd = (logQQ + NTL_ZZ_NBITS - 1) / NTL_ZZ_NBITS;
static const long bignum = 0xfffffff;
static const ZZ Q = power2_ZZ(logQ);
static const ZZ QQ = power2_ZZ(logQQ);


class ToyRingMultiplier {
  public:
    uint64_t *pVec = new uint64_t[nprimes];
    uint64_t *prVec = new uint64_t[nprimes];
    uint64_t *pInvVec = new uint64_t[nprimes];
    uint64_t **scaledRootPows = new uint64_t *[nprimes];
    uint64_t **scaledRootInvPows = new uint64_t *[nprimes];
    uint64_t *scaledNInv = new uint64_t[nprimes];
    /* _ntl_general_rem_one_struct* red_ss_array[nprimes]; */
    /* mulmod_precon_t* coeffpinv_array[nprimes]; */

    /* ZZ* pProd = new ZZ[nprimes]; */
    /* ZZ* pProdh = new ZZ[nprimes]; */
    /* ZZ** pHat = new ZZ*[nprimes]; */
    /* uint64_t** pHatInvModp = new uint64_t*[nprimes]; */

    ToyRingMultiplier();

    bool primeTest(uint64_t p);

    void NTT(uint64_t *a, long index);

    void INTT(uint64_t *a, long index);

    /* void CRT(uint64_t* rx, ZZ* x, const long np); */

    /* void addNTTAndEqual(uint64_t* ra, uint64_t* rb, const long np); */

    /* void reconstruct(ZZ* x, uint64_t* rx, long np, const ZZ& QQ); */

    /* void mult(ZZ* x, ZZ* a, ZZ* b, long np, const ZZ& QQ); */

    /* void multNTT(ZZ* x, ZZ* a, uint64_t* rb, long np, const ZZ& QQ); */

    /* void multDNTT(ZZ* x, uint64_t* ra, uint64_t* rb, long np, const ZZ& QQ); */

    /* void multAndEqual(ZZ* a, ZZ* b, long np, const ZZ& QQ); */

    /* void multNTTAndEqual(ZZ* a, uint64_t* rb, long np, const ZZ& QQ); */

    /* void square(ZZ* x, ZZ* a, long np, const ZZ& QQ); */

    /* void squareNTT(ZZ* x, uint64_t* ra, long np, const ZZ& QQ); */

    /* void squareAndEqual(ZZ* a, long np, const ZZ& QQ); */

    void mulMod(uint64_t &r, uint64_t a, uint64_t b, uint64_t p);

    /* void mulModBarrett(uint64_t& r, uint64_t a, uint64_t b, uint64_t p, uint64_t pr); */

    void butt(uint64_t &a, uint64_t &b, uint64_t W, uint64_t p, uint64_t pInv);

    void ibutt(uint64_t &a, uint64_t &b, uint64_t W, uint64_t p, uint64_t pInv);

    void idivN(uint64_t &a, uint64_t NScale, uint64_t p, uint64_t pInv);

    uint64_t invMod(uint64_t x, uint64_t p);

    uint64_t powMod(uint64_t x, uint64_t y, uint64_t p);

    uint64_t inv(uint64_t x);

    uint64_t pow(uint64_t x, uint64_t y);

    uint32_t bitReverse(uint32_t x);

    void findPrimeFactors(vector<uint64_t> &s, uint64_t number);

    uint64_t findPrimitiveRoot(uint64_t m);

    uint64_t findMthRootOfUnity(uint64_t M, uint64_t p);

    void printParams();
};

#endif
