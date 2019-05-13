/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#include "TestScheme.h"

#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>

#include "Ciphertext.h"
#include "EvaluatorUtils.h"
#include "Profile.h"
#include "Ring.h"
#include "Scheme.h"
#include "SchemeAlgo.h"
#include "SecretKey.h"
#include "StringUtils.h"
#include "TimeUtils.h"
#include "SerializationUtils.h"

using namespace std;
using namespace NTL;


#define MAX_REPEAT 10000
#define MAX_REPEAT_BT 30
double timeResult[MAX_PROFILE_FUNC] = {0};
int dummy;
int dummy1;

//----------------------------------------------------------------------------------
//   STANDARD TESTS
//----------------------------------------------------------------------------------

int startProfile(int ctl) __attribute__((noinline));
int endProfile(int ctl) __attribute__((noinline));
complex<double> normalAdd(complex<double> a, complex<double> b) __attribute__((noinline));
complex<double> normalMult(complex<double> a, complex<double> b) __attribute__((noinline));

// int startProfile(int ctl) { cout << "startProfile! in TestProgram : " << ctl << endl; }

// int endProfile(int ctl) { cout << "endProfile! in TestProgram : " << ctl << endl; }

// int startProfile(int ctl) { dummy = ctl + 2; }

// int endProfile(int ctl) { dummy = ctl + 3; }

int startProfile(int ctl) { dummy1 = dummy + ctl + 2; }

int endProfile(int ctl) { dummy = dummy1 + ctl + 3; }

complex<double> normalAdd(complex<double> a, complex<double> b) { return a + b; }

complex<double> normalMult(complex<double> a, complex<double> b) { return a * b; }

void initTimeResult() {
    for (int i = 0; i < MAX_PROFILE_FUNC; i++)
        timeResult[i] = 0;
}
void printTimeResult() {
    for (int i = 0; i < MAX_PROFILE_FUNC; i++) {
        if (timeResult[i] == 0)
            continue;

        if (i == NORMAL_ADD) cout << "NORMAL_ADD : " ;
        if (i == HE_ADD) cout << "HE_ADD : ";
        if (i == RING_addAndEqual) cout << "RING_addAndEqual : " ;
        if (i == NORMAL_MULT) cout << "NORMAL_MULT : ";
        if (i == HE_MULT) cout << "HE_MULT : ";
        if (i == RING_CRT) cout << "RING_CRT : ";
        if (i == RING_multDNTT) cout << "RING_multDNTT : ";
        if (i == RING_addNTTAndEqual) cout << "RING_addNTTAndEqual : ";
        if (i == RING_rightShiftAndEqual) cout << "RING_rightShiftAndEqual : ";
        if (i == RING_subAndEqual) cout << "RING_subAndEqual : ";
        if (i == RING_leftRotate) cout << "RING_leftRotate : ";
        if (i == RING_MULT_CRT) cout << "RING_MULT_CRT : ";
        if (i == RING_MULT_NTT) cout << "RING_MULT_NTT : ";
        if (i == RING_MULT_multDNTT) cout << "RING_MULT_multDNTT : ";
        if (i == RING_MULT_INTT) cout << "RING_MULT_INTT : ";
        if (i == RING_MULT_reconstruct) cout << "RING_MULT_reconstruct : ";
        if (i == RING_MULT_addNTTAndEqual) cout << "RING_MULT_addNTTAndEqual : ";
        if (i == RING_MULT_multNTT) cout << "RING_MULT_multNTT : ";
        if (i == RING_MULT_squareNTT) cout << "RING_MULT_squareNTT : ";
        if (i == BOOT_STRAP) cout << "BOOT_STRAP : ";
        if (i == BOOT_STRAP_SUB_SUM) cout << "BOOT_STRAP_SUB_SUM : ";
        if (i == BOOT_STRAP_CoeffToSlot) cout << "BOOT_STRAP_CoeffToSlot : ";
        if (i == BOOT_STRAP_EvalExp) cout << "BOOT_STRAP_EvalExp : ";
        if (i == BOOT_STRAP_SlotToCoeff) cout << "BOOT_STRAP_SlotToCoeff : ";

        cout << timeResult[i] << endl;
    }
}

void printFeature() {
#ifdef _BS_PROFILE_
    cout << "_BS_PROFILE_ Feature Yes!" << endl;
#else
    cout << "_BS_PROFILE_ Feature No!" << endl;
#endif
#ifdef _PROFILE_
    cout << "_PROFILE_ Feature Yes!" << endl;
#else
    cout << "_PROFILE_ Feature No!" << endl;
#endif
#ifdef _BS_TIME_PERFORM_
    cout << "_BS_TIME_PERFORM_ Feature Yes!" << endl;
#else
    cout << "_BS_TIME_PERFORM_ Feature No!" << endl;
#endif
#ifdef _TIME_PERFORM_
    cout << "_TIME_PERFORM_ Feature Yes!" << endl;
#else
    cout << "_TIME_PERFORM_ Feature No!" << endl;
#endif


}

void TestScheme::testProfileIns(long logq, long logp, long logn, char *cmd) {
    cout << "!!! START TEST Profile Ins !!!(1)" << endl;
    printFeature();

    srand(time(NULL));
    SetNumThreads(1);
    TimeUtils timeutils;
    long n = (1 << logn);
    complex<double> *base = EvaluatorUtils::randomComplexArray(n);
    complex<double> *operand = EvaluatorUtils::randomComplexArray(n);
    Ciphertext obase, cbase, coperand;

    if (!strcmp(cmd, "NORMAL_ADD")) {
#ifdef _TIME_PERFORM_
        initTimeResult();
        timeutils.start(__func__);
#endif
#ifdef _PROFILE_
        startProfile(NORMAL_ADD);
#endif
        for (long i = 0; i < n; i++) {
            normalAdd(base[i], operand[i]);
        }
#ifdef _PROFILE_
        endProfile(NEGATIVE(NORMAL_ADD));
#endif
#ifdef _TIME_PERFORM_
        timeResult[NORMAL_ADD] += timeutils.stop(__func__);
#endif
    }

    if (!strcmp(cmd, "NORMAL_MULT")) {
#ifdef _TIME_PERFORM_
        initTimeResult();
        timeutils.start(__func__);
#endif
#ifdef _PROFILE_
       startProfile(NORMAL_MULT);
#endif
        for (long i = 0; i < n; i++) {
            normalMult(base[i], operand[i]);
        }

#ifdef _PROFILE_
        endProfile(NEGATIVE(NORMAL_MULT));
#endif
#ifdef _TIME_PERFORM_
        timeResult[NORMAL_MULT] += timeutils.stop(__func__);
#endif
    }

    if (!strcmp(cmd, "HE_ADD")) {
        Ring ring;
        SecretKey secretKey(ring);
        Scheme scheme(secretKey, ring);

        scheme.encrypt(obase, base, n, logp, logq);
        scheme.encrypt(cbase, base, n, logp, logq);
        scheme.encrypt(coperand, operand, n, logp, logq);

        cout << "!!! PREPARED Profile Ins !!!" << endl;
#ifdef _TIME_PERFORM_
        initTimeResult();
        timeutils.start(__func__);
#endif

#ifdef _PROFILE_
        startProfile(HE_ADD);
#endif
        scheme.addAndEqual(cbase, coperand);

#ifdef _PROFILE_
        endProfile(NEGATIVE(HE_ADD));
#endif
#ifdef _TIME_PERFORM_
        timeResult[HE_ADD] += timeutils.stop(__func__);
#endif

    }

    if (!strcmp(cmd, "HE_MULT")) {
        Ring ring;
        SecretKey secretKey(ring);
        Scheme scheme(secretKey, ring);

        scheme.encrypt(obase, base, n, logp, logq);
        scheme.encrypt(cbase, base, n, logp, logq);
        scheme.encrypt(coperand, operand, n, logp, logq);

        cout << "!!! PREPARED Profile Ins !!!" << endl;

#ifdef _TIME_PERFORM_
        initTimeResult();
        timeutils.start(__func__);
#endif
#ifdef _PROFILE_
        startProfile(HE_MULT);
#endif
        scheme.multAndEqual(obase, coperand);

#ifdef _PROFILE_
        endProfile(NEGATIVE(HE_MULT));
#endif
#ifdef _TIME_PERFORM_
        timeResult[HE_MULT] += timeutils.stop(__func__);
#endif

    }

    cout << "!!! END TEST Profile Ins !!!" << endl << endl;

#ifdef _TIME_PERFORM_
    printTimeResult();
#endif
}

void TestScheme::testProfileTime(long logq, long logp, long logn, int repeat) {
    cout << "!!! START TEST Profile !!!(1)" << endl;
    printFeature();

    srand(time(NULL));
    SetNumThreads(1);

    if (repeat > MAX_REPEAT) {
        cout << "repeat count is too large" << endl;
        return;
    } else {
        cout << "repeat : " << repeat << endl;
    }

    TimeUtils timeutils;
    Ring ring;
    SecretKey secretKey(ring);
    Scheme scheme(secretKey, ring);

    long n = (1 << logn);
    complex<double> addR[n];
    complex<double> multR[n];
    complex<double> *base = EvaluatorUtils::randomComplexArray(n);
    complex<double> *operand = EvaluatorUtils::randomComplexArray(n);
    Ciphertext obase, cbase, coperand;

    for (int i = 0; i < n; i++) {
        addR[i] = base[i];
        multR[i] = base[i];
    }

    // base.real(1.0);
    // base.imag(0.0);
    // operand.real(2.0);
    // operand.imag(0.0);

    scheme.encrypt(obase, base, n, logp, logq);
    scheme.encrypt(cbase, base, n, logp, logq);
    scheme.encrypt(coperand, operand, n, logp, logq);

    timeutils.start("Normal Add");
    for (int i = 0; i < MAX_REPEAT; i++) {
        for (int j = 0; j < n; j++) {
            addR[j] += operand[j];
        }
    }
    timeutils.stop("Normal Add");

    timeutils.start("HE Add");
    scheme.addAndEqual(cbase, coperand);
    for (int i = 0; i < MAX_REPEAT; i++) {
        scheme.addAndEqual(cbase, coperand);
    }
    timeutils.stop("HE Add");

    complex<double>* dmult = scheme.decrypt(secretKey, cbase);
    StringUtils::compare(addR, dmult, n, "mult");

    timeutils.start("Normal Mult");
    for (int i = 0; i < repeat; i++) {
        for (int j = 0; j < n; j++) {
            multR[j] *= operand[j];
        }
    }
    timeutils.stop("Normal Mult");

    timeutils.start("HE Mult");
    for (int i = 0; i < repeat; i++) {
        scheme.multAndEqual(obase, coperand);
    }
    timeutils.stop("HE Mult");

    complex<double> *dmult2 = scheme.decrypt(secretKey, obase);
    StringUtils::compare(multR, dmult2, n, "mult");

    cout << "!!! END TEST Profile !!!" << endl;
}

void TestScheme::testTrace(long logq, long logp, long logn) {
    cout << "!!! START TEST TRACE !!!" << endl;
    srand(time(NULL));
    SetNumThreads(8);
    TimeUtils timeutils;
    Ring ring;
    SecretKey secretKey(ring);
    Scheme scheme(secretKey, ring);

    long n = (1 << logn);
    complex<double> *mvec1 = EvaluatorUtils::randomComplexArray(n);
    complex<double> *mvec2 = EvaluatorUtils::randomComplexArray(n);

#ifdef _TRACE_
    // Set mvec with simple value
    mvec1[0].real(1.0);
    mvec1[0].imag(0);
    mvec1[1].real(2.0);
    mvec1[1].imag(0);
    mvec2[0].real(3.0);
    mvec2[0].imag(0);
    mvec2[1].real(4.0);
    mvec2[1].imag(0);
#endif

    Ciphertext cipher1, cipher2;
    scheme.encrypt(cipher1, mvec1, n, logp, logq);
    //scheme.encrypt(cipher2, mvec2, n, logp, logq);

    complex<double> *dvec = scheme.decrypt(secretKey, cipher1);

    StringUtils::compare(mvec1, dvec, n, "trace");
    cout << "!!! END TEST TRACE !!!" << endl;
}


void TestScheme::testEncrypt(long logq, long logp, long logn) {
	cout << "!!! START TEST ENCRYPT !!!" << endl;
	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long n = (1 << logn);
	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	Ciphertext cipher;

	timeutils.start("Encrypt");
	scheme.encrypt(cipher, mvec, n, logp, logq);
	timeutils.stop("Encrypt");

	timeutils.start("Decrypt");
	complex<double>* dvec = scheme.decrypt(secretKey, cipher);
	timeutils.stop("Decrypt");

	StringUtils::compare(mvec, dvec, n, "val");

	cout << "!!! END TEST ENCRYPT !!!" << endl;
}

void TestScheme::testEncryptSingle(long logq, long logp) {
	cout << "!!! START TEST ENCRYPT SINGLE !!!" << endl;
	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	complex<double> mval = EvaluatorUtils::randomComplex();
	Ciphertext cipher;

	timeutils.start("Encrypt Single");
	scheme.encryptSingle(cipher, mval, logp, logq);
	timeutils.stop("Encrypt Single");

	complex<double> dval = scheme.decryptSingle(secretKey, cipher);

	StringUtils::compare(mval, dval, "val");

	cout << "!!! END TEST ENCRYPT SINGLE !!!" << endl;
}

void TestScheme::testAdd(long logq, long logp, long logn) {
	cout << "!!! START TEST ADD !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long n = (1 << logn);
	complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(n);
	complex<double>* madd = new complex<double>[n];

	for(long i = 0; i < n; i++) {
		madd[i] = mvec1[i] + mvec2[i];
	}

	Ciphertext cipher1, cipher2;
	scheme.encrypt(cipher1, mvec1, n, logp, logq);
	scheme.encrypt(cipher2, mvec2, n, logp, logq);

	timeutils.start("Addition");
	scheme.addAndEqual(cipher1, cipher2);
	timeutils.stop("Addition");

	complex<double>* dadd = scheme.decrypt(secretKey, cipher1);

	StringUtils::compare(madd, dadd, n, "add");

	cout << "!!! END TEST ADD !!!" << endl;
}

void TestScheme::testMult(long logq, long logp, long logn) {
	cout << "!!! START TEST MULT !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long n = (1 << logn);
	complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mmult = new complex<double>[n];
	for(long i = 0; i < n; i++) {
		mmult[i] = mvec1[i] * mvec2[i];
	}

	Ciphertext cipher1, cipher2;
	scheme.encrypt(cipher1, mvec1, n, logp, logq);
	scheme.encrypt(cipher2, mvec2, n, logp, logq);

	timeutils.start("Multiplication");
	scheme.multAndEqual(cipher1, cipher2);
	timeutils.stop("Multiplication");

	complex<double>* dmult = scheme.decrypt(secretKey, cipher1);

	StringUtils::compare(mmult, dmult, n, "mult");

	cout << "!!! END TEST MULT !!!" << endl;
}

void TestScheme::testiMult(long logq, long logp, long logn) {
	cout << "!!! START TEST i MULTIPLICATION !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long n = (1 << logn);

	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	complex<double>* imvec = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		imvec[i].real(-mvec[i].imag());
		imvec[i].imag(mvec[i].real());
	}

	Ciphertext cipher;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start("Multiplication by i");
	scheme.imultAndEqual(cipher);
	timeutils.stop("Multiplication by i");

	complex<double>* idvec = scheme.decrypt(secretKey, cipher);

	StringUtils::compare(imvec, idvec, n, "imult");

	cout << "!!! END TEST i MULTIPLICATION !!!" << endl;
}


//----------------------------------------------------------------------------------
//   ROTATE & CONJUGATE
//----------------------------------------------------------------------------------


void TestScheme::testRotateFast(long logq, long logp, long logn, long logr) {
	cout << "!!! START TEST ROTATE FAST !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	long n = (1 << logn);
	long r = (1 << logr);
	scheme.addLeftRotKey(secretKey, r);
	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	Ciphertext cipher;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start("Left Rotate Fast");
	scheme.leftRotateFastAndEqual(cipher, r);
	timeutils.stop("Left Rotate Fast");

	complex<double>* dvec = scheme.decrypt(secretKey, cipher);

	EvaluatorUtils::leftRotateAndEqual(mvec, n, r);
	StringUtils::compare(mvec, dvec, n, "rot");

	cout << "!!! END TEST ROTATE BY POWER OF 2 BATCH !!!" << endl;
}

void TestScheme::testConjugate(long logq, long logp, long logn) {
	cout << "!!! START TEST CONJUGATE !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	scheme.addConjKey(secretKey);

	long n = (1 << logn);

	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mvecconj = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mvecconj[i] = conj(mvec[i]);
	}

	Ciphertext cipher;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start("Conjugate");
	scheme.conjugateAndEqual(cipher);
	timeutils.stop("Conjugate");

	complex<double>* dvecconj = scheme.decrypt(secretKey, cipher);
	StringUtils::compare(mvecconj, dvecconj, n, "conj");

	cout << "!!! END TEST CONJUGATE !!!" << endl;
}


//----------------------------------------------------------------------------------
//   POWER & PRODUCT TESTS
//----------------------------------------------------------------------------------


void TestScheme::testPowerOf2(long logq, long logp, long logn, long logdeg) {
	cout << "!!! START TEST POWER OF 2 !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n = 1 << logn;
	long degree = 1 << logdeg;
	complex<double>* mvec = new complex<double>[n];
	complex<double>* mpow = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mvec[i] = EvaluatorUtils::randomCircle();
		mpow[i] = pow(mvec[i], degree);
	}

	Ciphertext cipher, cpow;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start("Power of 2");
	algo.powerOf2(cpow, cipher, logp, logdeg);
	timeutils.stop("Power of 2");

	complex<double>* dpow = scheme.decrypt(secretKey, cpow);
	StringUtils::compare(mpow, dpow, n, "pow2");

	cout << "!!! END TEST POWER OF 2 !!!" << endl;
}

//-----------------------------------------

void TestScheme::testPower(long logq, long logp, long logn, long degree) {
	cout << "!!! START TEST POWER !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n = 1 << logn;
	complex<double>* mvec = EvaluatorUtils::randomCircleArray(n);
	complex<double>* mpow = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mpow[i] = pow(mvec[i], degree);
	}

	Ciphertext cipher, cpow;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start("Power");
	algo.power(cpow, cipher, logp, degree);
	timeutils.stop("Power");

	complex<double>* dpow = scheme.decrypt(secretKey, cpow);
	StringUtils::compare(mpow, dpow, n, "pow");

	cout << "!!! END TEST POWER !!!" << endl;
}


//----------------------------------------------------------------------------------
//   FUNCTION TESTS
//----------------------------------------------------------------------------------


void TestScheme::testInverse(long logq, long logp, long logn, long steps) {
	cout << "!!! START TEST INVERSE !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n = 1 << logn;
	complex<double>* mvec = EvaluatorUtils::randomCircleArray(n, 0.1);
	complex<double>* minv = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		minv[i] = 1. / mvec[i];
	}

	Ciphertext cipher, cinv;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start("Inverse");
	algo.inverse(cinv, cipher, logp, steps);
	timeutils.stop("Inverse");

	complex<double>* dinv = scheme.decrypt(secretKey, cinv);
	StringUtils::compare(minv, dinv, n, "inv");

	cout << "!!! END TEST INVERSE !!!" << endl;
}

void TestScheme::testLogarithm(long logq, long logp, long logn, long degree) {
	cout << "!!! START TEST LOGARITHM !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n = 1 << logn;
	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n, 0.1);
	complex<double>* mlog = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mlog[i] = log(mvec[i] + 1.);
	}

	Ciphertext cipher, clog;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start(LOGARITHM);
	algo.function(clog, cipher, LOGARITHM, logp, degree);
	timeutils.stop(LOGARITHM);

	complex<double>* dlog = scheme.decrypt(secretKey, clog);
	StringUtils::compare(mlog, dlog, n, LOGARITHM);

	cout << "!!! END TEST LOGARITHM !!!" << endl;
}

void TestScheme::testExponent(long logq, long logp, long logn, long degree) {
	cout << "!!! START TEST EXPONENT !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n = 1 << logn;
	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mexp = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mexp[i] = exp(mvec[i]);
	}

	Ciphertext cipher, cexp;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start(EXPONENT);
	algo.function(cexp, cipher, EXPONENT, logp, degree);
	timeutils.stop(EXPONENT);

	complex<double>* dexp = scheme.decrypt(secretKey, cexp);
	StringUtils::compare(mexp, dexp, n, EXPONENT);

	cout << "!!! END TEST EXPONENT !!!" << endl;
}

void TestScheme::testExponentLazy(long logq, long logp, long logn, long degree) {
	cout << "!!! START TEST EXPONENT LAZY !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n = 1 << logn;
	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	complex<double>* mexp = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		mexp[i] = exp(mvec[i]);
	}
	Ciphertext cipher, cexp;
	scheme.encrypt(cipher, mvec, n, logp, logQ);

	timeutils.start(EXPONENT + " lazy");
	algo.functionLazy(cexp, cipher, EXPONENT, logp, degree);
	timeutils.stop(EXPONENT + " lazy");

	complex<double>* dexp = scheme.decrypt(secretKey, cexp);
	StringUtils::compare(mexp, dexp, n, EXPONENT);

	cout << "!!! END TEST EXPONENT LAZY !!!" << endl;
}

//-----------------------------------------

void TestScheme::testSigmoid(long logq, long logp, long logn, long degree) {
	cout << "!!! START TEST SIGMOID !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n = 1 << logn;

	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	complex<double>* msig = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		msig[i] = exp(mvec[i]) / (1. + exp(mvec[i]));
	}

	Ciphertext cipher, csig;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start(SIGMOID);
	algo.function(csig, cipher, SIGMOID, logp, degree);
	timeutils.stop(SIGMOID);

	complex<double>* dsig = scheme.decrypt(secretKey, csig);
	StringUtils::compare(msig, dsig, n, SIGMOID);

	cout << "!!! END TEST SIGMOID !!!" << endl;
}

void TestScheme::testSigmoidLazy(long logq, long logp, long logn, long degree) {
	cout << "!!! START TEST SIGMOID LAZY !!!" << endl;

	srand(time(NULL));
//	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	SchemeAlgo algo(scheme);

	long n = 1 << logn;
	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	complex<double>* msig = new complex<double>[n];
	for (long i = 0; i < n; ++i) {
		msig[i] = exp(mvec[i]) / (1. + exp(mvec[i]));
	}

	Ciphertext cipher, csig;
	scheme.encrypt(cipher, mvec, n, logp, logq);

	timeutils.start(SIGMOID + " lazy");
	algo.functionLazy(csig, cipher, SIGMOID, logp, degree);
	timeutils.stop(SIGMOID + " lazy");

	complex<double>* dsig = scheme.decrypt(secretKey, csig);
	StringUtils::compare(msig, dsig, n, SIGMOID);

	cout << "!!! END TEST SIGMOID LAZY !!!" << endl;
}


void TestScheme::testWriteAndRead(long logq, long logp, long logSlots) {
	cout << "!!! START TEST WRITE AND READ !!!" << endl;

	cout << "!!! END TEST WRITE AND READ !!!" << endl;
}


void TestScheme::testBootstrap(long logq, long logp, long logSlots, long logT) {
	cout << "!!! START TEST BOOTSTRAP !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	timeutils.start("Key generating");
	scheme.addBootKey(secretKey, logSlots, logq + 4);
	timeutils.stop("Key generated");

	long slots = (1 << logSlots);
	complex<double>* mvec = EvaluatorUtils::randomComplexArray(slots);

	Ciphertext cipher;
	scheme.encrypt(cipher, mvec, slots, logp, logq);

	cout << "cipher logq before: " << cipher.logq << endl;

	scheme.modDownToAndEqual(cipher, logq);
	scheme.normalizeAndEqual(cipher);
	cipher.logq = logQ;
	cipher.logp = logq + 4;

	Ciphertext rot;
	timeutils.start("SubSum");
	for (long i = logSlots; i < logNh; ++i) {
		scheme.leftRotateFast(rot, cipher, (1 << i));
		scheme.addAndEqual(cipher, rot);
	}
	scheme.divByPo2AndEqual(cipher, logNh);
	timeutils.stop("SubSum");

	timeutils.start("CoeffToSlot");
	scheme.coeffToSlotAndEqual(cipher);
	timeutils.stop("CoeffToSlot");

	timeutils.start("EvalExp");
	scheme.evalExpAndEqual(cipher, logT);
	timeutils.stop("EvalExp");

	timeutils.start("SlotToCoeff");
	scheme.slotToCoeffAndEqual(cipher);
	timeutils.stop("SlotToCoeff");

	cipher.logp = logp;
	cout << "cipher logq after: " << cipher.logq << endl;

	complex<double>* dvec = scheme.decrypt(secretKey, cipher);

	StringUtils::compare(mvec, dvec, slots, "boot");

	cout << "!!! END TEST BOOTSRTAP !!!" << endl;
}

void TestScheme::testBootstrapProfile(long logq, long logp, long logSlots, long logT) {
    cout << "!!! START TEST BOOTSTRAP Profile(1) !!!" << endl;
    printFeature();

    srand(time(NULL));
    SetNumThreads(1);
    TimeUtils timeutils;
    Ring ring;
    SecretKey secretKey(ring);
    Scheme scheme(secretKey, ring);

    timeutils.start("Key generating");
    scheme.addBootKey(secretKey, logSlots, logq + 4);
    timeutils.stop("Key generated");

    long slots = (1 << logSlots);
    complex<double> *mvec = EvaluatorUtils::randomComplexArray(slots);

    Ciphertext cipher;
    scheme.encrypt(cipher, mvec, slots, logp, logq);

    cout << "cipher logq before: " << cipher.logq << endl;

#ifdef _BS_TIME_PERFORM_
    initTimeResult();
    timeutils.start(__func__);
#endif

#ifdef _BS_PROFILE_
    cout << "Start BS profile" << endl;
    startProfile(BOOT_STRAP);
#endif
    scheme.modDownToAndEqual(cipher, logq);
    scheme.normalizeAndEqual(cipher);
    cipher.logq = logQ;
    cipher.logp = logq + 4;

    Ciphertext rot;

#ifdef _BS_TIME_PERFORM_
        timeutils.startMid(__func__);
#endif

#ifdef _BS_PROFILE_
        startProfile(BOOT_STRAP_SUB_SUM);
#endif

    for (long i = logSlots; i < logNh; ++i) {
        scheme.leftRotateFast(rot, cipher, (1 << i));
        scheme.addAndEqual(cipher, rot);
    }
    scheme.divByPo2AndEqual(cipher, logNh);

#ifdef _BS_PROFILE_
        endProfile(NEGATIVE(BOOT_STRAP_SUB_SUM));
#endif
#ifdef _BS_TIME_PERFORM_
        timeResult[BOOT_STRAP_SUB_SUM] += timeutils.stopMid(__func__);
#endif

#ifdef _BS_TIME_PERFORM_
        timeutils.startMid(__func__);
#endif

#ifdef _BS_PROFILE_
        startProfile(BOOT_STRAP_CoeffToSlot);
#endif

        scheme.coeffToSlotAndEqual(cipher);

#ifdef _BS_PROFILE_
        endProfile(NEGATIVE(BOOT_STRAP_CoeffToSlot));
#endif
#ifdef _BS_TIME_PERFORM_
        timeResult[BOOT_STRAP_CoeffToSlot] += timeutils.stopMid(__func__);
#endif



#ifdef _BS_TIME_PERFORM_
        timeutils.startMid(__func__);
#endif

#ifdef _BS_PROFILE_
        startProfile(BOOT_STRAP_EvalExp);
#endif

    scheme.evalExpAndEqual(cipher, logT);

#ifdef _BS_PROFILE_
        endProfile(NEGATIVE(BOOT_STRAP_EvalExp));
#endif
#ifdef _BS_TIME_PERFORM_
        timeResult[BOOT_STRAP_EvalExp] += timeutils.stopMid(__func__);
#endif

#ifdef _BS_TIME_PERFORM_
        timeutils.startMid(__func__);
#endif

#ifdef _BS_PROFILE_
        startProfile(BOOT_STRAP_SlotToCoeff);
#endif

    scheme.slotToCoeffAndEqual(cipher);

#ifdef _BS_PROFILE_
        endProfile(NEGATIVE(BOOT_STRAP_SlotToCoeff));
#endif
#ifdef _BS_TIME_PERFORM_
        timeResult[BOOT_STRAP_SlotToCoeff] += timeutils.stopMid(__func__);
#endif


    cipher.logp = logp;

#ifdef _BS_PROFILE_
    endProfile(NEGATIVE(BOOT_STRAP));
#endif
#ifdef _BS_TIME_PERFORM_
    timeResult[BOOT_STRAP] += timeutils.stop(__func__);
#endif

#ifdef _BS_TIME_PERFORM_
    printTimeResult();
#endif


    cout << "cipher logq after: " << cipher.logq << endl;

    complex<double> *dvec = scheme.decrypt(secretKey, cipher);

    StringUtils::compare(mvec, dvec, slots, "boot");

    cout << "!!! END TEST BOOTSRTAP Profile !!!" << endl;
}

void TestScheme::testBootstrapSingleReal(long logq, long logp, long logT) {
	cout << "!!! START TEST BOOTSTRAP SINGLE REAL !!!" << endl;

	srand(time(NULL));
	SetNumThreads(8);
	TimeUtils timeutils;
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);

	timeutils.start("Key generating");
	scheme.addBootKey(secretKey, 0, logq + 4);
	timeutils.stop("Key generated");

	double mval = EvaluatorUtils::randomReal();

	Ciphertext cipher;
	scheme.encryptSingle(cipher, mval, logp, logq);

	cout << "cipher logq before: " << cipher.logq << endl;
	scheme.modDownToAndEqual(cipher, logq);
	scheme.normalizeAndEqual(cipher);
	cipher.logq = logQ;

	Ciphertext rot, cconj;
	timeutils.start("SubSum");
	for (long i = 0; i < logNh; ++i) {
		scheme.leftRotateFast(rot, cipher, 1 << i);
		scheme.addAndEqual(cipher, rot);
	}
	scheme.conjugate(cconj, cipher);
	scheme.addAndEqual(cipher, cconj);
	scheme.divByPo2AndEqual(cipher, logN);
	timeutils.stop("SubSum");

	timeutils.start("EvalExp");
	scheme.evalExpAndEqual(cipher, logT);
	timeutils.stop("EvalExp");

	cout << "cipher logq after: " << cipher.logq << endl;

	cipher.logp = logp;
	complex<double> dval = scheme.decryptSingle(secretKey, cipher);

	StringUtils::compare(mval, dval.real(), "boot");

	cout << "!!! END TEST BOOTSRTAP SINGLE REAL !!!" << endl;
}
