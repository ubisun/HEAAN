/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include <iostream>
#include <string>
#include "ToyRingMultiplier.h"

using namespace std;

#define SIZE 8

void printArray(string name, uint64_t* buf, int size);
void basic_test();

int main(int argc, char **argv) {

    basic_test();
    return 0;
}


void printArray(string name, uint64_t* buf, int size) {
    cout << name << " : (";
    for (int i = 0; i < size; i++) {
        cout << " " << buf[i];
    }
    cout << " )" << endl << endl;
}

void basic_test() {
    ToyRingMultiplier trm;
    uint64_t X[SIZE] = {4, 1, 4, 12, 1, 13, 8, 9};
    //uint64_t X[SIZE] = {4, 1, 4, 7};

    cout << endl << endl;
    printArray("X", X, SIZE);
    trm.NTT(X, 0);
    printArray("NTT(X)", X, SIZE);
    trm.INTT(X, 0);
    printArray("INTT(X)", X, SIZE);

}

