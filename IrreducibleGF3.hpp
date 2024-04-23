// sobolGF3.hpp 		modification for GF(3)
#pragma once

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cstdint>

using namespace std;

//--------------------------------- constants

#define SOBOLGFNSEQLENGTH	10	// 3^10 = 59049 		3^20 = 3486784401

#define N_SOBOL_INIT_TAB_ENTRIES    48 	// IN initIrreducibleGF3.dat

//--------------------------------- global variables
typedef int32_t int32;

int32 sobol_dj[N_SOBOL_INIT_TAB_ENTRIES];
int32 sobol_sj[N_SOBOL_INIT_TAB_ENTRIES];
int32 sobol_aj[N_SOBOL_INIT_TAB_ENTRIES];
int32 sobol_mk[N_SOBOL_INIT_TAB_ENTRIES][32];

//--------------------------------- structures
typedef uint8_t uchar;

unsigned int pow3tab[21] = {1,3,9,27,81,243,729,2187,6561,19683,59049,177147,531441,1594323,4782969,14348907,43046721,129140163,387420489,1162261467,3486784401};
int32 convertToGF3[3] = {0,2,1};

//--------------------------------- Sobol-related routines
inline vector<int32> IntegerDigits(int32 val, const int base, const int len) {
  vector<int32> digits;
  for (int i = 0; i < len; i++) {
    digits.push_back(val % base);
    val = (val / base);
  }
  return digits;
}

inline int32 FromDigits(const vector<int32>& digits, const int base, const int len) {
  int32 pow = 1, res = 0;
  for (int i = 0; i < len; i++) {
    res += pow * digits[i];
    pow = pow*base;
  }
  return res;
}

inline int32 multiplyByFactorInGFN(const int32 x, const int32 factor, const int base, const int len) {
  vector<int32>  digits = IntegerDigits(x,base,len);
  for (int i = 0; i < len; i++)
    digits[i] = (digits[i] * factor) % base;
  return FromDigits(digits, base, len);
}

inline int32 BitXorGFN(const int base, const vector<int32>& lst, const int len, const int polynomialDegree) {
  int32 digits[SOBOLGFNSEQLENGTH][SOBOLGFNSEQLENGTH];
  for (int i = 0; i < SOBOLGFNSEQLENGTH; i++)
    for (int j = 0; j < SOBOLGFNSEQLENGTH; j++)
      digits[i][j] = 0;
  for (int i = 0; i <= polynomialDegree; i++) {
    vector<int32>  d = IntegerDigits(lst[i],base,len);
    for(int j = 0; j < len; j++)  digits[i][j] = d[j];
  }
  vector<int32> final_digits(SOBOLGFNSEQLENGTH+1,0);
  for (int i = 0; i < len; i++) {
    final_digits[i] = 0;
    for (int j = 0; j <= polynomialDegree; j++) {
      final_digits[i] += digits[j][i];
    }
    final_digits[i] = final_digits[i] % base;
  }
  return FromDigits(final_digits, base, len);
}       //BitXorGFN

inline
void generate_mkGF3(const int32 ipolynomial, const int32 polynomialDegree, int32* msobol, const int base) {	// msobol is suppossed to be pre-filled
  vector<int32> polynomial = IntegerDigits(ipolynomial,base,polynomialDegree+1);
  for (int i = polynomialDegree+1; i <= SOBOLGFNSEQLENGTH; i++) {
    vector<int32> lst;
    lst.push_back(msobol[i-polynomialDegree-1]);
    for (int j = 1; j < polynomialDegree+1; j++) {
      lst.push_back(pow3tab[j] * multiplyByFactorInGFN(msobol[i-j-1], convertToGF3[polynomial[polynomialDegree-j]], base, SOBOLGFNSEQLENGTH));
    }
    msobol[i-1] = BitXorGFN(base, lst, i, polynomialDegree);
  }
}

inline int32 load_mk(const string& filename, const bool dbg_flag) {
  int size = 0, dim_from = 1;
  if (dbg_flag) cout << "Loading J&K file " << filename << " ... dim_from=" << dim_from << endl;
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cout << "load_mk: " << filename << " not found." << "\n";
    exit(1);
  };
  char c = file.get();
  if( c == 'd' ) {
    file.ignore(256, '\n');
  } else {
    file.putback(c);
  }
  c = file.get();
  if( c == 'd' ) {
    file.ignore(256, '\n');
  } else {
    file.putback(c);
  }
  int32 index = dim_from;
  while (file.good() && index < N_SOBOL_INIT_TAB_ENTRIES) {
    int32 d, sj, aj;
    file >> d >> sj >> aj;
    sobol_aj[index] = aj;
    sobol_sj[index] = sj;
    sobol_dj[index] = d;
    if (dbg_flag) std::cout << index << " : " << d << " " << sj << " "  << aj << " \t " ;
    for (int i = 0; i < sj; ++i) {
      file >> sobol_mk[index][i];
      if (dbg_flag) std::cout << sobol_mk[index][i] << " " ;
    }
    if (dbg_flag) std::cout << std::endl;
    
    index++;
  }
  return index-2;
}	// load_mk


//--------------------------------- end of Sobol-related routines
