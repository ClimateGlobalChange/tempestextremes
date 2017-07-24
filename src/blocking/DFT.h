////////////////////////////////////////////////
///
///      DFT.h
///      Author: Marielle Pinheiro
///      Version 1.0 May 8th, 2017
///
/////////////////////////////////////////////////

#ifndef _DFT_H_
#define _DFT_H_

#include <cstdlib>
#include <cmath>
#include <complex>
#include <vector>
//#include "/Users/mariellep/tempestextremes/src/base/DataVector.h"
//#include "/Users/mariellep/tempestextremes/src/base/Exception.h"
//#include "DataVector.h"
#include "Exception.h"

std::vector<std::complex<double> > DFT(std::vector<double> inputVals,
         int numCoefs
);

std::vector <double> IDFT(std::vector<std::complex<double> > FFTvals);

#endif
