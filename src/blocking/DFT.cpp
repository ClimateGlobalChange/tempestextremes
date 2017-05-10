/////////////////////////////////////////////////////////
///
///         DFT.cpp
///         Author: Marielle Pinheiro
///         Version 1.0 May 8, 2017
///
/////////////////////////////////////////////////////////

/*
This is a discrete Fourier transform function that will perform 
*/

#include <cstdlib>
#include <cmath>
#include <complex>
#include <vector>
//#include "/Users/mariellep/tempestextremes/src/base/DataVector.h"
//#include "/Users/mariellep/tempestextremes/src/base/Exception.h"
#include "DFT.h"
//#include "DataVector.h"
#include "Exception.h"

std::vector<std::complex<double> > DFT(std::vector<double> inputVals,
         int numCoefs
         ){
  double pi = std::atan(1.0)*4.0;
  std::complex <double> compi(0.,1.);

  int N = inputVals.size();

  if (numCoefs > N){
    _EXCEPTIONT("Number of specified coefficients exceeds length of input vector.");
  }

  //Declare the output array for the Fourier coefficients
  std::vector <std::complex<double> > FourierCoefs(N);

  //Couple of values for calculations
  double Ndiv = 1./N;
  std::complex <double> expCoef(0.,0.);
  std::complex <double> sumVals;
  //Begin calculating the coefficients
  for (int k=0; k<N; k++){
    double fk = float(k);
    expCoef = -2.*compi*pi*fk*Ndiv;
    sumVals = std::complex<double>(0.,0.);
    for (int n=0; n<N; n++){
      double fn = float(n);
      sumVals += inputVals[n]*std::exp(fn*expCoef);
    }
    FourierCoefs[k] = sumVals;
  }
  if (numCoefs<N){
    //Zero out the higher wavenumbers
    for (int i=numCoefs; i<(N-numCoefs); i++){
      FourierCoefs[i] = std::complex<double>(0.,0.);
    }
  }
  return(FourierCoefs);
}

std::vector <double> IDFT(std::vector<std::complex<double> > FFTvals){
  double pi = std::atan(1.)*4.;
  std::complex <double> compi(0.,1.);

  int N = FFTvals.size();
  //Create an output data vector (only the real part will be returned)
  std::vector <std::complex<double> > outputs(N);
  double Ndiv = 1./N;
  std::complex <double> expCoef(0.,0.);
  std::complex <double> sumVals;

  for (int k=0; k<N; k++){
    double fk = float(k);
    expCoef = 2.*compi*pi*fk*Ndiv;
    sumVals = std::complex<double>(0.,0.);
    for (int n=0; n<N; n++){
      double fn = float(n);
      sumVals += FFTvals[n]*std::exp(fn*expCoef);
    }
    outputs[k] = sumVals*Ndiv;
  }
  std::vector<double> realOutputs(N);
  for (int i=0; i<N; i++){
    realOutputs[i] = outputs[i].real();
  }
  return(realOutputs);
}
