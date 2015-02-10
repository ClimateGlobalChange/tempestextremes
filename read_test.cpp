/////////////////////////////////////////
///
///    \file read_test.cpp
//     \author Marielle Pinheiro
///    \version February 3, 2015
///


//#include "CLIVAR_blocks.h"
#include "Exception.h"
#include "NetCDFUtilities.h"
#include "netcdfcpp.h"

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>

int main(){
  NcFile readin("ERA_2013_ptblock.nc", NcFile::ReadOnly);
  //Gather basic information about file
  int nDims = readin.num_dims();
  int nVars = readin.num_vars();
  for (int a=0; a++; a<nDims){
    NcDim *readDim = readin.get_dim(a);
    std::string dimName = readDim->name();
    long dimSize = readDim->size();
    std::cout<< "The name of dimension "<< a << " is " << dimName.c_str()<<". It is " << dimSize << " long.";
  }
  for (int b=0; b++; b<nVars){
    NcVar *readVar = readin.get_var(b);
    std::string varName = readVar->name();
    int dimsVar = readVar->num_dims();
    int atts = readVar->num_atts();
    std::cout<< "There are "<< atts <<" attributes.";
    std::cout << "Variable " << varName.c_str() << " has attributes ";
    for (int c=0; c++; c<atts){
      NcAtt *varAtt = readVar->get_att(c);
      std::string attName = varAtt->name();
      std::cout << attName.c_str() <<", ";
    }
    std::cout<< "and it has "<< dimsVar << " dimensions."; 
  }
  return 0;
}
