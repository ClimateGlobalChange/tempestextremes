/////////////////////////////////////////
///
///    \file CLIVAR_block.h
//     \author Marielle Pinheiro
///    \version February 3, 2015
///

#ifndef _CLIVAR_BLOCK_H_
#define _CLIVAR_BLOCK_H_


NcVar* interpolate_lev(
  NcVar* var, 
  NcVar* hyam, 
  NcVar* hybm, 
  NcVar* ps 
);


#endif
