////////////////////////////////////////////////////////
///     \file blockingGHAnom.cpp
///     \author Marielle Pinheiro
///     \version November 30, 2015

/*This code calculates anomalies of 500 mb geopotential heights
based on Dole and Gordon 1983, with 2 modifications: 1) re-using
existing averaging code due to time constraints and 2) defining
the anomaly magnitude as 100m rather than 1 sd (perhaps modify
in the future?)
*/
#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "TimeObj.h"
