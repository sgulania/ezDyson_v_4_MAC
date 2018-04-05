#ifndef _dyson_main_
#define _dyson_main_

/*! \file dyson_main.h
\brief  the "main" for dyson job
*/

#include "cklm.h"
#include "eikr.h"
#include "orbital.h"
#include "aobasis.h"
#include "xyzgrid.h"
#include "klmgrid.h"
//#include "radialfns.h"
#include "readwrite.h"
#include "ylm.h"
#include "sph.h"
#include "pad.h"
#include "rotnmatr.h"
#include "anglegrid.h"

#include "simple_xml_parser.h"
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>

//! main dyson code
bool dyson_main(const char* xmlFileName);

#endif
