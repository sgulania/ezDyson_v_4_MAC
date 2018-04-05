#ifndef _readwrite_h
#define _readwrite_h

#include "tools.h"
#include "gauss.h"
#include "complexno.h"
#include "rotnmatr.h"
#include "orbital.h"

#include "simple_xml_parser.h"
#include <cstdlib>
#include <iostream>

//! Which MOs in AO basis to read from QChem molden_format output?
void MONumbers(int *monumbers, const char* xmlFileName, int nomos);

//! --
void ReadOrbitalFromFile(Orbital& orb, int nb, const char* xmlFileName);
//! Read MO to plot from qchem.
void ReadMOFromFile(Orbital &mo, int nb, int monumber, const char* xmlFileName);
//!Read info about averaging over molec orientations in lab frame.
MolAvg ReadMolecOrientAvg(const char* xmlFileName);
//!x,y,z lab frame components of ionization laser polarization and IE.
void ReadRIonz(double &ie, double *rioniz, const char* xmlFileName);


void PrintOrbitalsScanZ(Orbital *mostoplot, int nomos, XYZGrid& labgrid);


//!x,y,z lab frame components of excitation laser polarization.
//VM: void ReadRExc(double *rex, FILE *input);
//!Molecular frame components of excitation dipole moment(x,y,z-Left,x,y,z-Right).
//VM: void ReadMExc(double *mex, FILE *input);

#endif


