#ifndef _AOBasis_h
#define _AOBasis_h

#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "tools.h"
#include "gauss.h"

#include "simple_xml_parser.h"

//!AO basis information
class AOBasis
{
  //! number of AOs
  int nbasis;
  int natoms;
  //! Array [nbasis] of contracted gaussians
  Gauss *basis;
  double *xgeom;    
  double *ygeom;    
  double *zgeom;    
  //! Number of basis functions per atom
  int *gtos_atom;
  int purecart[4];      

  bool IfAlloc() const { return nbasis*natoms; }
  void Alloc(int nb, int na);
  void Free();

  AOBasis() : nbasis(0), natoms(0) {} 
  //  virtual AOBasis& operator=(const AOBasis& other) = 0;

 public:

  ~AOBasis() ;
  
  //VM: void IniAOBasis(FILE *infile);
  void IniAOBasis(const char* xmlFileName);
  //void IniAOBasis(int nb, int na, Gauss *bas, double *xg, double *yg, double *zg, 
  //int *gtos, const int *purcar);
  void PrintAOBasis(FILE *fout=stdout);

  int NAtoms() const { return natoms; }
  int NBasis() const { return nbasis; }
  Gauss& GetGauss(int i) const { return basis[i]; }
  double CalcGaussValue(int ng, int nat, double x, double y, double z) const
    {
      //check(!(ng<nbasis && nat < natoms),"AOBasis::GetGaussValue() : invalid ng/na");
      return basis[ng].GaussValue(x,y,z,xgeom[nat],ygeom[nat],zgeom[nat]);
    }

  double GetGTOsAtom(int i) const { return gtos_atom[i]; }
  //Shift the basis in the new coordinate system 
  void ShiftGeom(double *newcenter);

  friend AOBasis& TheAOBasis();
};

extern  AOBasis& TheAOBasis();


#endif
