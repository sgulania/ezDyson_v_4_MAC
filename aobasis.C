#include "aobasis.h"

void AOBasis::Alloc(int nb, int na) 
{
  check(!(nb>0),"AOBasis::Alloc() : Invalid nbasis\n");
  check(!(na>0),"AOBasis::Alloc() : Invalid natoms\n");

  Free();

  nbasis=nb;
  natoms=na;
  basis=new Gauss[nbasis];
  xgeom=new double[natoms];
  ygeom=new double[natoms];
  zgeom=new double[natoms];
  gtos_atom=new int[natoms];
}

void AOBasis::Free()
{
  if(IfAlloc())
   {
     delete [] basis;
     delete [] xgeom;;
     delete [] ygeom;
     delete [] zgeom;
     delete [] gtos_atom;
     nbasis=0;
     natoms=0;
   }
}

AOBasis::~AOBasis() 
{ 
  Free();
}


void AOBasis::ShiftGeom(double *newcenter)
{
  printf("\nAOBasis::ShiftGeom(): Shifting molecular geometry to the new center\n");
  fflush(stdout);
  for(int i=0; i<natoms; i++)
    {
      xgeom[i]-=newcenter[X];
      ygeom[i]-=newcenter[Y];
      zgeom[i]-=newcenter[Z];
    }
  printf("\nAOBasis::ShiftGeom(): New molecular geometry is:\n");
  printf(" atom         X             Y             Z\n");
  for (int i=0; i<natoms; i++)
    printf("%4d %13lf %13lf %13lf\n",i+1,xgeom[i]/ANGS_TO_BOHR,ygeom[i]/ANGS_TO_BOHR,zgeom[i]/ANGS_TO_BOHR);
}

void AOBasis::PrintAOBasis(FILE *fout)
{
  // Print purecart.//
  fprintf(fout,"purecart %d%d%d%d\n",purecart[0],purecart[1],purecart[2],purecart[3]);

  //Print standard geom.
  fprintf(fout,"\nStandard geometry\n");
  fprintf(fout," atom         X             Y             Z\n");
  for (int i=0; i<natoms; i++)
    fprintf(fout,"%4d %13lf %13lf %13lf\n",i+1,xgeom[i]/ANGS_TO_BOHR,ygeom[i]/ANGS_TO_BOHR,zgeom[i]/ANGS_TO_BOHR);
  
  // Print AO basis information.
  fprintf(fout,"\nAO basis\n");
  fprintf(fout," atom  AO    NctG ctrG  L    ia   ja   ka       exp           coef        coeff_norm\n");
  for (int n=0,k=0,i=0; i<nbasis; i++)
    {
      for (int j=0; j<basis[i].NContr(); j++)
	fprintf(fout,"%5d%5d%5d%5d%5d%5d%5d%5d%15lE%15lE%15lE\n",n+1,i+1,j+1,basis[i].NContr(),basis[i].AngMom(),
		basis[i].Ia(),basis[i].Ja(),basis[i].Ka(),basis[i].Alpha(j),basis[i].Coeff(j),basis[i].CoeffNorm(j));
      basis[i].CheckAngMom();
      k++;
      if (k == gtos_atom[n])
	{
	  k=0;
         n++;
	}
     }   

  // Print no contracted GTO's per atom
  fprintf(fout,"\n atom   no_basis_fns\n");
  for (int n=0; n<natoms; n++)
    fprintf(fout,"%4d %10d\n",n+1,gtos_atom[n]);
  fprintf(fout,"-----------------------------------------------------------------------------------\n");
}

AOBasis& TheAOBasis()
{
  static AOBasis theAOBasis;
  //if(!theAOBasis)
  //theAOBasis=new AOBasis;
  
  //Note there is no check if theAOBasis was initialized
  return theAOBasis;
}


// ======================================================================
//  AOBasis::IniAOBasis()
// ======================================================================

void AOBasis::IniAOBasis(const char* xmlFileName)
{
  static bool if_ini=0;
  check(if_ini,"AOBasis::IniAOBasis() : Already initialized\n");
  if_ini=1;

 
  simpleXMLparser xmlF;
  xmlF.assignFile(xmlFileName); 
  //VM:int na = GetNoAtoms(infile);
  xmlF.reset().node("root").node("geometry");
  int na = xmlF.getIntValue("n_of_atoms");

  xmlF.reset().node("root").node("basis");

  int nb = xmlF.getIntValue("n_of_basis_functions");

  Alloc(nb,na);

  std::string Str;

  //-------------------- GetAtomicXYZ --------------------
  xmlF.reset().node("root").node("geometry");
  My_istringstream geom_iStr( xmlF.value("text") );


  for (int i=0; i<natoms; i++)
    {
      //discard an atomic name:
      geom_iStr.getNextWord(Str);
      //get coordinates & convert to a.u.
      xgeom[i] = geom_iStr.getNextDouble()*ANGS_TO_BOHR;
      ygeom[i] = geom_iStr.getNextDouble()*ANGS_TO_BOHR;
      zgeom[i] = geom_iStr.getNextDouble()*ANGS_TO_BOHR;
      xmlF.exitOnFormatError(geom_iStr.fail());
    }

  //-------------------- ReadPureCart --------------------
  xmlF.reset().node("root").node("basis");
// SG: changed purecart from string to integer (then converted to string). Avoids problem with reading purecart
//  xmlF.getWordValue(Str, "purecart");

  int tempStr = xmlF.getIntValue("purecart");
  std::ostringstream convert;
  convert << tempStr;
  Str = convert.str();

  if (Str.size()>4)
    {
      std::cout << "Wrong format of \"purecart\": " << Str <<". Purecart should be of no more than four digits and each one should be either 1 or 2.\n";
      xmlF.exitOnFormatError(true);
    }

  int n=0; //purecrt index 0..3;
  for (int i=0; i<Str.size(); i++, n++)
    if (Str[i]=='1')
	purecart[i]=1;
    else if (Str[i]=='2')
      purecart[i]=2;
    else
      {
	std::cout << "Wrong format: \"purecart["<<i+1<<"]\"="<<Str[i]<<" in \"" << Str <<"\"\n"
		  <<"Purecart should be of no more than four digits and each one should be either 1 or 2.\n";
	xmlF.exitOnFormatError(true);
      }
 
  //fill the rest of purecart (default length=4) with the defalt values {1,1,1,1};
  for (int i=n; i<4; i++)
    {
      if (i==3) purecart[i]=1;
      else purecart[i]=1;
    }

   std::string AOordering;
   int AOorder = 0;
   xmlF.reset().node("root");
   AOordering=xmlF.node("basis").value("AO_ordering");
   if (AOordering=="Q-Chem")
     AOorder = 1; // Q-Chem AO order
   else if (AOordering=="Molden")
     AOorder = 2; // Molden AO order
   else {
     std::cout << "AO_ordering must be 1 (for Q-Chem) or 2 (for Molden).\n" ;
     xmlF.exitOnFormatError(true);
   }

  //-------------------- ReadAOBasis --------------------
  xmlF.reset().node("root").node("basis");
  //! number of gaussians in a contracted b.fn.
  int nc; 
  //! amplitude
  double sc; 
  //! arrays of exponents and coefficients of the gaussians in each contraction
  double *exps=NULL, *coefs=NULL, *pcoefs=NULL; 
  
  //! total number of AO b.fns. loaded for the last contracted b.fn.
  int gtos_loaded; 
  //! total number of AO b.fns. on the current atom
  int gtos_on_current_atom; 
  //! total number of AO b.fns. (for the whole molecule)
  int gtos_count=0; 

  //! if a basis set on the current atom is over
  bool ifEndOfBasisOnAtom;

  //! for each atom:
  for (int atom_number=0; atom_number<natoms; atom_number++)
    {
      //! read the basis on the n-th atom
      xmlF.node("atom", atom_number+1);
      My_istringstream basis_iStr( xmlF.value("text") );

      //! read and discard the atom name:
      basis_iStr.getNextWord(Str);
 
      //! read and discard the coef.
      basis_iStr.getNextDouble();

      ifEndOfBasisOnAtom=false;
      gtos_on_current_atom=0;

      //! read first orbital type to Str:
      basis_iStr.getNextWord(Str);
 
     //! load b.fns. on the current atom (**** is the end)
      while ( not(ifEndOfBasisOnAtom) and not(basis_iStr.fail()) )
	{
	  //! read the number of contractions:
	  nc = basis_iStr.getNextInt();
	  sc = basis_iStr.getNextDouble();
	  xmlF.exitOnFormatError(basis_iStr.fail());

	  exps=new double[nc]; // exponents
	  coefs=new double[nc]; // coefficients
	  pcoefs=new double[nc]; // coefficients for SP

	  //! load AO b.fns.
	  if (Str=="S")
	    {
	      for (int j=0; j<nc; j++)
		{
		  exps[j] = basis_iStr.getNextDouble();
		  coefs[j]= basis_iStr.getNextDouble();
		}
	      basis[gtos_count].IniGauss(0,nc,0,0,0,0,0,0,0,0,0,exps,coefs,sc);
	      gtos_loaded=1;
	    }

	  else if (Str=="P")
	    {
	      for (int j=0; j<nc; j++)
		{
		  exps[j] = basis_iStr.getNextDouble();
		  coefs[j]= basis_iStr.getNextDouble();
		}

	      basis[gtos_count].IniGauss(1,nc,1,0,0,0,0,0,0,0,0,exps,coefs,sc);
	      basis[gtos_count+1].IniGauss(1,nc,0,1,0,0,0,0,0,0,0,exps,coefs,sc);
	      basis[gtos_count+2].IniGauss(1,nc,0,0,1,0,0,0,0,0,0,exps,coefs,sc);
	      gtos_loaded=3;
	    }

	  else if (Str=="D")
	    {
	      for (int j=0; j<nc; j++)
		{
		  exps[j] = basis_iStr.getNextDouble();
		  coefs[j]= basis_iStr.getNextDouble();
		}

	      if (purecart[3] == 1 && AOorder == 1)
		{
		  basis[gtos_count].IniGauss(2,nc,1,1,0,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+1].IniGauss(2,nc,0,1,1,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+2].IniGauss(2,nc,0,0,0,0,2,0,0,0,0,exps,coefs,sc);                 
		  basis[gtos_count+3].IniGauss(2,nc,1,0,1,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+4].IniGauss(2,nc,0,0,0,2,0,0,0,0,0,exps,coefs,sc);
		  gtos_loaded=5;
		}
	      if (purecart[3] == 2 && AOorder == 1)
		{
		  basis[gtos_count].IniGauss(2,nc,2,0,0,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+1].IniGauss(2,nc,1,1,0,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+2].IniGauss(2,nc,0,2,0,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+3].IniGauss(2,nc,1,0,1,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+4].IniGauss(2,nc,0,1,1,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+5].IniGauss(2,nc,0,0,2,0,0,0,0,0,0,exps,coefs,sc);
		  gtos_loaded=6;
		}
              if (purecart[3] == 1 && AOorder == 2)
                {
                  basis[gtos_count].IniGauss(2,nc,0,0,0,0,2,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+1].IniGauss(2,nc,1,0,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+2].IniGauss(2,nc,0,1,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+3].IniGauss(2,nc,0,0,0,2,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+4].IniGauss(2,nc,1,1,0,0,0,0,0,0,0,exps,coefs,sc);
                  gtos_loaded=5;
                }
              if (purecart[3] == 2 && AOorder == 2)
                {
                  basis[gtos_count].IniGauss(2,nc,2,0,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+1].IniGauss(2,nc,0,2,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+2].IniGauss(2,nc,0,0,2,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+3].IniGauss(2,nc,1,1,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+4].IniGauss(2,nc,1,0,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+5].IniGauss(2,nc,0,1,1,0,0,0,0,0,0,exps,coefs,sc);
                  gtos_loaded=6;
                }

	    }
	
	  else if (Str=="F")
	    {
	      for (int j=0; j<nc; j++)
		{
		  exps[j] = basis_iStr.getNextDouble();
		  coefs[j]= basis_iStr.getNextDouble();
		}

	      if (purecart[2] == 1 && AOorder == 1)
		{
		  basis[gtos_count].IniGauss(3,nc,0,0,0,0,0,3,1,0,0,exps,coefs,sc);
		  basis[gtos_count+1].IniGauss(3,nc,1,1,1,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+2].IniGauss(3,nc,0,0,0,0,0,3,3,0,0,exps,coefs,sc);
		  basis[gtos_count+3].IniGauss(3,nc,0,0,0,0,0,3,4,0,0,exps,coefs,sc);
		  basis[gtos_count+4].IniGauss(3,nc,0,0,0,0,0,3,5,0,0,exps,coefs,sc);
		  basis[gtos_count+5].IniGauss(3,nc,0,0,0,0,0,3,6,0,0,exps,coefs,sc);
		  basis[gtos_count+6].IniGauss(3,nc,0,0,0,0,0,3,7,0,0,exps,coefs,sc);
		  gtos_loaded=7;
		}
	      if (purecart[2] == 2 && AOorder == 1)
		{
		  basis[gtos_count].IniGauss(3,nc,3,0,0,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+1].IniGauss(3,nc,2,1,0,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+2].IniGauss(3,nc,1,2,0,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+3].IniGauss(3,nc,0,3,0,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+4].IniGauss(3,nc,2,0,1,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+5].IniGauss(3,nc,1,1,1,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+6].IniGauss(3,nc,0,2,1,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+7].IniGauss(3,nc,1,0,2,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+8].IniGauss(3,nc,0,1,2,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+9].IniGauss(3,nc,0,0,3,0,0,0,0,0,0,exps,coefs,sc);
		  gtos_loaded=10;
		}
              if (purecart[2] == 1 && AOorder == 2)
                {
                  basis[gtos_count].IniGauss(3,nc,0,0,0,0,0,3,4,0,0,exps,coefs,sc);
                  basis[gtos_count+1].IniGauss(3,nc,0,0,0,0,0,3,5,0,0,exps,coefs,sc);
                  basis[gtos_count+2].IniGauss(3,nc,0,0,0,0,0,3,3,0,0,exps,coefs,sc);
                  basis[gtos_count+3].IniGauss(3,nc,0,0,0,0,0,3,6,0,0,exps,coefs,sc);
                  basis[gtos_count+4].IniGauss(3,nc,1,1,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+5].IniGauss(3,nc,0,0,0,0,0,3,7,0,0,exps,coefs,sc);
                  basis[gtos_count+6].IniGauss(3,nc,0,0,0,0,0,3,1,0,0,exps,coefs,sc);
                  gtos_loaded=7;
                }
              if (purecart[2] == 2 && AOorder == 2)
                {
                  basis[gtos_count].IniGauss(3,nc,3,0,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+1].IniGauss(3,nc,0,3,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+2].IniGauss(3,nc,0,0,3,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+3].IniGauss(3,nc,1,2,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+4].IniGauss(3,nc,2,1,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+5].IniGauss(3,nc,2,0,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+6].IniGauss(3,nc,1,0,2,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+7].IniGauss(3,nc,0,1,2,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+8].IniGauss(3,nc,0,2,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+9].IniGauss(3,nc,1,1,1,0,0,0,0,0,0,exps,coefs,sc);
                  gtos_loaded=10;
                }
	    }

          else if (Str=="G")
            {
              for (int j=0; j<nc; j++)
                {
                  exps[j] = basis_iStr.getNextDouble();
                  coefs[j]= basis_iStr.getNextDouble();
                }

              if (purecart[1] == 1 && AOorder == 1)
                {
                  basis[gtos_count].IniGauss(4,nc,0,0,0,0,0,0,0,4,1,exps,coefs,sc);
                  basis[gtos_count+1].IniGauss(4,nc,0,0,0,0,0,0,0,4,2,exps,coefs,sc);
                  basis[gtos_count+2].IniGauss(4,nc,0,0,0,0,0,0,0,4,3,exps,coefs,sc);
                  basis[gtos_count+3].IniGauss(4,nc,0,0,0,0,0,0,0,4,4,exps,coefs,sc);
                  basis[gtos_count+4].IniGauss(4,nc,0,0,0,0,0,0,0,4,5,exps,coefs,sc);
                  basis[gtos_count+5].IniGauss(4,nc,0,0,0,0,0,0,0,4,6,exps,coefs,sc);
                  basis[gtos_count+6].IniGauss(4,nc,0,0,0,0,0,0,0,4,7,exps,coefs,sc);
                  basis[gtos_count+7].IniGauss(4,nc,0,0,0,0,0,0,0,4,8,exps,coefs,sc);
                  basis[gtos_count+8].IniGauss(4,nc,0,0,0,0,0,0,0,4,9,exps,coefs,sc);
                  gtos_loaded=9;
                }
              if (purecart[1] == 2 && AOorder == 1)
                {
                  basis[gtos_count].IniGauss(4,nc,4,0,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+1].IniGauss(4,nc,3,1,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+2].IniGauss(4,nc,2,2,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+3].IniGauss(4,nc,1,0,3,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+4].IniGauss(4,nc,0,4,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+5].IniGauss(4,nc,3,0,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+6].IniGauss(4,nc,2,1,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+7].IniGauss(4,nc,1,2,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+8].IniGauss(4,nc,0,3,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+9].IniGauss(4,nc,2,0,2,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+10].IniGauss(4,nc,1,1,2,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+11].IniGauss(4,nc,0,2,2,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+12].IniGauss(4,nc,1,0,3,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+13].IniGauss(4,nc,0,1,3,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+14].IniGauss(4,nc,0,0,4,0,0,0,0,0,0,exps,coefs,sc);
                  gtos_loaded=15;
                }
              if (purecart[1] == 1 && AOorder == 2)
                {
                  basis[gtos_count].IniGauss(4,nc,0,0,0,0,0,0,0,4,5,exps,coefs,sc);
                  basis[gtos_count+1].IniGauss(4,nc,0,0,0,0,0,0,0,4,6,exps,coefs,sc);
                  basis[gtos_count+2].IniGauss(4,nc,0,0,0,0,0,0,0,4,4,exps,coefs,sc);
                  basis[gtos_count+3].IniGauss(4,nc,0,0,0,0,0,0,0,4,7,exps,coefs,sc);
                  basis[gtos_count+4].IniGauss(4,nc,0,0,0,0,0,0,0,4,3,exps,coefs,sc);
                  basis[gtos_count+5].IniGauss(4,nc,0,0,0,0,0,0,0,4,8,exps,coefs,sc);
                  basis[gtos_count+6].IniGauss(4,nc,0,0,0,0,0,0,0,4,2,exps,coefs,sc);
                  basis[gtos_count+7].IniGauss(4,nc,0,0,0,0,0,0,0,4,9,exps,coefs,sc);
                  basis[gtos_count+8].IniGauss(4,nc,0,0,0,0,0,0,0,4,1,exps,coefs,sc);
                  gtos_loaded=9;
                }
              if (purecart[1] == 2 && AOorder == 2)
                {
                  basis[gtos_count].IniGauss(4,nc,4,0,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+1].IniGauss(4,nc,0,4,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+2].IniGauss(4,nc,0,0,4,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+3].IniGauss(4,nc,3,1,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+4].IniGauss(4,nc,3,0,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+5].IniGauss(4,nc,1,0,3,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+6].IniGauss(4,nc,0,3,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+7].IniGauss(4,nc,1,0,3,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+8].IniGauss(4,nc,0,1,3,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+9].IniGauss(4,nc,2,2,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+10].IniGauss(4,nc,2,0,2,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+11].IniGauss(4,nc,0,2,2,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+12].IniGauss(4,nc,2,1,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+13].IniGauss(4,nc,1,2,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+14].IniGauss(4,nc,1,1,2,0,0,0,0,0,0,exps,coefs,sc);
                  gtos_loaded=15;
                }

            }
	  else if (Str=="SP")
	    {
	      for (int j=0; j<nc; j++)
		{
		  exps[j]  = basis_iStr.getNextDouble();
		  coefs[j] = basis_iStr.getNextDouble();
		  pcoefs[j]= basis_iStr.getNextDouble();
		}

	      basis[gtos_count].IniGauss(0,nc,0,0,0,0,0,0,0,0,0,exps,coefs,sc);
	      basis[gtos_count+1].IniGauss(1,nc,1,0,0,0,0,0,0,0,0,exps,pcoefs,sc);
	      basis[gtos_count+2].IniGauss(1,nc,0,1,0,0,0,0,0,0,0,exps,pcoefs,sc);
	      basis[gtos_count+3].IniGauss(1,nc,0,0,1,0,0,0,0,0,0,exps,pcoefs,sc);
	      gtos_loaded=4;
	    }
	  else
	    {
	      std::cout << "Unknown l of atomic orbital: \""<< Str <<"\"\n";
	      xmlF.exitOnFormatError(true);
	    }
	  gtos_on_current_atom+=gtos_loaded;
	  gtos_count+=gtos_loaded;

	  delete [] exps; exps=NULL;
	  delete [] coefs; coefs=NULL;
	  delete [] pcoefs; pcoefs=NULL;

	  //! check if the end of the basis for the atom (marked by "****")
	  basis_iStr.getNextWord(Str);
	  if ( Str=="****" )
	    ifEndOfBasisOnAtom=true;
	  //else: Str is the next orbital type
	}
      gtos_atom[atom_number] = gtos_on_current_atom; 
      xmlF.stepBack();
    }

  //-----------------------------------------------------

  PrintAOBasis();
}


// ======================================================================



