#include "klmgrid.h"

void KLMPoints::Alloc(int nk, int nlm)
{
  check(!(nk*nlm>0),"SphHarm::Alloc() : Invalid nsphwaves\n");

  Free();
  double gsize=0.0;

  nkv = nk;
  nlmv = nlm;
  kv = new double[nkv];
  gsize+=sizeof(double)*nkv;
  check(!(kv>NULL), "KLMPoints::Alloc() : Allocation failed");
  lv = new int[nlmv];
  gsize+=sizeof(int)*nlmv;
  check(!(lv>NULL), "KLMPoints::Alloc() : Allocation failed");
  mv = new int[nlmv];
  gsize+=sizeof(int)*nlmv;
  check(!(mv>NULL), "KLMPoints::Alloc() : Allocation failed");
  printf("\nKLMPoints::Alloc:  %lf MB allocated\n",gsize/MEGABITE);
}

void KLMPoints::Free()
{
  if (IfAlloc())
   {
     if (NULL!=kv)
      {
        delete [] kv;
        kv=NULL;
      }
     if (NULL!=lv)
      {
        delete [] lv;
        lv=NULL;
      }
     if (NULL!=mv)
      {
        delete [] mv;
        mv=NULL;
      }
     nkv = 0;
     nlmv = 0;
   }
}

void KLMPoints::Set()
{
  nkv=0;
  nlmv=0;
  kv=NULL;
  lv=NULL;
  mv=NULL;
}


KLMPoints::KLMPoints()
{
  Set();
}

KLMPoints::KLMPoints(const char* xmlFileName)
{
  Set();
  int nkv=0, nlm=0;

  simpleXMLparser xmlF;
  xmlF.assignFile(xmlFileName); 

  xmlF.reset().node("root").node("free_electron");

  lmax = xmlF.getIntValue("l_max");
  xmlF.node("k_grid");
  nkv = xmlF.getIntValue("n_points");
  kmin = xmlF.getDoubleValue("min") * EV_TO_HAR;
  kmax = xmlF.getDoubleValue("max") * EV_TO_HAR;
  extern SPH& theSPH() ;
  check(!(lmax<theSPH().LMax()),"Invalid lmax");

  nlm = (lmax+1)*(lmax+1);
  Alloc(nkv,nlm);
  CalcKGrid();
  CalcLMPoints();
}


KLMPoints::~KLMPoints() 
{ 
  Free();
}

//!Generate k grid: kv[nkv] from energy points
void KLMPoints::CalcKGrid()
{
  if (nkv==1)
  {
    kv[0] = kmin;
    kv[0] = sqrt(2.0*kv[0]);
  }
  else
   {
    double dk = (kmax-kmin)/(nkv-1);
    for (int i=0; i<nkv; i++)
      {
        kv[i] = kmin + i*dk;
        kv[i] = sqrt(2.0*kv[i]);
      }
   }
}

//!From l and m matrices. lv matrix does not have unique values: 
//! same l for each m value,e.g. lv[0 1 1 1], mv[0 -1 0 1]
void KLMPoints::CalcLMPoints()
{
  for (int n=0,l=0; l<=lmax; l++)
    for (int m=-l; m<=l; m++)
      {
	lv[n]=l;
	mv[n]=m;
	n++;
      }
}

//! Print k,l,m grid info.
void KLMPoints::PrintGridInfo(FILE *outfile) const
{
  fprintf(outfile,"\nKLM Grid\n");
  fprintf(outfile,"k-grid   %d pts  %lf  %lf eV\n",nkv,kmin/EV_TO_HAR,kmax/EV_TO_HAR);
  fprintf(outfile,"max l:  %d\n",lmax);
  fprintf(outfile,"no. l,m values:  %d\n",nlmv);
  fflush(outfile);
}

//! Get k and l value for the n-th spherical wave, by absolute address.
void KLMPoints::GetKLPoint(int absaddr, double &kpt, double &lpt)
{
  int ik, jl;
  ik = (absaddr-1)/nlmv;
  jl = absaddr-1-ik*nlmv;
  kpt = kv[ik];
  lpt = lv[jl];
}


