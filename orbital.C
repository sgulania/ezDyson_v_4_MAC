#include "orbital.h"


void Orbital::Alloc(int nb)
{
  check(!(nb>0),"Orbital::Alloc() : Invalid nbasis\n");
  
  Free();
  
  nbasis=nb;
  Lcoeff=new double[nbasis];
  Rcoeff=new double[nbasis];
  zeromap=new short[nbasis];
}

void Orbital::Free()
{
  if(IfAlloc())
   {
     if (NULL!=Lcoeff)
      {
        delete [] Lcoeff;
        Lcoeff=NULL;
      }
     if (NULL!=Rcoeff)
      {  
        delete [] Rcoeff;
        Rcoeff=NULL;
      }
     if (NULL!=zeromap)
       {  
	 delete [] zeromap;
	 zeromap=NULL;
      }
     nbasis=0;
   }
}

void Orbital::Set()
{
  nbasis=0;
  Lcoeff=NULL;
  Rcoeff=NULL;
  zeromap=NULL;
}


void Orbital::AllocAndCopyData(int nb, bool ifr, int tran, double ln, double rn, double *Lc, double *Rc, short *zm)
{
  Alloc(nb);
  ifrestr = ifr;
  transno = tran;
  Lnorm = ln;
  Rnorm = rn;
  memcpy(Lcoeff,Lc,sizeof(double)*nbasis);
  memcpy(Rcoeff,Rc,sizeof(double)*nbasis);
  memcpy(zeromap,zm,sizeof(short)*nbasis);
}

void Orbital::SetZeroMap(double thresh)
{
  for(int i=0; i<nbasis; i++)
    zeromap[i]=(fabs(Lcoeff[i])<thresh && fabs(Rcoeff[i])<thresh) ? 0 : 1;
}

Orbital::Orbital() 
{
  Set();
}

Orbital::Orbital(const Orbital& other)
{
  Set();
  AllocAndCopyData(other.nbasis,other.ifrestr,other.transno,other.Lnorm,other.Rnorm,other.Lcoeff,other.Rcoeff,other.zeromap);
}

Orbital::Orbital(int nb, bool ifr, int tran, double ln, double rn, double *Lc, double *Rc)
{
  Set();
  Alloc(nb);
  ifrestr = ifr;
  transno = tran;
  Lnorm = ln;
  Rnorm = rn;
  memcpy(Lcoeff,Lc,sizeof(double)*nbasis);
  memcpy(Rcoeff,Rc,sizeof(double)*nbasis);
  //AllocAndCopyData(nb,ifr,tran,Lc,Rc);
  SetZeroMap();
}

Orbital::~Orbital() 
{ 
  Free();
}

Orbital& Orbital::operator=(const Orbital& other)
{
  if(this!=&other)
    AllocAndCopyData(other.nbasis,other.ifrestr,other.transno,other.Lnorm,other.Rnorm,other.Lcoeff,other.Rcoeff,other.zeromap);
  return *this; 
}


//! Print Orbital coeffs in AO basis.
void Orbital::Print(FILE *fout) const
{
  fprintf(fout,"Left norm (sqrt|dyson_l x dyson_l|)=%lf\n",Lnorm);
  fprintf(fout,"Right norm (sqrt|dyson_r x dyson_r|)=%lf\n",Rnorm);
  fprintf(fout,"LeftxRight=%lf\n",Lnorm*Rnorm);

  fprintf(fout,"\nDysonMOs in AO basis\n");
  fprintf(fout,"  AO     [LEFT]              [RIGHT]        ZERO_MAP\n");
  for (int i=0; i<nbasis; i++)
    fprintf(fout,"%3d %15.8lf     %15.8lf      %d\n",i+1,Lcoeff[i],Rcoeff[i],zeromap[i]);
  fprintf(fout,"-----------------------------------------------------------------------------------\n");
  fflush(fout);
}

//! Add up all basis fns with proper coefs to get orbital value at x,y,z.
void Orbital::CalcDyson_xyz(double &Ldys_value, double &Rdys_value, double x, double y, double z) const
{
  static AOBasis *pbasis=&(TheAOBasis());
  static int nbasis = pbasis->NBasis();

  double onegto;

  //Initialize to zero.
  Ldys_value = 0.0;
  Rdys_value = 0.0;

  for(int n=0,k=0,i=0; i<nbasis; i++,k++)
    {
      if(k==pbasis->GetGTOsAtom(n))
       {
         n++;
         k=0;
       }
     
      if(zeromap[i]) //Do calcs only if non-zero coefficient
	{
	  onegto = pbasis->CalcGaussValue(i,n,x,y,z);
	  Ldys_value += Lcoeff[i]*onegto;
	  Rdys_value += Rcoeff[i]*onegto;
	}
    }

  //Slater test
/*    double rsq = x*x+y*y+z*z;
    Ldys_value = 0.23*z*(2.57*exp(-20.96*rsq)+2.41*exp(-4.8*rsq)+1.865*exp(-1.46*rsq));
    Ldys_value +=0.38*z*0.575*exp(-0.483*rsq)+0.43*z*0.128*exp(-0.146*rsq);
    Ldys_value +=0.22*z*0.02856*exp(-0.0483*rsq)+0.05*z*0.00637*exp(-0.01319*rsq);

    Ldys_value = z*exp(-1.2*BOHR_TO_ANGS*sqrt(rsq));
    Rdys_value = Ldys_value;  */
}

double Orbital::GetOrbitalNormAndCenter(double *xyzcenter, const XYZGrid& grid) const
{
  double ldys,rdys,ldyson_sq,rdyson_sq;
  double lnorm=0.0,rnorm=0.0,lrnorm=0.0, dV=grid.DXYZ(X)*grid.DXYZ(Y)*grid.DXYZ(Z);
  double xi, yi, zi;
  int i;

  memset(xyzcenter,0,sizeof(double)*XYZ);

  for(i=0; i<grid.NXYZ(); i++)
    {
      grid.GetPoint(i,xi,yi,zi,X,Y,Z);
      CalcDyson_xyz(ldys,rdys,xi,yi,zi);
      ldyson_sq=ldys*ldys;
      rdyson_sq=rdys*rdys;
      lnorm+=ldyson_sq;
      rnorm+=rdyson_sq;
      lrnorm+=ldys*rdys;
      
      xyzcenter[X] += xi*ldyson_sq;
      xyzcenter[Y] += yi*ldyson_sq;
      xyzcenter[Z] += zi*ldyson_sq;
    }
  
  for (i=0; i<XYZ; i++)
    xyzcenter[i]/=lnorm;

  lnorm*=dV;
  rnorm*=dV;
  lrnorm*=dV;
  lnorm=sqrt(lnorm);
  rnorm=sqrt(rnorm);
  lrnorm=sqrt(lrnorm);

  fprintf(stdout,"\nNorm of the left Dyson orbital integrated on the grid: %15.8lf\n", lnorm);
  fprintf(stdout,"Norm of the right Dyson orbital integrated on the grid: %15.8lf\n", rnorm);
  fprintf(stdout,"Left-Right Norm of the Dyson orbital integrated on the grid: %15.8lf\n", lrnorm);
  fprintf(stdout,"Xavg Yavg Zavg:  %lf  %lf  %lf a.u.\n",xyzcenter[X],xyzcenter[Y],xyzcenter[Z]);
  fprintf(stdout,"-----------------------------------------------------------------------------------\n");

  //Now compute X^2 Y^2 and Z^2
  double second_mom[XYZ];
  second_mom[X]=second_mom[Y]=second_mom[Z]=0.0;

  for(i=0; i<grid.NXYZ(); i++)
    {
      grid.GetPoint(i,xi,yi,zi,X,Y,Z);
      CalcDyson_xyz(ldys,rdys,xi,yi,zi);
      ldyson_sq=ldys*ldys;
      lnorm+=ldyson_sq;

      xi-=xyzcenter[X];
      yi-=xyzcenter[Y];
      zi-=xyzcenter[Z];

      second_mom[X] += xi*xi*ldyson_sq;
      second_mom[Y] += yi*yi*ldyson_sq;
      second_mom[Z] += zi*zi*ldyson_sq;
    }

  double size=0.0;

  for (i=0; i<XYZ; i++)
    {
      second_mom[i]/=lnorm;
      second_mom[i]*=(BOHR_TO_ANGS*BOHR_TO_ANGS);
      size+=second_mom[i];
    }

  fprintf(stdout,"<X^2> <Y^2> <Z^2> <R^2>:  %lf  %lf  %lf %lf Angs\n",
	  second_mom[X],second_mom[Y],second_mom[Z],size);
  fprintf(stdout,"-----------------------------------------------------------------------------------\n");
  fflush(stdout);
  return lnorm;
}



