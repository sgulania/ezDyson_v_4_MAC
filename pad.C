#include "pad.h"
#include "ylm.h"

PAD::PAD(int nth, int nk) : ntheta(nth), nkv(nk) 
{
  pad=new double[nkv*ntheta];
  crosspad=new double[nkv*ntheta];
  totalpad=new double[nkv*ntheta];  
  memset(pad,0,nkv*ntheta*sizeof(double));
  memset(crosspad,0,nkv*ntheta*sizeof(double));
  memset(totalpad,0,nkv*ntheta*sizeof(double));

}

PAD::~PAD()
{
  delete [] pad;
  delete [] crosspad;
  delete [] totalpad;
}


//! Caution: this is unscaled PAD
// SG: PAD calculated for averaged Cklm
void PAD::CalcPad(CklmCoeff& allCklm, MolAvg molavg) const
{
  printf("\nCalculating PADs and cross sections\nntheta=%d nkv=%d\n", ntheta,nkv); fflush(stdout);

  double delta=DeltaTheta(), theta, tmp, tmp2, coeff;
  int lmax=allCklm.GetLMax();

  for (int k=0; k<nkv; k++) //LOOP BY k
    for (int t=0; t<ntheta; t++) //Loop by theta
      {
	theta = t*delta;
	//Loops by i=(l,m) 
	for (int i=0, l=0; l<=lmax; l++)
	  for (int m=-l; m<=l; m++, i++)
	    {
	      tmp = ThetaYlm(theta,l,m);
	      for (int l2=0; l2<=lmax; l2++)
		{
		  if (l2==l) //diagonal contributons
		    {
                      coeff = allCklm.GetCklmSq(k,i);
		      pad[k*ntheta+t] += coeff*tmp*tmp;//*costsq;
		    }
		  else //cross-terms
		    {
		      tmp2 = ThetaYlm(theta,l2,m);
                      coeff = allCklm.GetCrossCklm(k,i,l2);
		      crosspad[k*ntheta+t] += coeff*tmp*tmp2;//*costsq;
		    }
		}
	    }
	totalpad[k*ntheta+t] = (pad[k*ntheta+t] + crosspad[k*ntheta+t]);
      }
  printf("PAD: done\n");
}

void PAD::Print(const KLMPoints& klmgrid) const
{
  //AIK: need to rewrite appropriately
  FILE *outfile = fopen("pad.dat","w");
  double delta=DeltaTheta(), theta;
  
  for (int k=0,no=0; k<nkv; k++) //LOOP BY k
    {
      double kwave=klmgrid.GetKV(k);
      double ene=kwave*kwave/2.*HAR_TO_EV;
      fprintf(outfile,"Energy=%lf eV\n",ene); 
      fprintf(outfile,"theta,rad   sq_contrib(lm)   cross_contrib(lml’m)   total_PAD\n");
      for (int t=0; t<ntheta; t++,no++) //Loop by theta
	{
	  theta = t*delta;
	  fprintf(outfile," %13.6lf %17.6lE %17.6lE %17.6lE\n",theta,pad[no],crosspad[no],totalpad[no]);
	}
      fprintf(outfile,"__________________________________\n");
    }
  fclose(outfile);
}

//! Calculates total xsec by numeric integration, need it for consistency check
double PAD::TotalXSec(int kv) const
{
  double xsec=0.0, delta=DeltaTheta();

  int t=0;
  xsec+=totalpad[kv*ntheta+t]*sin(t*delta)*delta/2.;   
  
  for(t=1; t<ntheta-1; t++)
    xsec+=totalpad[kv*ntheta+t]*sin(t*delta)*delta; 
  
  xsec+=totalpad[kv*ntheta+t]*sin(t*delta)*delta/2.; 
  
  return xsec*2.*M_PI; //f integration 2PI;
}

//! returns parallel x-sec, i.e. theta=0
double PAD::XSec_par(int kv) const
{
  return totalpad[kv*ntheta+0]; 
}

//! returns perpendicular x-sec, i.e. theta=pi/2
double PAD::XSec_perp(int kv) const
{
  double delta=DeltaTheta();
  int midpoint=ntheta/2;
  double theta=midpoint*delta;

  //printf("PAD::XSec_Perp Ntheta=%d midpoint=%d value/(pi/2)=%lf",ntheta, midpoint,theta/(M_PI/2.));
  
  //AIK: not sure about 2PI
  if(fabs(theta/(M_PI/2.)-1.0)<0.00001)
    return totalpad[kv*ntheta+midpoint];
  else
    {
      printf("PAD::XSec_Perp() : Warning: The midpoint is not pi/2\n");
      return 0.5*(totalpad[kv*ntheta+midpoint]+totalpad[kv*ntheta+midpoint+1]);
    }
}
