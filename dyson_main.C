#include "dyson_main.h"

bool dyson_main(const char* xmlFileName)
{
  std::cout << "Starting ezDyson version 4.0\n\n";

  std::cout << "Reading the input file\n\n";
  simpleXMLparser xmlF;
  xmlF.assignFile(xmlFileName); 

  time_t start,end;
  time (&start);
  
  //! Read atoms and AO basis information from QChem output.
  TheAOBasis().IniAOBasis(xmlFileName); 
  printf("AO BASIS done\n\n"); fflush(stdout);

  //! Initialize/read info about Dyson MOs.
  Orbital dysorb;
  ReadOrbitalFromFile(dysorb, TheAOBasis().NBasis(), xmlFileName);
  printf("DYSON ORBITAL done\n\n"); fflush(stdout);

  //! Read MOs to plot.
  xmlF.reset().node("root").node("job_parameters");
  int nomos = xmlF.getIntValue("number_of_MOs_to_plot");
  int *monumbers = new int[nomos];
  MONumbers(monumbers,xmlFileName,nomos);
  Orbital *mostoplot = new Orbital[nomos];
  for (int i=0; i<nomos; i++)
    ReadMOFromFile(mostoplot[i],TheAOBasis().NBasis(),monumbers[i],xmlFileName);
  delete [] monumbers;

  //!Read lab frame x,y,z grid info.
  XYZGrid labgrid(xmlFileName);
  labgrid.PrintGridInfo(stdout);
  //labgrid.Print("labgrid.chk");

  //! Read averaging method from input.
  MolAvg molavg = ReadMolecOrientAvg(xmlFileName);

  //! Read direction of ionization laser and ionization energy (read in eV, convert to au)
  double RIonz[XYZ], IE;
  ReadRIonz(IE, RIonz, xmlFileName);

  //! Calc center of (un)averaged DysonMO(left one).
  double dyscenter[XYZ], dnorm;
  dnorm=dysorb.GetOrbitalNormAndCenter(dyscenter,labgrid);
  TheAOBasis().ShiftGeom(dyscenter);

  //! Plot MOs and Z*MOs after labframe center shift.  
  if (nomos)
    {
      PrintOrbitalsScanZ(mostoplot,nomos,labgrid);
      delete [] mostoplot;
    }

  //! Read info for generating the k,l,m points of the spherical waves.
  KLMPoints klmgrid(xmlFileName);
  klmgrid.PrintGridInfo();

  //! Pre-compute all Rkl at labgrid points.
  //SG: Old implementation: required high mem... uncomment next line for high mem version
  //RadialFunction allRkl(xmlFileName,labgrid,klmgrid);
  //allRkl.Print("rkl.chk");
  //allRkl.NormMatrix(100.0,1001,klmgrid);
  int nkv = klmgrid.NKV();
  if (molavg == NUM)
    {
       NumEikr allEikr;
       std::cout << "Data is loaded from the input file " << xmlFileName << " \n\n";

       AngleGrid anggrid;

       allEikr.IniNumEikr(labgrid,dysorb,anggrid,klmgrid);
       xmlF.reset().node("root").node("job_parameters");
       double spinfactor, orbfactor;
       spinfactor=xmlF.getDoubleValue("spin_degeneracy");
       orbfactor=xmlF.getDoubleValue("orbital_degeneracy");

       double norm=(1.)/(2*M_PI*C_AU);
       double dyson_norm=dysorb.GetLRNorm();
       //Calculate the conversion factors for the xsections
       double au_to_cm_sq = pow(BOHR_TO_ANGS*pow(10.,-8),2);
       printf("\nConversion factors for x-sections:\n");
       printf("Xsec(cm^-2) = Xsec(au^-2) * %10lE\n", au_to_cm_sq);
       printf("Xsec(au^-2) = Xsec(cm^-2) * %10lE\n", 1./au_to_cm_sq);
     
     
       printf("\n\nPADs:\n");
       printf("NOTE: PADs are only correct in molecular frame for Z-polarized light\n");
       printf("____________________________________________________________\n");
       printf("E=IE+E_k,eV   Sigma_par     Sigma_perp    Sigma_tot     beta\n");
       for (int k=0; k<nkv; k++)
         {
           double kwave=klmgrid.GetKV(k);
           double ene=kwave*kwave/2. + IE;  //energy of laser, in AU
     
         //Bethe salpter.
           double scale=norm*ene*kwave*dyson_norm*spinfactor*orbfactor;
     
           printf("%lf", ene*HAR_TO_EV);
           double sigma_par=allEikr.GetCPar(k)*scale;
           double sigma_perp=allEikr.GetCPerp(k)*scale;
           double sigma_tot=(sigma_par+2.*sigma_perp)*4.0*M_PI/3.;
           double beta=2.*(sigma_par-sigma_perp)/(sigma_par+2.*sigma_perp);
     
           printf(" %13.6lf %13.6lf %13.6lf %13.6lf\n",sigma_par,sigma_perp,sigma_tot,beta);
         }
     }
   else if (molavg==AVG || molavg==NOAVG)
    {
      //SG: Radial Function from input
       xmlF.reset().node("root").node("free_electron");
       double radfn;
       radfn=xmlF.getDoubleValue("charge_of_ionized_core");

       std::cout << "Charge of ionized core = " << radfn <<'\n';

       //! Calc avg |Cklm|^2.
       CklmCoeff allCklm; 
       std::cout << "Data is loaded from the input file " << xmlFileName << " \n\n";

       if (NOAVG==molavg)
         {
           allCklm.IniCklmCoeff(labgrid,dysorb,RIonz,radfn,klmgrid);
         }
        else // if (AVG==molavg)
         {
           allCklm.IniCklmCoeff(labgrid,dysorb,radfn,klmgrid);
         }

       int lmax = klmgrid.LMax();
       int ntheta = 101;
       PAD totalpad(ntheta,nkv);
       totalpad.CalcPad(allCklm, molavg);

       xmlF.reset().node("root").node("job_parameters");
       double spinfactor, orbfactor;
       spinfactor=xmlF.getDoubleValue("spin_degeneracy");
       orbfactor=xmlF.getDoubleValue("orbital_degeneracy");

       double *xsectot = new double[nkv];
       double *xsec1 = new double[nkv];
       memset(xsec1,0,sizeof(double)*nkv);
       memset(xsectot,0,sizeof(double)*nkv);
       for (int k=0; k<nkv; k++) //LOOP BY k
        {
          xsec1[k]=allCklm.GetSumCklm(k);
        }
       //! Bethe-Salpeter factor of 4*pi^2/c devided by factor from Sakurai of (2*pi). See updated derivation
       double norm=(8*M_PI*M_PI*4)/(3*C_AU);


       printf("\n\nTotal crossection vs E(eV), in au, using norm=(1/C_AU)*E*k\n");
       double dyson_norm=dysorb.GetLRNorm();

       printf("[E*2*Pi/(c*k) * |Cklm|^2]\n");
       printf("E=IE+E_k,eV   xsec,a.u.\n");
       for (int k=0; k<nkv; k++)
         {
           double kwave=klmgrid.GetKV(k);

           double ene=kwave*kwave/2. + IE;  //energy of laser, in AU

           double scale=norm*kwave*spinfactor*orbfactor;
           xsec1[k] *= scale*dyson_norm;
           xsectot[k] = (xsec1[k]);

           printf("%lf", ene*HAR_TO_EV);
             printf(" %13.6lf\n",xsectot[k]*ene);
         }
       printf("\n\nTotal cross-section / E (in a.u.) vs. Ek (in eV):\n");
       printf("To be used with xsecFCF script for incorporating FCFs\n");
       printf("E_k,eV        xsec/E,a.u. \n");

       for (int k=0; k<nkv; k++)
         {
           double kwave=klmgrid.GetKV(k);
           double ene=kwave*kwave/2.;  //energy of the electron

           printf("%lf", ene*HAR_TO_EV);
           printf("%14.6lf\n",xsectot[k]);
         }
       //Calculate the conversion factors for the xsections
       double au_to_cm_sq = pow(BOHR_TO_ANGS*pow(10.,-8),2);
       printf("\nConversion factors for x-sections:\n");
       printf("Xsec(cm^-2) = Xsec(au^-2) * %10lE\n", au_to_cm_sq);
       printf("Xsec(au^-2) = Xsec(cm^-2) * %10lE\n", 1./au_to_cm_sq);
     
     
       printf("\n\nPADs:\n");
       printf("NOTE: PADs are only correct in molecular frame for Z-polarized light\n");
       printf("____________________________________________________________\n");
       printf("E=IE+E_k,eV   Sigma_par     Sigma_perp    Sigma_tot     beta\n");
       for (int k=0; k<nkv; k++)
         {
           double kwave=klmgrid.GetKV(k);
           double ene=kwave*kwave/2. + IE;  //energy of laser, in AU
     
         //Bethe salpter.
           double scale=norm*ene*kwave*dyson_norm*spinfactor*orbfactor;
     
           printf("%lf", ene*HAR_TO_EV);
           double sigma_par=totalpad.XSec_par(k)*scale;
           double sigma_perp=totalpad.XSec_perp(k)*scale;
           double sigma_tot=(sigma_par+2.*sigma_perp)*4.0*M_PI/3.;
           double beta=2.*(sigma_par-sigma_perp)/(sigma_par+2.*sigma_perp);

           printf(" %13.6lf %13.6lf %13.6lf %13.6lf\n",sigma_par,sigma_perp,sigma_tot,beta);
         }
           //Write PADs on disk
           totalpad.Print(klmgrid);
          
            delete [] xsec1;
            delete [] xsectot;

    }

  time (&end);
  double dif = difftime (end,start);
  printf("\nJob time: %lf seconds\n",dif);
  
  return true;
}



