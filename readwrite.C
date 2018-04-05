#include "readwrite.h"

//! Which MOs in AO basis to read from QChem molden_format output?
void MONumbers(int *monumbers, const char* xmlFileName, int nomos)
{
  simpleXMLparser xmlF;
  xmlF.assignFile(xmlFileName); 

  xmlF.reset().node("root").node("job_parameters");

  std::istringstream tmp_iStr;

  if (nomos>0)
    {
      tmp_iStr.str( xmlF.value("MOs_to_plot") );
      for (int i=0; i<nomos; i++)
	{
	  tmp_iStr >> monumbers[i];
	  xmlF.exitOnFormatError( tmp_iStr.fail() );
	}

      std::cout << nomos << " orbitals to plot with the following numbers: ";
      for (int i=0; i<nomos; i++)
	std::cout << monumbers[i] << ' ';
      std::cout << '\n';
    } 
  else
    std::cout << "No orbitals were requested to plot.\n";
}


//!
void ReadOrbitalFromFile(Orbital& orb, int nb, const char* xmlFileName)
{
  simpleXMLparser xmlF;
  xmlF.assignFile(xmlFileName); 

  xmlF.reset().node("root").node("job_parameters");
//  bool ifrestr = !( xmlF.getBoolValue("unrestricted") );

/*! SG: Here I will override the above and set ifrestr = true always
 * This means regardless of what is in the input we will always read only two
 * dyson orbitals (left-right) and not four (left/alpha-left/beta-right/alpha-right/beta).
 * This is because in ground to IP cases we only ionize one electron (alpha)
 * and so beta Dyson orbitals are always 0.
*/
  bool ifrestr = true;

  xmlF.reset().node("root").node("job_parameters");
  int transno = xmlF.getIntValue("Dyson_MO_transitions");

  double ln, rn;

  //! Read dyson Norm
  xmlF.reset().node("root").node("dyson_molecular_orbitals");
  
  //transno = 1,2...;
  xmlF.node("DMO",transno*2-1);
  ln=xmlF.getDoubleValue("norm");
  xmlF.stepBack();

  if (ifrestr==true)
    xmlF.node("DMO",transno*2);
  else
    xmlF.node("DMO",transno*2+1);
  rn=xmlF.getDoubleValue("norm");
  xmlF.stepBack();

  double *tmpR=new double[nb];
  double *tmpL=new double[nb];

  //! Read coeffs of Dyson orbs in AO basis for n-th transition only.

  //position =0,1,...
  int position = 2*(transno-1);   //left-right
  if (ifrestr==false)
    position = 2*position;        //alpha and beta

  std::stringstream tmp_iStr;
  xmlF.reset().node("root").node("dyson_molecular_orbitals");

  xmlF.node("DMO",position+1);
  tmp_iStr.str( xmlF.value("text") );
  for (int j=0; j<nb; j++)
    tmp_iStr >> tmpL[j];
  xmlF.stepBack();

  position++;
  if (ifrestr==false) position++;

  xmlF.node("DMO",position+1);
  tmp_iStr.str( xmlF.value("text") );
  for (int j=0; j<nb; j++)
    tmp_iStr >> tmpR[j];

  orb=Orbital(nb,ifrestr,transno,ln,rn,tmpL,tmpR);
  delete [] tmpR;
  delete [] tmpL;

  orb.Print();
}


//! Read MO to plot from qchem.
void ReadMOFromFile(Orbital &mo, int nb, int monumber, const char* xmlFileName)
{
  simpleXMLparser xmlF;
  xmlF.assignFile(xmlFileName); 

  xmlF.reset().node("root").node("job_parameters");
  bool ifrestr = !( xmlF.getBoolValue("unrestricted") );
  
  double *tmpA= new double[nb];
  double *tmpB= new double[nb];
 
  //! Read MOs in AO basis coeffs from qchem.
  xmlF.reset().node("root").node("molecular_orbitals").node("alpha_MOs").node("MO",monumber);
  std::stringstream tmp_iStr;
  tmp_iStr.str( xmlF.value("text") );
  for (int k=0; k<nb; k++)
    tmp_iStr >> tmpA[k];
  if (ifrestr==false)
    {
      xmlF.reset().node("root").node("molecular_orbitals").node("beta_MOs").node("MO",monumber);
      tmp_iStr.str( xmlF.value("text") );
      for (int k=0; k<nb; k++)
	tmp_iStr >> tmpB[k];
    }

  mo=Orbital(nb,ifrestr,monumber,1.0,1.0,tmpA,tmpB);
  delete [] tmpA;
  delete [] tmpB;

  mo.Print();
}


//!Read info about averaging over molec orientations in lab frame.
//SG: Now we only give a choice to the user to use avg or not
MolAvg ReadMolecOrientAvg(const char* xmlFileName)
{
  simpleXMLparser xmlF;
  xmlF.assignFile(xmlFileName); 

  std::string tmpStr;
  MolAvg molavg;

  xmlF.reset().node("root");
  if ( xmlF.CheckSubNode("averaging") )
    {
      tmpStr=xmlF.node("averaging").value("method");
      if (tmpStr=="noavg")
	molavg=NOAVG;
      else if (tmpStr=="avg")
	molavg=AVG;
      else if (tmpStr=="num")
        molavg=NUM;
      else
	xmlF.exitOnFormatError(true);
    }
  else
    molavg=AVG;

  fprintf(stdout,"\nAveraging %d\n",molavg); fflush(stdout);

  return molavg;
}


//!x,y,z lab frame components of ionization laser polarization.
void ReadRIonz(double &ie, double *rioniz, const char* xmlFileName)
{
  simpleXMLparser xmlF;
  xmlF.assignFile(xmlFileName); 

  xmlF.reset().node("root").node("laser");

  // read and convert to a.u.
  ie = xmlF.getDoubleValue("ionization_energy") / HAR_TO_EV;

  // SG: read laser polarization (now only applicable for fixed molecular orientation)
  xmlF.node("laser_polarization");
  rioniz[0]=xmlF.getDoubleValue("x");
  rioniz[1]=xmlF.getDoubleValue("y");
  rioniz[2]=xmlF.getDoubleValue("z");

  std::cout << "IE=" << ie*HAR_TO_EV <<"eV;  Polarization=(" << rioniz[0] << ',' << rioniz[1] << ',' << rioniz[2] << ")\n";
}


//! Print
void PrintOrbitalsScanZ(Orbital *mostoplot, int nomos, XYZGrid& labgrid)
{
  FILE *file_orbs=fopen("mosplot.dat","w");
  double *gridptr_z=labgrid.GetGridPtr(Z);
  double xyz[XYZ], Ldys_value, Rdys_value;

  xyz[X]=xyz[Y]=0.0;

  fprintf(file_orbs, "Z,A RDys_val RD*Z \n");

  for(int i=0; i<labgrid.NPoints(Z); i++)
    {
      xyz[Z]=gridptr_z[i];
      fprintf(file_orbs,"%lE",xyz[Z]*BOHR_TO_ANGS);
      for (int j=0; j<nomos; j++)
        {
         mostoplot[j].CalcDyson_xyz(Ldys_value,Rdys_value,0.0,0.0,xyz[Z]);
         fprintf(file_orbs," %lE %lE   ", Rdys_value, Rdys_value*xyz[Z]);
        }
      fprintf(file_orbs,"\n");
    }

  fclose(file_orbs);
}











//-------------------------------------------------------------------------------------------------------------

/* VM:

//!x,y,z lab frame components of excitation laser polarization.
void ReadRExc(double *rex, FILE *input)
{
  rewind(input);
  char str[256];

  while (strcmp(str,"exc_laser")!=0 && !feof(input))
      {
        fscanf(input,"%s\n",str);
        LowCase(str);
      }
  fscanf(input,"%lf",&rex[0]);
  fscanf(input,"%lf",&rex[1]);
  fscanf(input,"%lf",&rex[2]);
  fprintf(stdout,"Rexc %lf %lf %lf\n",rex[0],rex[1],rex[2]); fflush(stdout);
}


//!Molecular frame components of excitation dipole moment(x,y,z-Left,x,y,z-Right).
void ReadMExc(double *mex, FILE *input)
{
//Read excitation dipole moment.
  rewind(input);
  char str[256];

  while (strcmp(str,"exc_dipole_lr")!=0 && !feof(input))
      {
        fscanf(input,"%s\n",str);
        LowCase(str);
      }
  fscanf(input,"%lf",&mex[0]);
  fscanf(input,"%lf",&mex[1]);
  fscanf(input,"%lf",&mex[2]);
  fscanf(input,"%lf",&mex[3]);
  fscanf(input,"%lf",&mex[4]);
  fscanf(input,"%lf",&mex[5]);
  fprintf(stderr,"Mexc %lf %lf %lf %lf %lf %lf\n",mex[0],mex[1],mex[2],mex[3],mex[4],mex[5]);
}

*/
