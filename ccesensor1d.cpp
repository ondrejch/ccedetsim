/* ROOT macro to 1D case detector simualtion with CCE sum		*
 * g++ ccesensor1d.C -o ccesensor1d.exe `root-config --cflags --libs`
 *      Ondrej, April 6 2013                                            */
#include <vector>
#include <iostream>
#include <fstream>
#include <string.h>
#include <TF1.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TFile.h>
#include <TNtuple.h>
using namespace std;

const int debug=1;

/* Default detector parameters, can be modified by parameter input file	*
Format: plain ASCII file with the following entries, NUMBERS ONLY
<W>	<Nsigma>
<Efield>	<mu_e>	<t_e>	<mu_h>	<t_h>
<Nccebins>	<ccemin>	<ccemax>					*/

double W 	= 0.01;//178; 	// sensor length [cm]
double Nsigma 	= 0.1; 		// beam attenuation coef in the detector \sigma*N [cm^-1]
// CCE parameters
double Efield = 10000.;		// 10000 V/cm
double mu_e = 5e3;		// electron mobilitymu = 5000  cm^2/V s 
double t_e  = 1e-5;		// electron collection time = 10^-6 s
double mu_h = mu_e; //5e3;		// hole mobility mu = 5000  cm^2/V s
double t_h  = t_e;  //1e-6;		// hole collection time = 10^-6 s
double ccemin = .0;		// low cutoff of the CCE histogram [MeV]
double ccemax = 4.8;		// high cutoff of the CCE histogram [MeV]
int Nccebins = (ccemax-ccemin)*1e2;		// bins for CCE sums

TF1 *fbeamatt; 		// beam attenuation function
TF1 *fcce; 		// CCE weighting function
TH1D *hcce;		// histogram for CCE sums
TNtuple *nt;		// detailed info about events 
TFile *rout;		// ROOT output file for the histogram, ntuples, etc.
TRandom3 *rndgen;	// fast random number generator, period 2^19937-1
class dpt { 		// datapoint class for Bragg curves
 public:
  double r, E;		// radius and energy deposit
  dpt(double _r, double _E):r(_r),E(_E){};
  void print() { cout << r << " " << E <<endl;};
};
vector<dpt> braggA;	// Bragg curve for alphas
vector<dpt> braggT;	// Bragg curve for tritons
int maxA, maxT; 	// sizes of vectors

void initBragg() {
  const char* filetritons = "iTRITONS.DAT";
  const double Etriton = 2.73;	// triton energy [MeV]
  const char* filealphas  = "iALPHAS.DAT"; 
  const double Ealpha = 2.05; 	// alpha energy  [MeV]
  double r,e,e2,edepsum=0;
// read in triton Bragg curve
  ifstream fin(filetritons);
  while (!fin.eof()) {
    fin >> r >> e;// >> e2;    
    //e += e2; 		// add recoils to ionizations
    r *= 1e-8;	        // Angstrom -> cm
    dpt v(r,e);
    braggT.push_back(v);
    edepsum += e;
  }
  fin.close();
  maxT = braggT.size();
  // renormalize 
  double scalefactor = Etriton/edepsum;
  for(int j=0; j<maxT; j++) braggT[j].E *= scalefactor;
// read in alpha bragg curve
  fin.open(filealphas);
  edepsum=0;
  while (!fin.eof()) {
    fin >> r >> e;// >> e2;
//    e += e2;              // add recoils to ionizations
    r *= 1e-8;            // Angstrom -> cm
    dpt v(r,e);
    braggA.push_back(v);
    edepsum += e;
  }
  fin.close();  
  maxA = braggA.size();
  scalefactor = Ealpha/edepsum;
  for(int j=0; j<maxA; j++) braggA[j].E *= scalefactor;
  cout << "Bragg curves read in. " << endl;
}

void initPars(char* parfname="mydet.par") {
  ifstream fin;
  fin.open(parfname);
  if(fin.is_open()) { 	// if parameter file exists, read it in
     fin >> W >> Nsigma;
     fin >> Efield >> mu_e >> t_e >> mu_h >> t_h;
     fin >> Nccebins >> ccemin >> ccemax;
     cout << "Detector parameters read in from file " << parfname <<endl;
  }
  fin.close();
  cout << "* Detector parameters, 1D case *" << endl;
  cout << "* W = " << W << " cm" <<endl;
  cout << "* Beam attenuation coef = "<<Nsigma<<" cm^-1" << endl;
  cout << "* Detector electric field = "<<Efield<<" V/cm"<< endl;
  cout << "* mobility for e & h = " << mu_e << " & " << mu_h << " cm^2/V s" << endl;
  cout << "* collection time for e & h = " << t_e << " & " << t_h << " s" << endl;
  cout << "* Histogram for deposit sums = " << Nccebins <<" bins" << endl;
  cout << "* Histogram min & max = " << ccemin << " & " << ccemax << " MeV"<<endl;
}
  
void ccesensor1d(long Nneutrons=10, char* outfname="mydet.root", char* parfname="mydet.par") {
  cout << "Running Sensor1D with "<<Nneutrons<<" neutrons" <<endl;
  rout = new TFile(outfname, "RECREATE");
  cout << "Writing output to file " << outfname <<endl;
  initPars(parfname);
// create and initialize objects we need
  rndgen = new TRandom3(0);	// start with unique seed every time
  fbeamatt = new TF1("fbeamatt","TMath::Exp(-x*[0])", 0, W);
  fbeamatt->SetParameter(0,Nsigma);
  hcce = new TH1D("hcce",Form("Sensor1D amplitudes, det W = %f cm",W), Nccebins, ccemin, ccemax);
  hcce->SetXTitle("Energy [MeV]");
  hcce->SetYTitle("dN/dE");
  fcce  = new TF1("fcce","[1]*(1-exp(-x/[1]))/[0] + [2]*(1-exp(-([0]-x)/[2]))/[0]",0,W);
  fcce->SetParameter(0,W); 
  fcce->SetParameter(1,mu_e*t_e*Efield); 
  fcce->SetParameter(2,mu_h*t_h*Efield); 
  initBragg(); 			// reads in Bragg curves for T and alphas
  cout << "Triton Bragg curve points: " << maxT << endl;
  cout << "Alpha Bragg curve points: "  << maxA << endl;
  if(debug>2) for (int i=0; i<10; i++)   braggT[i].print();
  if(debug>0) nt = new TNtuple("nt","Detailed event info","x:angle:cosangle:edep:cce");
  
  double x;   	// neutron interaction X coordinate
  double angle;	// decay angle of triton w.r.t X coordinate
  for (int i=1; i<=Nneutrons; i++) { 	// fire neutrons!
    double ccesum = 0, edep = 0;
    x = fbeamatt->GetRandom(); 		
    angle = rndgen->Rndm()*TMath::TwoPi(); // random angle
    const double cosangle = TMath::Cos(angle); // avoid recalculating cos 
    if(debug>2) cout << "i x angle cosangle    " << i<<" "<< x<<" "<<angle<<" "<<cosangle<<endl;
    // process triton response
    for(int j=0; j<maxT; j++) { 
      double xcosr = x + braggT[j].r * cosangle;
      if(xcosr>=0 && xcosr<=W) 
      { // triton is in the detector, add energy deposit 
	if(debug>0) edep += braggT[j].E;
	ccesum += braggT[j].E * fcce->Eval(xcosr);
      }
    }
    // process alpha response
    for(int j=0; j<maxA; j++) { 
      double xcosr = x - braggA[j].r * cosangle;
      if(xcosr>=0 && xcosr<=W) 
      { // alpha is in the detector, add energy deposit 
	if(debug>0) edep += braggA[j].E;
	ccesum += braggA[j].E * fcce->Eval(xcosr);
      }
    }
    hcce->Fill(ccesum);
    if(debug>0) nt->Fill(x,angle,cosangle,edep,ccesum);
    if(!(i%10000)) cout << " ... " << i <<endl;
  }
  cout << "All neutrons processed." << endl;
}


int main(int argc, char *argv[]) {
  long n = 10;
  char outfname[200] = "mydet.root";
  char parfname[200] = "mydet.par";
  if(argc>1) n = atol(argv[1]);
  if(argc>2) strcpy(outfname, argv[2]);
  if(argc==4) strcpy(parfname, argv[3]);
  if(argc>4) { 
    cout << "This program takes max 3 arguments: number of neutrons, the output file name, the parameter file name." << endl;
    return 99;
  }
  ccesensor1d(n, outfname, parfname);  
  fbeamatt->Write();
  fcce->Write();
  hcce->Write();
  if(nt) nt->Write();
  rout->Close();
  cout << "Output file written, finished!" << endl;
  return 0;
}
