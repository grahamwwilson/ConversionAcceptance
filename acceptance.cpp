// Conversion geometric acceptance calculation.
//
// Here assume primary vertex is a Gaussian in z with specified (x,y,z).
// Extrapolate photon in a straight line and see if it 
// satisfies the |z| < 25 cm cut at the specified radius. 
// here hard-coded as R=16.2 cm. 
//
// If so histogram the corresponding 1/sin(theta) factor and 
// the resulting z_PC distribution.
// 
// Code currently assumes that the photon eta distribution is 
// uniform. This is a good first approximation but in principle should 
// be refined a bit.
// 
//              Graham W. Wilson, June 15th, 2020.
//
#include <iostream> 
#include <algorithm> //std::sort
#include <cmath>
#include <TRandom3.h>
#include <TMath.h>   //TMath::Prob
#include <TH1D.h>
#include <TFile.h>
#include <vector>    
#include <cstdlib>
//using namespace std;

int main(){

   double PIE=4.0*atan(1.0);
 
   const double ZPV = -0.58;
   const double ZPVRES = 3.394;
//   const double ETAMAX = 4.0;
   const double ETAMAX = -log(tan(0.5*acos(0.85)));
   const double X0 =  0.096;
   const double Y0 = -0.064 ;

   unsigned int seed = 4359;
   std::cout << "Using seed " << seed << std::endl;
   std::cout << "Using ETAMAX " << ETAMAX << std::endl;

   const unsigned int NTHROWS = 1000000;
   std::cout << "Number of test particles " << NTHROWS << std::endl;
   std::vector<double> v;

   TRandom3 *rg1 = new TRandom3(seed);
   TRandom3 *rg2 = new TRandom3(seed+1);
   TRandom3 *rg3 = new TRandom3(seed+2);
   TFile *f = new TFile("acceptance.root","RECREATE");
   TH1D* etahist = new TH1D("etahist","etahist",150,-5.0,5.0);
   TH1D* phihist = new TH1D("phihist","phihist",100,0.0,2.0*PIE);
   TH1D* zpvhist = new TH1D("zpvhist","zpvhist",100,-20.0,20.0);
   TH1D* zhist   = new TH1D("zhist","zhist",200,-100.0,100.0);
   TH1D* factorhist = new TH1D("factorhist","factorhist",201,0.9975,2.0025);

// Simulate track acceptance

   const double R = 16.2;
//   const double R = 11.0;
//   const double R = 2.9;

   for (int i=0; i<NTHROWS; i++){
      double zpv = rg1->Gaus(ZPV,ZPVRES);
      zpvhist->Fill(zpv);
      double eta = rg2->Uniform(-ETAMAX,ETAMAX);
      etahist->Fill(eta);
      double phi = rg3->Uniform(0.0,2.0*PIE);
      phihist->Fill(phi);
// Calculation time
      double theta = 2.0*atan(exp(-eta));
// Starting position is (X0, Y0, zpv)
// direction cosines for each coordinate
      double vx = sin(theta)*cos(phi);
      double vy = sin(theta)*sin(phi);
      double vz = cos(theta);
// Need to solve a quadratic equation for photon straight-line propagation
      double a = vx*vx + vy*vy;
      double b = 2.0*(vx*X0 + vy*Y0);
      double c = X0*X0 + Y0*Y0 - R*R;
 
// There are two solutions, one is in the same direction as phi that is +ve, 
// and the other is opposite and is -ve.
      double det = b*b-4.0*a*c;
      double Lplus =  (-b + sqrt(det))/(2.0*a);
      double Lminus = (-b - sqrt(det))/(2.0*a);
      double L = std::max(Lplus, Lminus);
      double xp = X0 + L*vx;
      double yp = Y0 + L*vy;
      double zp = zpv + L*vz;
      zhist->Fill(zp);
      if(i<100){
         std::cout << " i: " << i << " zpv: " << zpv << " phi: " << phi << " eta: " << eta 
                   << " theta: " << theta << " costh " << cos(theta) 
                   << " 1/sinth " << 1.0/sin(theta) << " L+: " << Lplus << " " << " L-:  " << Lminus << std::endl;
         std::cout << "Extrapolated to " << xp << " " << " " << yp << " " << zp << " with R = " << sqrt(xp*xp+yp*yp) << std::endl;
         std::cout << " " << std::endl;
      }
      if(abs(zp)<25.0)factorhist->Fill(1.0/sin(theta));
   }
   zpvhist->Draw();
   etahist->Draw();
   phihist->Draw();
   f->Write();
}
