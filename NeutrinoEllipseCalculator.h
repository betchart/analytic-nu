// C++ implementation of arXiv:1305.1872
// Analytic solutions for neutrino momenta in decay of top quarks
// Gala Nicolas Kaufman (gnn4@cornell.edu)

#include <memory>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>

#include <TMatrix.h>
#include <TMatrixD.h>
#include <TArrayD.h>
#include <TMath.h>
#include "Math/VectorUtil.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"

using namespace std;
using namespace ROOT::Math;

class NeutrinoEllipseCalculator{

 private:
  //inputs: lepton and b-jet momenta
  TLorentzVector bJet_;
  TLorentzVector lepton_;

  double   bJetBeta_,  bJetBeta2_,  bJetGamma_,  bJetGamma2_;
  double leptonBeta_,leptonBeta2_,leptonGamma_,leptonGamma2_;

  //particle masses
  double mW_, mt_, mnu_;
  double mW2_, mt2_, mnu2_;

  //numbers
  double x0_, x0p_;
  double Sx_, Sy_;
  double epsilon2_;

  double c_, s_; //cosine and sine of theta_{b,mu}

  double omega_;
  double Omega_;
  double x1_, y1_;
  double Z2_;

  //matrices
  TMatrixD Ab_;
  TMatrixD Al_;

  TMatrixD Htilde_;
  TMatrixD H_;
  TMatrixD Hperp_;

  TMatrixD Nperp_;

 public:
  NeutrinoEllipseCalculator();
  NeutrinoEllipseCalculator(TLorentzVector& , TLorentzVector& , double& , double& , double& );
  ~NeutrinoEllipseCalculator();

  void setbJet(const TLorentzVector& );
  void setLepton(const TLorentzVector& );

  void bJetRelativisticFactors();
  void leptonRelativisticFactors();

  void setTopMass(double& mTop){ mt_=mTop; };
  void setWBosonMass(double& mW){ mW_=mW; };
  void setNeutrinoMass(double& mNu){ mnu_=mNu; };
  void setMasses(double& , double& , double& );
  
  void fillAngles();

  void initializeMatrices();

  TMatrixD rotationMatrix(int, double);

  void Wsurface();
  void bJetEllipsoid();
  void leptonEllipsoid();
  void neutrinoSolution();
  void labSystemTransform();

  TMatrixD getNeutrinoEllipse();

};
