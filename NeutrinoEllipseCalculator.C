#include "NeutrinoEllipseCalculator.h"

NeutrinoEllipseCalculator::NeutrinoEllipseCalculator()
{
  //cout << "constructor" << endl;
}

NeutrinoEllipseCalculator::~NeutrinoEllipseCalculator()
{
  //cout << "destructor" << endl;
}

NeutrinoEllipseCalculator::NeutrinoEllipseCalculator(TLorentzVector& bJetLorentzVector, 
						     TLorentzVector& leptonLorentzVector,
						     double& mTop, double& mW, double& mNu) :
  Ab_(4,4),
  Al_(4,4),
  Htilde_(3,3),
  H_(3,3),
  Hperp_(3,3),
  Nperp_(3,3)
{
  setbJet(bJetLorentzVector);
  setLepton(leptonLorentzVector);

  bJetRelativisticFactors();
  leptonRelativisticFactors();

  setMasses(mTop,mW,mNu);

  fillAngles();

  initializeMatrices();
}

void NeutrinoEllipseCalculator::setMasses(double& mTop, double& mW, double& mNu)
{
  setTopMass(mTop);
  setWBosonMass(mW);
  setNeutrinoMass(mNu);

  mW2_=mW_*mW_;
  mt2_=mt_*mt_;
  mnu2_=mnu_*mnu_;
}

void NeutrinoEllipseCalculator::setbJet(const TLorentzVector& bJetLorentzVector)
{
  bJet_=bJetLorentzVector;
}

void NeutrinoEllipseCalculator::setLepton(const TLorentzVector& leptonLorentzVector)
{
  lepton_=leptonLorentzVector;
}

void NeutrinoEllipseCalculator::bJetRelativisticFactors()
{
  bJetBeta_=bJet_.Beta();
  bJetBeta2_=bJetBeta_*bJetBeta_;
  bJetGamma_=bJet_.Gamma();
  bJetGamma2_=bJetGamma_*bJetGamma_;
}

void NeutrinoEllipseCalculator::leptonRelativisticFactors()
{
  leptonBeta_=lepton_.Beta();
  leptonBeta2_=leptonBeta_*leptonBeta_;
  leptonGamma_=lepton_.Gamma();
  leptonGamma2_=leptonGamma_*leptonGamma_;
}

void NeutrinoEllipseCalculator::fillAngles()
{
  c_=ROOT::Math::VectorUtil::CosTheta(lepton_,bJet_);
  s_=sqrt(1.-c_*c_);
}

void NeutrinoEllipseCalculator::initializeMatrices()
{
  Ab_.Zero();
  Al_.Zero();
  Htilde_.Zero();
  H_.Zero();
  Hperp_.Zero();
  Nperp_.Zero();
}

void NeutrinoEllipseCalculator::Wsurface()
{
  x0p_=-(0.5/bJet_.E())*(mt2_-mW2_-bJet_.M2());
  x0_=-(0.5/lepton_.E())*(mW2_-lepton_.M2()-mnu2_);
  Sx_=(1./leptonBeta2_)*(x0_*leptonBeta_-lepton_.P()*(1.-leptonBeta2_));
  epsilon2_=(1.-leptonBeta2_)*(mW2_-mnu2_);
}

void NeutrinoEllipseCalculator::bJetEllipsoid()
{
  Ab_[0][0]=1-c_*c_*bJetBeta2_;
  Ab_[1][0]=-c_*s_*bJetBeta2_;
  Ab_[2][0]=0;
  Ab_[3][0]=c_*x0p_*bJetBeta_;
    
  Ab_[0][1]=-c_*s_*bJetBeta2_;
  Ab_[1][1]=1-s_*s_*bJetBeta2_;
  Ab_[2][1]=0;
  Ab_[3][1]=s_*x0p_*bJetBeta_;
    
  Ab_[0][2]=0;
  Ab_[1][2]=0;
  Ab_[2][2]=1;
  Ab_[3][2]=0;
    
  Ab_[0][3]=c_*x0p_*bJetBeta_;
  Ab_[1][3]=s_*x0p_*bJetBeta_;
  Ab_[2][3]=0;
  Ab_[3][3]=mW2_-x0p_*x0p_;
}

void NeutrinoEllipseCalculator::leptonEllipsoid()
{
  Al_[0][0]=1.-leptonBeta2_;
  Al_[1][0]=0;
  Al_[2][0]=0;
  Al_[3][0]=Sx_*leptonBeta2_;
    
  Al_[0][1]=0;
  Al_[1][1]=1;
  Al_[2][1]=0;
  Al_[3][1]=0;
    
  Al_[0][2]=0;
  Al_[1][2]=0;
  Al_[2][2]=1;
  Al_[3][2]=0;
    
  Al_[0][3]=Sx_*leptonBeta2_;
  Al_[1][3]=0;
  Al_[2][3]=0;
  Al_[3][3]=mW2_-x0_*x0_-epsilon2_;
}

void NeutrinoEllipseCalculator::neutrinoSolution()
{
  Sy_=(1./s_)*(x0p_/bJetBeta_-c_*Sx_);
  omega_=(1./s_)*(leptonBeta_/bJetBeta_-c_); //only the positive slope
  Omega_=sqrt(max(0.,omega_*omega_+1.-leptonBeta2_));
  double Omega2=Omega_*Omega_;
  x1_=Sx_-(1./Omega2)*(Sx_+omega_*Sy_);
  y1_=Sy_-(1./Omega2)*omega_*(Sx_+omega_*Sy_);
  Z2_=x1_*x1_*Omega2-(Sy_-omega_*Sx_)*(Sy_-omega_*Sx_)-(mW2_-x0_*x0_-epsilon2_);
  double Z=sqrt(max(0.,Z2_));

  Htilde_[0][0]=Z/Omega_;
  Htilde_[0][1]=0;
  Htilde_[0][2]=x1_-lepton_.P();
	
  Htilde_[1][0]=omega_*Z/Omega_;
  Htilde_[1][1]=0;
  Htilde_[1][2]=y1_;
  	
  Htilde_[2][0]=0;
  Htilde_[2][1]=Z;
  Htilde_[2][2]=0;
}

TMatrixD NeutrinoEllipseCalculator::rotationMatrix(int axis, double angle)
{
  TMatrixD r(3,3);
  r.Zero();
  if (axis!=0 && axis!=1 && axis!=2) return r;
  
  for( int i=0; i<3; i++ ) 
    {
      r[i][i]=cos(angle);
    }

  for( int i=-1; i<=1; i++ )  
    {
      double row=(axis-i)%3; if(row<0) row+=3;
      double col=(axis+i)%3; if(col<0) col+=3;
      r[row][col]=i*sin(angle)+(1-i*i);
    }

  return r;
}

void NeutrinoEllipseCalculator::labSystemTransform()
{
  //rotate Htilde to H
  TMatrixD R(3,3);
  R.Zero();
  TMatrixD Rz=rotationMatrix(2,-lepton_.Phi());
  TMatrixD Ry=rotationMatrix(1,0.5*M_PI-lepton_.Theta());
  double bJetP[3]={bJet_.Px(),bJet_.Py(), bJet_.Pz()};
  TMatrixD bJet_xyz(3,1,bJetP);
  TMatrixD rM(Ry,TMatrixD::kMult,TMatrixD(Rz,TMatrixD::kMult,bJet_xyz));
  double* rA=rM.GetMatrixArray();
  double phi=-TMath::ATan2(rA[2],rA[1]);
  TMatrixD Rx=rotationMatrix(0,phi);
  R=TMatrixD(Rz,TMatrixD::kTransposeMult,TMatrixD(Ry,TMatrixD::kTransposeMult,Rx.T()));
  H_=TMatrixD(R,TMatrixD::kMult,Htilde_);

  //calculate Hperp
  double Hvalues[9]={H_[0][0],H_[0][1],H_[0][2],H_[1][0],H_[1][1],H_[1][2],0,0,1};
  TArrayD Harray(9,Hvalues);
  Hperp_.SetMatrixArray(Harray.GetArray());

  //calculate Nperp
  TMatrixD HperpInv(Hperp_);
  HperpInv.Invert();
  TMatrixD U(3,3);
  U.Zero();
  U[0][0]=1;
  U[1][1]=1;
  U[2][2]=-1;
  Nperp_=TMatrixD(HperpInv,TMatrixD::kTransposeMult,TMatrixD(U,TMatrixD::kMult,HperpInv));
}

TMatrixD NeutrinoEllipseCalculator::getNeutrinoEllipse()
{
  TMatrixD neutrinoEllipse(3,3); 
  neutrinoEllipse.Zero();
  Wsurface();
  leptonEllipsoid();
  bJetEllipsoid();
  neutrinoSolution();
  labSystemTransform();
  neutrinoEllipse=Nperp_;
  return neutrinoEllipse;
}
