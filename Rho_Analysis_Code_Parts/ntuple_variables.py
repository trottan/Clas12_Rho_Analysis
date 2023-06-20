import ROOT, sys
import numpy as np


ROOT.gStyle.SetGridColor(17)
ROOT.gStyle.SetPadGridX(1)
ROOT.gStyle.SetPadGridY(1)

ROOT.gStyle.SetPadRightMargin(0.01)
ROOT.gStyle.SetPadLeftMargin(0.075)
ROOT.gStyle.SetPadBottomMargin(0.12)

ROOT.gStyle.SetOptFit(1)

ROOT.TH1.AddDirectory(0)



fname_in =  "/work/clas12/kenjo/trains/root/eFPimFPip/lvl1_epimpip.nSidis_005*"
fname_out = "/work/clas12/trottan/rhoAnalysis/outb_data/epPipPim.outb.qa.skim_epimpip.nSidis_005*"

#fname = "outbending*root"
#fname = "outbending*root"
rdf = ROOT.RDataFrame("h22",{fname_in,fname_out})



vals = "pe,the,fie,ppip,thpip,fipip,ppim,thpim,fipim,angPimPip,phi_t"
vals += ",angPimPip_LF,angPimPip_CM"
vals += ",Q2,xb,ww,mt"
vals += ",mmpro,misse"
vals += ",mrho,rhoPhi_CM"
vals += ",x0,y0,x1,y1,x01,x4,y4,x5,y5"



rdf = rdf.Define("vals","""
double E0 = 10.6041, Mpro = 0.938272, Mele = 0.000511, Mpi = 0.13957, Mpi0 = 0.1349769;
TLorentzVector beam, targ;
beam.SetXYZM(0,0,E0,Mele);
targ.SetXYZM(0,0,0,Mpro);


TLorentzVector ele, pip, pim;
ele.SetXYZM(ex, ey, ez, Mele);
//pro.SetXYZM(prox, proy, proz, Mpro);
pip.SetXYZM(pipx, pipy, pipz, Mpi);
pim.SetXYZM(pimx, pimy, pimz, Mpi);

double pe = ele.P(), the=ele.Theta()*TMath::RadToDeg(), fie=ele.Phi()*TMath::RadToDeg();
double ppip = pip.P(), thpip=pip.Theta()*TMath::RadToDeg(), fipip=pip.Phi()*TMath::RadToDeg();
double ppim = pim.P(), thpim=pim.Theta()*TMath::RadToDeg(), fipim=pim.Phi()*TMath::RadToDeg();




double mrho = (pip+pim).M();
//double mdelta = (pro+pip).M();

//double mm2pip = (beam+targ-ele-pro-pim).M2();
//double mm2pim = (beam+targ-ele-pro-pip).M2();
double mmpro = (beam+targ-ele-pip-pim).M();
double misse = (beam+targ-ele-pip-pim).E();
//double mepx = (beam+targ-ele-pro).M();


auto q = (beam-ele);
double Q2 = - (beam-ele).M2();
double ww = (beam+targ-ele).M();
double xb = Q2/(ww*ww - targ.M2() + Q2);
double mt = - ((beam+targ-ele-pip-pim)-targ).M2();


double Theta_q = q.Theta();
double Phi_el = ele.Phi();

auto Rot_Matrix = [&](TLorentzVector vector, int Lab2CM_or_CM2Lab, double Theta_Rot, double Phi_Rot){ 

        double Rot_X1 = vector.X();
        double Rot_Y1 = vector.Y();
        double Rot_Z1 = vector.Z();

        double Rot_X = Rot_X1;
        double Rot_Y = Rot_Y1;
        double Rot_Z = Rot_Z1;



        // Lab2CM_or_CM2Lab is a parameter which determines if you rotating from the lab frame to the CM frame, or if you are rotating back in the opposite direction
        // Lab2CM_or_CM2Lab = -1 gives a rotation to the CM frame (from the lab frame)
        // Lab2CM_or_CM2Lab = +1 gives a rotation to the lab frame (from the CM frame)


        Theta_Rot = -1*Theta_Rot;   // Always give the angle of rotation Theta as the value given by .Theta()
                                    // This subroutine will handle the fact that the matrix rotation wants the negative of the angle of rotation


        // Rotation to Lab Frame
        if(Lab2CM_or_CM2Lab == -1){
            Rot_X = Rot_X1*TMath::Cos(Theta_Rot)*TMath::Cos(Phi_Rot) - Rot_Z1*TMath::Sin(Theta_Rot) + Rot_Y1*TMath::Cos(Theta_Rot)*TMath::Sin(Phi_Rot);
            Rot_Y = Rot_Y1*TMath::Cos(Phi_Rot) - Rot_X1*TMath::Sin(Phi_Rot);
            Rot_Z = Rot_Z1*TMath::Cos(Theta_Rot) + Rot_X1*TMath::Cos(Phi_Rot)*TMath::Sin(Theta_Rot) + Rot_Y1*TMath::Sin(Theta_Rot)*TMath::Sin(Phi_Rot);
        }


        // Rotation to CM Frame
        if(Lab2CM_or_CM2Lab == 1){
            Rot_X = Rot_X1*TMath::Cos(Theta_Rot)*TMath::Cos(Phi_Rot) + Rot_Z1*TMath::Cos(Phi_Rot)*TMath::Sin(Theta_Rot) - Rot_Y1*TMath::Sin(Phi_Rot);
            Rot_Y = Rot_Y1*TMath::Cos(Phi_Rot) + Rot_X1*TMath::Sin(Phi_Rot)*TMath::Cos(Theta_Rot) + Rot_Z1*TMath::Sin(Theta_Rot)*TMath::Sin(Phi_Rot);
            Rot_Z = Rot_Z1*TMath::Cos(Theta_Rot) - Rot_X1*TMath::Sin(Theta_Rot);
        }



        TLorentzVector vector_Rotated(Rot_X, Rot_Y, Rot_Z, vector.E());

        return vector_Rotated;


    };


auto angPimPip_LF = pim.Vect().Angle(pip.Vect())*TMath::RadToDeg();

auto rho = pip+pim;


// Center of Mass Frame
auto beam_Clone = Rot_Matrix(beam, -1, Theta_q, Phi_el);
auto targ_Clone = Rot_Matrix(targ, -1, Theta_q, Phi_el);
auto ele_Clone  = Rot_Matrix(ele,  -1, Theta_q, Phi_el);
auto pip_Clone = Rot_Matrix(pip,  -1, Theta_q, Phi_el);
auto pim_Clone = Rot_Matrix(pim,  -1, Theta_q, Phi_el);
auto q_Clone = Rot_Matrix(q,  -1, Theta_q, Phi_el);
auto rho_Clone = Rot_Matrix(rho,  -1, Theta_q, Phi_el);






auto fCM = q_Clone + targ_Clone;
auto boost = -(fCM.BoostVector());


auto qlv_Boost(q_Clone);
auto ele_Boost(ele_Clone);
auto beamBoost(beam_Clone);
auto targBoost(targ_Clone);
auto pipBoost(pip_Clone);
auto pimBoost(pim_Clone);
auto rhoBoost(rho_Clone);




qlv_Boost.Boost(boost);
ele_Boost.Boost(boost);
beamBoost.Boost(boost);
targBoost.Boost(boost);
pipBoost.Boost(boost);
pimBoost.Boost(boost);
rhoBoost.Boost(boost);





auto rhoPhi_CM = rhoBoost.Vect().Phi()*TMath::RadToDeg();



auto angPimPip_CM = pimBoost.Vect().Angle(pipBoost.Vect())*TMath::RadToDeg();




TVector3 v0, v1;
v0 = qlv_Boost.Vect().Cross(ele_Boost.Vect());
v1 = qlv_Boost.Vect().Cross(pipBoost.Vect()+pimBoost.Vect());
Double_t c0, c1, c2, c3;
c0 = v0.Dot(pipBoost.Vect()+pimBoost.Vect());
c1 = v0.Dot(v1);
c2 = v0.Mag();
c3 = v1.Mag();


double phi_t = (c0/TMath::Abs(c0)) * TMath::ACos(c1 /(c2*c3))*TMath::RadToDeg();
if(phi_t <= 0){
    phi_t = phi_t +360;
}








// Rho Frame
auto Phi_rho =  rhoBoost.Phi();
auto Theta_rho = rhoBoost.Theta();

TLorentzVector rhoBoost_r;
rhoBoost_r.SetXYZM(rhoBoost.Vect().x(), rhoBoost.Vect().y(), rhoBoost.Vect().z(), 0.7754);

auto pip_Clone2 = Rot_Matrix(pipBoost,  -1, Theta_rho, Phi_rho);
auto pim_Clone2 = Rot_Matrix(pimBoost,  -1, Theta_rho, Phi_rho);
auto rho_Clone2 = Rot_Matrix(rhoBoost_r,  -1, Theta_rho, Phi_rho);


auto fCM2 = -rho_Clone2;
auto boost2 = -(fCM2.BoostVector());

auto pipBoost2(pip_Clone2);
auto pimBoost2(pim_Clone2);

pipBoost2.Boost(boost2);
pimBoost2.Boost(boost2);
auto angPimPip = pimBoost2.Vect().Angle(pipBoost2.Vect())*TMath::RadToDeg();

double x0=0.15, y0=2.63;
double x1=0.45, y1=4.5;
double x01 = 0.34;
double x4=0.45, y4=4.8;
double x5=0.666, y5=6.2;

//ihel = static_cast<Int_t> (ihel);
//ihel = static_cast<Float_t> (ihel);

ihel = 1.00*ihel;

return vector<double> {"""+vals+"};")



phcuts = "angEPim>10 && angEPip>10"

for iv, vname in enumerate(vals.split(',')):
    rdf = rdf.Define(vname, "vals[{}]".format(iv))
    
    
    
    

    

    
    


    
    
rdf = rdf.Filter("Q2>0.8 && ww>1.8")
    



