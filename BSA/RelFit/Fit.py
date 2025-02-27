import numpy as np
import ROOT
#Function to have the same mean and sig values for each fit function
ipar1 = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12], dtype=np.int32) 
ipar2 = np.array([13,1,2,14,4,5,6,7,8,15,16,11,12],dtype=np.int32)


#ipart = np.array([[0,1,2,3,4,5,6,7,8,9,10,11,12],[13,1,2,14,4,5,6,7,8,15,16,11,12]], dtype=np.int32) 

par_Num = []
for i in range(0,2):
    if i == 0:
        n = [0,1,2,3,4,5,6,7,8,9]
    else:
        j = 10 + (i-1)*3
        #n = [j,1,2,j+1,4,5,j+2,7,8,j+3,10,11]
        n = [j,1,2,j+1,4,5,6,j+2,8,9]
    #print(n)
    par_Num.append(n)
print((par_Num))
ipart = np.array(par_Num, dtype=np.int32) 


#Custum class that makes a chi2 for both fits

class GlobalChi2(object):
    def __init__(self, f1,numOfHist):
        self._f1 = f1
        #self._f2 = f2
        self._numOfHist = numOfHist
    def __call__(self, par):
        # parameter vector is first background (in common 1 and 2) and then is
        # signal (only in 2)
 
        # the zero-copy way to get a numpy array from a double *
        par_arr = np.frombuffer(par, dtype=np.float64, count=self._numOfHist*len(ipar1))
 

        
        #p1 = par_arr[ipar1]
        #p2 = par_arr[ipar2]
        #p = [p1,p2]
        
        p = []
        for i in range(0,self._numOfHist):
            #print(par_arr[i])
            p.append(par_arr[ipart[i]])
        
        
        tot = 0
        #for i in self._f1:
        #    tot += i(p1)
        for i in range(0,(self._numOfHist)):
            tot += self._f1[i](p[i])

        
        return tot

# Define custom root fit functions for the Relativisit Breit Wigner and Background
ROOT.gInterpreter.Declare("""
    #include <cmath>
            
            
    #include "TMath.h"
    ROOT::Math::Functor foo(const std::function<double(double const *)> &x, int num) { return ROOT::Math::Functor(x,num); }
    //ROOT::Math::Functor foo(GlobalChi2 x, int num) { return ROOT::Math::Functor(x,num); }


        class FunctionSet {
            public:
            double breit_wigner_1(double *x, double *par) {
                double m = x[0];
                double amp = par[0];
                double m0 = par[1];
                double G0 = par[2];
                double mpi = 0.139; // Mass of pi+
    
                double p = TMath::Sqrt(m*m-4*mpi*mpi)/2;
                double p0 = TMath::Sqrt(m0*m0-4*mpi*mpi)/2;
                double Gamma = G0*TMath::Power(p / p0, 3);
                double gamma = TMath::Sqrt(m*m*(m0*m0 + Gamma*Gamma));
    
                double k1 = 2 * TMath::Sqrt(2) *m0* Gamma * gamma;
                double k2 = TMath::Pi() * TMath::Sqrt(m0*m0 + gamma);
                double k = k1 / k2;
    
                double bw1 = amp * k / ((m*m - m0*m0)*(m*m - m0*m0) + (Gamma*Gamma*m0*m0));
                return bw1;
            }
    
                double breit_wigner_2(double *x, double *par) {
                    double m = x[0];
                    double amp = par[7];
                    double mu = par[8];
                    double gamma = par[9];
    
                    return amp*TMath::BreitWigner(m, mu, gamma);
    
    
                }
    
            double background(double *x, double *par) {
                double m = x[0];
                double amp_bkg = par[3];
                double xi = par[4];
                double m0_bkg = par[5];
                double m1_bkg = par[6];
    
                double A1 = m1_bkg-m;
                double A2 = (m1_bkg-m0_bkg)*(m1_bkg-m0_bkg);
                double A = A1/A2;
    
                double B1 = (m1_bkg-m)*(m1_bkg-m);
                double B2 = (m1_bkg-m0_bkg)*(m1_bkg-m0_bkg);
                double B = B1/B2;
    
                double C1 = (m1_bkg-m)*(m1_bkg-m);
                double C2 = (m1_bkg-m0_bkg)*(m1_bkg-m0_bkg);
                double C = C1/C2;
    
                double func;
                //if(B < 1){
                    func =xi*xi*xi*A*TMath::Sqrt(1-B)*TMath::Exp(-1/2*xi*xi*(1-C));
                //} else {
                //    func = 0;
                //}
                return  amp_bkg*func;
            }
    
            double background_only(double *x, double *par) {
                double m = x[0];
                double amp_bkg = par[0];
                double xi = par[1];
                double m0_bkg = par[2];
                double m1_bkg = par[3];
    
                double A1 = m1_bkg-m;
                double A2 = (m1_bkg-m0_bkg)*(m1_bkg-m0_bkg);
                double A = A1/A2;
    
                double B1 = (m1_bkg-m)*(m1_bkg-m);
                double B2 = (m1_bkg-m0_bkg)*(m1_bkg-m0_bkg);
                double B = B1/B2;
    
                double C1 = (m1_bkg-m)*(m1_bkg-m);
                double C2 = (m1_bkg-m0_bkg)*(m1_bkg-m0_bkg);
                double C = C1/C2;
    
                double func;
                //if(B < 1){
                    func = xi*xi*xi*A*TMath::Sqrt(1-B)*TMath::Exp(-1/2*xi*xi*(1-C));
                //} else {
                //    func = 0;
                //}
                return amp_bkg*func;
            }
    
    
    
            double total_function(double *x, double *par) {
                double m = x[0];
                double amp1 = par[0];
                double m01 = par[1];
                double G01 = par[2];
    
    
                double amp_bkg = par[3];
                double xi = par[4];
                double m0_bkg = par[5];
                double m1_bkg = par[6];
    
                double amp2 = par[7];
                double m02 = par[8];
                double G02 = par[9];
    
                // Evaluate the first Breit-Wigner function
                double bw1 = breit_wigner_1(x, par);
    
                double bkg = background(x,par);
    
                double bw2 = breit_wigner_2(x, par);
    
                //std::cout << bw1 << " " << bkg << " " << bw2 << endl;
                // Return the sum of the two Breit-Wigner functions and the background
                return bw1+ bkg + bw2;
                }
            };



    
    
    """)

#Fit function Q2binNum = number for Q2bin, tbinNum = number for t bin, phiCount = phi bin number, maxValue = max value invariant mass in hist
def fit1d(h,Q2binNum,tbinNum,phiCount,maxValue):
   # 3 pieces to the over all fit
    fyp = []

    print("number of hist",len(h))
    # Create the TF1 object using the defined function
    #custom_tf1 = ROOT.TF1("custom_tf1", cppyy.gbl.breit_wigner_1, 0.3, 2.5, 3)
    funcs = cppyy.gbl.FunctionSet()


    for i in range(0,len(h)):
        fyp.append(ROOT.TF1("fy", funcs.total_function,0.3,maxValue,10))
        
        
    fy1p = []
    for i in range(0,len(h)):
        fy1p.append(ROOT.TF1("fy",funcs.breit_wigner_1,0.3,maxValue,3))
  
    fy2p = []
    for i in range(0,len(h)):
        fy2p.append(ROOT.TF1("fy",funcs.background_only,0.3,maxValue,4))
    fy3p = []
    for i in range(0,len(h)):
        fy3p.append(ROOT.TF1("fy","breitwigner(0)",0.3,maxValue))
        
    for i in fyp:
        i.SetParName(0,"#rho Amp")
        i.SetParName(1,"#rho Mean")
        i.SetParName(2,"#rho Sigma")
        
        i.SetParName(3,"bkg Amp")
        i.SetParName(4,"bkg xi")
        i.SetParName(5,"bkg m0")
        i.SetParName(6,"bkg m1")
        
        i.SetParName(7,"f2 Amp")
        i.SetParName(8,"f2 Mean")
        i.SetParName(9,"f2 Sigma")
        
    
    opt = ROOT.Fit.DataOptions()
    rang = ROOT.Fit.DataRange()

    rang.SetRange(0.3,maxValue)
    data = [None]*len(h)
    chi2 = []
    
    
    hp = []
    for i in range(0,len(h)):
        hp.append(ROOT.TH1D())
        hp[i] = h[i].GetPtr()
        
        wfy1 = ROOT.Math.WrappedMultiTF1(fyp[i], 1)
        data[i] = ROOT.Fit.BinData(opt,rang)
        ROOT.Fit.FillData(data[i],hp[i])
        
        chi2.append(ROOT.Fit.Chi2Function(data[i],wfy1))
       
    
    
    print("dog",chi2[0],chi2[0].NDim())
    
    
    #Making global Chi2 in custum class
    globalChi2 = GlobalChi2(chi2,len(h))
    
    
    
    
    fitter = ROOT.Fit.Fitter()
    
    Npar = 10 + (len(h)-1)*3
    
        
    #Using the fits from normal BW; parCount access the phi bin for each kinematic bin
    if Q2binNum == 0 and tbinNum == 3:
        parCount = 20
    elif Q2binNum == 1 and tbinNum == 0 and (phiCount == 5 or phiCount ==6):
        parCount = 4*20
    elif Q2binNum == 2 and (tbinNum == 1) and phiCount != 0:
        parCount = 0
    elif  Q2binNum == 2 and (tbinNum == 1) and phiCount == 0:
        parCount = 0
    elif Q2binNum == 2 and (tbinNum == 0):
        parCount = 5
    elif Q2binNum == 3 and tbinNum == 0 and (phiCount == 4 or phiCount ==5):
        parCount = 3*20
    elif Q2binNum == 4:# and tbinNum == 0 and (phiCount == 3):
        parCount = 1*20
    txtFile = "/w/hallb-scshelf2102/clas12/trottan/rhoAnalysis/rhoFits_andrey/Fits.fitted.hists.mmp-lt-110.bin" + str(Q2binNum) + str(tbinNum) + ".root.txt"
    print(txtFile)
    fitPar = []
   
    with open(txtFile,'r') as file_temp:
        line = file_temp.readlines()
        for i in range(parCount,parCount + 20):
            #print(float(line[i]))
            fitPar.append(float(line[i]))
            parCount +=1
    print(fitPar)
    
    parTemp = []
    
    
   #parTemp is my initialization of parameters
    parTemp.append(fitPar[0])
    parTemp.append(fitPar[1])
    parTemp.append(fitPar[2])
    
    if (Q2binNum == 3 and tbinNum == 1):
        parTemp.append(0)
    # elif (Q2binNum == 4):
    #     parTemp.append(40000)
    elif (Q2binNum == 2 and tbinNum == 0 and phiCount ==4):
        parTemp.append(1)
    elif (Q2binNum ==4 and tbinNum ==0):
        parTemp.append(fitPar[6]*5000)
    else:
        parTemp.append(fitPar[6])
    parTemp.append(fitPar[7])
    parTemp.append(fitPar[8])
    parTemp.append(fitPar[9])
    
    
    if (Q2binNum == 2 and tbinNum == 2 and phiCount != 0):
        parTemp.append(0)
    else:
        parTemp.append(fitPar[3])
    parTemp.append(fitPar[4])
    parTemp.append(fitPar[5])
    

    parTemp.append(fitPar[10])
    if (Q2binNum == 3 and tbinNum == 1) or (Q2binNum == 2 and tbinNum == 0 and phiCount ==4):
        parTemp.append(0)
    elif (Q2binNum == 2 and tbinNum == 0 and phiCount ==4):
        parTemp.append(1)
    elif(Q2binNum == 4 and tbinNum == 0):
        parTemp.append(fitPar[16]*5000)
    else:
        parTemp.append(fitPar[16])

    if Q2binNum == 2 and tbinNum == 2:
        parTemp.append(0)
    else:
        parTemp.append(fitPar[13])
    
    
    
    
 
    
#     parTemp = fitPar
#     parTemp = [2000, 0.777,0.15]
#     parTemp += [10,0.01,0.3,31]
#     parTemp += [0,1.275,0.3]
#     parTemp += [200,10,0]
    
    print(parTemp)
    #parTemp += [1,0.05,-0.05,0.5]
    #parTemp += [100,0.05,0.5]

    
    #fyp.SetParLimits(3,0,1e7)

        
    print("Partemp len", len(parTemp), "Npar", Npar) 
    par0 = np.array(parTemp)
    fitter.Config().SetParamsSettings(Npar, par0)
    
    #Here I am limiting the parameters
    fitter.Config().ParSettings(0).SetLimits(0,1e7) # 0,1,2
    if Q2binNum == 4:
        fitter.Config().ParSettings(1).SetLimits(0,1)
        fitter.Config().ParSettings(2).SetLimits(0.05,0.175)
    else:
        fitter.Config().ParSettings(1).SetLimits(0,10)
        fitter.Config().ParSettings(2).SetLimits(0.05,0.175)


    if Q2binNum == 4 and tbinNum == 0:
        print("Setting specfic background values")
        fitter.Config().ParSettings(3).SetLimits(0,1e7) # 3,4,5,6
        fitter.Config().ParSettings(4).SetLimits(6,10) #(1,2) # 3,4,5,6
        fitter.Config().ParSettings(5).SetLimits(2,10)#SetLimits(0.2,0.35) #.Fix() # 3,4,5,6
        fitter.Config().ParSettings(6).SetLimits(0,2)#SetLimits(1,10) # 3,4,5,6fitter.Config().ParSettings(5)
    else:
        fitter.Config().ParSettings(3).SetLimits(0,1e7) # 3,4,5,6
        fitter.Config().ParSettings(4).SetLimits(0,1e7) # 3,4,5,6
        fitter.Config().ParSettings(5).SetLimits(0,1) #.Fix() # 3,4,5,6
        fitter.Config().ParSettings(6).SetLimits(0,5) # 3,4,5,6fitter.Config().ParSettings(5)
    #Fix values if zero
    if fitPar[6] ==0:
        fitter.Config().ParSettings(3).Fix()
        fitter.Config().ParSettings(4).Fix()
        fitter.Config().ParSettings(5).Fix()
        fitter.Config().ParSettings(6).Fix()

    fixValues =  (Q2binNum == 2 and tbinNum == 0 and phiCount != 2) or (Q2binNum == 2 and tbinNum == 2) or (Q2binNum ==4 and tbinNum == 0)
    #fixValues += 
    if  (fitPar[3] ==0) or (Q2binNum == 2 and tbinNum == 2):
        print("fitPar == 3 and ",Q2binNum,tbinNum,phiCount)
        if fixValues:
            print("FIXING F2 values to zero")
            fitter.Config().ParSettings(7).Fix()
            fitter.Config().ParSettings(8).Fix()
            fitter.Config().ParSettings(9).Fix()
            fitter.Config().ParSettings(12).Fix()
        else:
            print("NOT FIXING F2 values to zero")
            fitter.Config().ParSettings(7).SetLimits(0,1e7) # 7,8,9
            fitter.Config().ParSettings(12).SetLimits(0,1e7) # 12,8,9
            if Q2binNum == 2 and (tbinNum == 1) and phiCount == 0:
                fitter.Config().ParSettings(8).SetLimits(0,1)
            else:
                 fitter.Config().ParSettings(8).SetLimits(0,2)
            fitter.Config().ParSettings(9).SetLimits(0,2)
    else:
        print("Skipping FIXING F2 values to zero")
        fitter.Config().ParSettings(7).SetLimits(0,1e7) # 7,8,9 
        fitter.Config().ParSettings(12).SetLimits(0,1e7) # 7,8,9 
        if (Q2binNum == 3 and tbinNum == 1):
            print("3 1 for f2 fix")
            fitter.Config().ParSettings(8).SetLimits(1,1.5)
            fitter.Config().ParSettings(9).SetLimits(0.1,0.5)
        else:
            fitter.Config().ParSettings(8).SetLimits(1,maxValue)
            fitter.Config().ParSettings(9).SetLimits(0,0.2)
    # if fitPar[3] ==0 or parTemp[6] == 0:
    #     fitter.Config().ParSettings(7).Fix()
    #     fitter.Config().ParSettings(8).Fix()
    #     fitter.Config().ParSettings(9).Fix()
    
    
    fitter.Config().ParSettings(10).SetLimits(0,1e7)
    fitter.Config().ParSettings(11).SetLimits(0,1e7) # 7,8,9
    if fitPar[6] ==0:
        fitter.Config().ParSettings(11).Fix()
        
    
    
    
    
    
    #Minimizing chi2 fits
    fitter.Config().MinimizerOptions().SetPrintLevel(0)
    fitter.Config().SetMinimizer("Minuit2", "Migrad")
    
    
    #in order to fit global Chi2 function it needs to be in a wrapper
    #Functor fails in pyroot for root version 6.24 or less
    #globalChi2Functor = ROOT.Math.Functor(globalChi2, Npar) Faster way of doing it for v6.26 and up 
    
    
    #calls ROOT functor in ROOT
    globalChi2Functor = ROOT.foo(globalChi2,Npar)
    
    print("functor ",globalChi2Functor.NDim())
    
    #making fit for our Global Chi2 function
    dataSum = 0
    for i in data:
        dataSum += i.Size()
    print(dataSum)

    fitter.FitFCN(globalChi2Functor, 0,int(dataSum), True)
    
 
    
    
    result = fitter.Result()
    parsN = result.GetParams()
    parsNError = result.GetErrors()
    print(parsN)
    print("old Fit", fitPar[0],fitPar[6],fitPar[3],fitPar[10],fitPar[16],fitPar[13])
    print("parsN ", parsN[0],parsN[3],parsN[7],parsN[10],parsN[11],parsN[12])
    print("parsNError ", parsNError[0],parsNError[3],parsNError[7],parsNError[10],parsNError[11],parsNError[12])
    #print("parsNError ", parsNError[0],parsNError[13],parsNError[17])
    
    for i in range(0,len(h)):
        fyp[i].SetFitResult(result, ipart[i])
        fyp[i].SetRange(rang().first, rang().second)
        
    
#     rhoAmp = [parsN[0]]
#     for i in range(1,len(h)):
#         j = 13 + (i-1)*4
#         rhoAmp.append(parsN[j])
        
#     rhoAmpError = [parsNError[0]]
#     for i in range(1,len(h)):
#         j = 13 + (i-1)*4
#         rhoAmpError.append(parsNError[j])
        
    #Fixes values for individual fitting functions.     
    for i in range(0,len(h)):
        if i == 0:
            fy1p[i].FixParameter(0,parsN[0])
            fy1p[i].FixParameter(1,parsN[1])
            fy1p[i].FixParameter(2,parsN[2])
        else:
            j = 10 + (i-1)*3
            fy1p[i].FixParameter(0,parsN[j])
            fy1p[i].FixParameter(1,parsN[1])
            fy1p[i].FixParameter(2,parsN[2])
        fy1p[i].SetLineColor(ROOT.kBlack)
        
    for i in range(0,len(h)):
        if i == 0:
            fy2p[i].FixParameter(0,parsN[3])
            fy2p[i].FixParameter(1,parsN[4])
            fy2p[i].FixParameter(2,parsN[5])
            fy2p[i].FixParameter(3,parsN[6])
        else:
            j = 10 + (i-1)*3
            fy2p[i].FixParameter(0,parsN[j+1])
            fy2p[i].FixParameter(1,parsN[4])
            fy2p[i].FixParameter(2,parsN[5])
            fy2p[i].FixParameter(3,parsN[6])
        fy2p[i].SetLineColor(ROOT.kGreen)
            
    for i in range(0,len(h)):
        if i == 0:
            fy3p[i].FixParameter(0,parsN[7])
            fy3p[i].FixParameter(1,parsN[8])
            fy3p[i].FixParameter(2,parsN[9])
        else:
            j = 10 + (i-1)*3
            fy3p[i].FixParameter(0,parsN[j+2])
            fy3p[i].FixParameter(1,parsN[8])
            fy3p[i].FixParameter(2,parsN[9])
        fy3p[i].SetLineColor(ROOT.kBlue)

    # for i in range(0,len(h)):
    #     hp[i].GetListOfFunctions().Add(fyp[i])
    #     hp[i].GetListOfFunctions().Add(fy1p[i])
    #     hp[i].GetListOfFunctions().Add(fy2p[i])
    #     hp[i].GetListOfFunctions().Add(fy3p[i])
        
        
    print(len(hp),hp)
    return parsN[0],parsN[10],parsNError[0],parsNError[10],fyp,fy1p,fy2p, h#,rhoAmp, rhoAmpErro

