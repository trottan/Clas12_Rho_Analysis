import numpy as np
ipar1 = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12], dtype=np.int32) 
ipar2 = np.array([13,1,2,14,4,5,6,7,8,15,16,11,12],dtype=np.int32)

#ipart = np.array([[0,1,2,3,4,5,6,7,8,9,10,11,12],[13,1,2,14,4,5,6,7,8,15,16,11,12]], dtype=np.int32) 



par_Num = []
for i in range(0,2*(135)):
    if i == 0:
        n = [0,1,2,3,4,5,6,7,8,9,10,11,12]
    else:
        j = 13 + (i-1)*4
        n = [j,1,2,j+1,4,5,6,7,8,j+2,j+3,11,12]
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
ROOT.gInterpreter.Declare("""
    ROOT::Math::Functor foo(const std::function<double(double const *)> &x, int num) { return ROOT::Math::Functor(x,num); }
    //ROOT::Math::Functor foo(GlobalChi2 x, int num) { return ROOT::Math::Functor(x,num); }
    
    """)


def fit1d(h):
   # 3 pieces to the over all fit
    fyp = []
    for i in range(0,len(h)):
        fyp.append(ROOT.TF1("fy","breitwigner(0)+breitwigner(3)+(x-2*0.1396)^[9]*exp(pol2(6))+breitwigner(10)",0.28,2.5))
        
        
    fy1p = []
    for i in range(0,len(h)):
        fy1p.append(ROOT.TF1("fy1","breitwigner(0)",0.28,2.5))
        
    fy2p = []
    for i in range(0,len(h)):
        fy2p.append(ROOT.TF1("fy2","breitwigner(0)",0.28,2.5))
        
    fy3p = []
    for i in range(0,len(h)):
        fy3p.append(ROOT.TF1("fy3","(x-2*0.1396)^[3]*exp(pol2(0))",0.28,2.5))
        
    fy4p = []
    for i in range(0,len(h)):
        fy4p.append(ROOT.TF1("fy4","breitwigner(0)",0.28,2.5))
        
        
    for i in fyp:
        i.SetParName(0,"#rhoAmp")
        i.SetParName(1,"#rho Mean")
        i.SetParName(2,"#rho Sigma")

        i.SetParName(3,"f2 Amp")
        i.SetParName(4,"f2 Mean")
        i.SetParName(5,"f2 Sigma")

        i.SetParName(6,"bg 0")
        i.SetParName(7,"bg 1")
        i.SetParName(8,"bg 2")
        i.SetParName(9,"bg n")
    
        i.SetParName(10,"f0 Amp")
        i.SetParName(11,"f0 Mean")
        i.SetParName(12,"f0 Sigma")
        
    
    opt = ROOT.Fit.DataOptions()
    rang = ROOT.Fit.DataRange()
   
    rang.SetRange(0.28,2.5)
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
    
    Npar = 13 + len(h)*4
    
    parTemp = [200,0.77549,0.1478,10,1.275,0.1867,1,0.05,-0.05,0.5,10,0.980,0.05,200,10,0.5,10]
    for  i in range(1,len(h)):
        parTemp.append(200)
        parTemp.append(10)
        parTemp.append(0.5)
        parTemp.append(10)
        
    print("Partemp len", len(parTemp), "Npar", Npar) 
    par0 = np.array(parTemp)
        

    
    fitter.Config().SetParamsSettings(Npar, par0)
    #fyp.FixParameter(4,1.275)
    #fyp.FixParameter(5,0.1867)
    #fitter.Config().ParSettings(4).Fix()
    fitter.Config().ParSettings(4).SetLimits(1,1.275)
    
    fitter.Config().ParSettings(5).SetLimits(0,0.1867)

     #fy.FixParameter(12,0.09)
    #fyp.SetParLimits(12,0.04,0.09)
    fitter.Config().ParSettings(12).SetLimits(0.04,0.09)
    
    #fyp.SetParLimits(2,0,1e7)
    #fitter.Config().ParSettings(2).Fix()
    fitter.Config().ParSettings(2).SetLimits(0,0.1478)

    #fyp.SetParLimits(3,0,1e7)
    fitter.Config().ParSettings(3).SetLimits(0,1e7)
        
    
    
    #fyp.SetParLimits(10,0,1e7)
    fitter.Config().ParSettings(10).SetLimits(0,1e7)
        
        
    #fyp.SetParLimits(11,0.8,1)
    fitter.Config().ParSettings(11).SetLimits(0.8,1)
    
    fitter.Config().ParSettings(0).SetLimits(0,1e7)
    
    
    for i in range(1,len(h)):
        j = 13 + (i-1)*4
        #print(j,j+1,j+2,j+3)
        fitter.Config().ParSettings(j).SetLimits(0,1e7)
        fitter.Config().ParSettings(j+1).SetLimits(0,1e7)
        fitter.Config().ParSettings(j+2).SetLimits(0,1)
        fitter.Config().ParSettings(j+3).SetLimits(0,1e7)
        
    
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
    
    print("parsN ", parsN[0],parsN[13],parsN[17])
    print("parsNError ", parsNError[0],parsNError[13],parsNError[17])
    
    for i in range(0,len(h)):
        fyp[i].SetFitResult(result, ipart[i])
        fyp[i].SetRange(rang().first, rang().second)
        
    
    rhoAmp = [parsN[0]]
    for i in range(1,len(h)):
        j = 13 + (i-1)*4
        rhoAmp.append(parsN[j])
        
    rhoAmpError = [parsNError[0]]
    for i in range(1,len(h)):
        j = 13 + (i-1)*4
        rhoAmpError.append(parsNError[j])
        
        
    for i in range(0,len(h)):
        if i == 0:
            fy1p[i].FixParameter(0,parsN[0])
            fy1p[i].FixParameter(1,parsN[1])
            fy1p[i].FixParameter(2,parsN[2])
        else:
            j = 13 + (i-1)*4
            fy1p[i].FixParameter(0,parsN[j])
            fy1p[i].FixParameter(1,parsN[1])
            fy1p[i].FixParameter(2,parsN[2])
        fy1p[i].SetLineColor(ROOT.kBlack)
        
        
    for i in range(0,len(h)):
        if i == 0:
            fy2p[i].FixParameter(0,parsN[3])
            fy2p[i].FixParameter(1,parsN[4])
            fy2p[i].FixParameter(2,parsN[5])
        else:
            j = 13 + (i-1)*4 +1
            fy2p[i].FixParameter(0,parsN[j])
            fy2p[i].FixParameter(1,parsN[4])
            fy2p[i].FixParameter(2,parsN[5])
        fy2p[i].SetLineColor(ROOT.kGreen)
        
        
    for i in range(0,len(h)):
        if i == 0:
            fy3p[i].FixParameter(0,parsN[6])
            fy3p[i].FixParameter(1,parsN[7])
            fy3p[i].FixParameter(2,parsN[8])
            fy3p[i].FixParameter(3,parsN[9])
        else:
            j = 13 + (i-1)*4 +2
            fy3p[i].FixParameter(0,parsN[6])
            fy3p[i].FixParameter(1,parsN[7])
            fy3p[i].FixParameter(2,parsN[8])
            fy3p[i].FixParameter(3,parsN[j])
        fy3p[i].SetLineColor(ROOT.kBlue)
        
    for i in range(0,len(h)):
        if i == 0:
            fy4p[i].FixParameter(0,parsN[10])
            fy4p[i].FixParameter(1,parsN[11])
            fy4p[i].FixParameter(2,parsN[12])
        else:
            j = 13 + (i-1)*4 +3
            fy4p[i].FixParameter(0,parsN[j])
            fy4p[i].FixParameter(1,parsN[11])
            fy4p[i].FixParameter(2,parsN[12])
        fy4p[i].SetLineColor(ROOT.kPink)

    for i in range(0,len(h)):
        hp[i].GetListOfFunctions().Add(fyp[i])
        hp[i].GetListOfFunctions().Add(fy1p[i])
        hp[i].GetListOfFunctions().Add(fy2p[i])
        hp[i].GetListOfFunctions().Add(fy3p[i])
        hp[i].GetListOfFunctions().Add(fy4p[i])
    print(len(hp),hp)
    return hp,rhoAmp, rhoAmpError