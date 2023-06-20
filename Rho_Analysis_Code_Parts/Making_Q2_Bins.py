import numpy as np


nbins = 9# //number of bins for BSA

binSize = 360/nbins# //Size of each bin using Phi* (times 2 for Positive and negtive)

#hphistar = rdf.Histo1D(("hphistar","Trento #phi;Trento #phi",200,-200,200), "phi_t")
#cut2 = rdf.Filter("misse < 2")
cut2 = rdf.Filter("mmpro < 1.2 && mmpro > 0.85")
#cut2 = cut2.Filter("mt < 1")
#cut2 = cut2.Filter("mrho > 0.65 && mrho < 0.9")




cdNumber = 1
Np1 = np.zeros(nbins)
Nm1 = np.zeros(nbins)

q2Cuts = ["xb<x1 && Q2>2 && Q2<(y0 + (xb-x0)/(x1-x0)*(y1-y0)) && xb<x01","xb<x1 && Q2>2 && Q2<(y0 + (xb-x0)/(x1-x0)*(y1-y0)) && xb>x01","xb<x1 && Q2>2","xb>x1 && Q2>2 && Q2< (y4 + (xb-x4)/(x5-x4)*(y5-y4))","xb>x1 && Q2>2"]    

tLimits = [np.linspace(0.42,0.9,num=4),np.linspace(0.42,0.9,num=4),np.linspace(0.46,0.9,4),np.linspace(0.62,0.98,4),np.linspace(0.82,1.26,4)] 

tCuts = []
for i in tLimits:
    #print("-t > " + str(i[0]) + "  && -t < "+ str(i[1]))
    #print("-t > " + str(i[1]) + "  && -t < "+ str(i[2]))
    tCuts.append("mt > " + str(i[0]) + "  && mt < "+ str(i[1]))
    tCuts.append("mt > " + str(i[1]) + "  && mt < "+ str(i[2]))
    tCuts.append("mt > " + str(i[2]) + "  && mt < "+ str(i[3]))
    

phiBins = []


tCount = 0
#c1.Print(pdfname+"[")
for k in q2Cuts:
    print(k)
    for j in range(tCount,len(tCuts)):
        print(tCuts[j])
        new = cut2.Filter(k).Filter(tCuts[j])
        for i in range(0,360,int(binSize)):
            temp = int((i+binSize)/binSize - 1)
            printTemp = "phi_t > " + str(i) + " && phi_t < " + str(i+ binSize)
            #print(printTemp)
            c2 = new.Filter(printTemp)
            phiBins.append(c2)
        tCount += 1
        if tCount % 3 ==0:
            break
    
