tbinSize = 3

hBSA =[]

for i in range(0,int(splitSize)):
    hBSA.append(ROOT.TGraphErrors())
#print(hBSA)


parsFit = np.zeros(int(splitSize))
parsFitError = np.zeros(int(splitSize))


count = 0
for i in BSA:
    tcount = 0
    for x in i:
        y = trePhi[tcount]
        yerr = eBSA[count][tcount]
        hBSA[count].SetPoint(hBSA[count].GetN(),y,x)
        #print(x,y,count,tCount)
        hBSA[count].SetPointError(hBSA[count].GetN()-1, 0, yerr)
        t1= ROOT.TF1("f1","[0]*sin(x*TMath::DegToRad())",0,360)
        #t1= ROOT.TF1("f1","([0]*sin(x*TMath::DegToRad()))/(1+[1]*cos(x*TMath::DegToRad())+[2]*cos(2*x*TMath::DegToRad()))",0,360)
        hBSA[count].Fit(t1,"QR")
        parsFit[count] = t1.GetParameter(0)
        parsFitError[count] = t1.GetParError(0)
        tcount += 1
    count += 1
        


c4 = ROOT.TCanvas("c1","c1", 1600,800)
c4.Divide(3,2,0,0)
c4.Draw()



count = 1
for i in range(0,6):
    c4.cd(i+1)
    title = "BSA Q2 bin " + str(count) + " -t bin " + str((i%tbinSize) +1)
    hBSA[i].SetTitle(title)
    hBSA[i].GetYaxis().SetRangeUser(-0.35,0.35)
    hBSA[i].SetMarkerStyle(20)
    hBSA[i].Draw("AP")
    if i % tbinSize == (tbinSize-1):
        count +=1
c2 = ROOT.TCanvas("c1","c1", 1600,800)
c2.Divide(3,3,0,0)
c2.Draw()



count = 3
for i in range(6,15):
    c2.cd(i-5)
    title = "BSA Q2 bin " + str(count) + " -t bin " + str((i%tbinSize) +1)
    hBSA[i].SetTitle(title)
    hBSA[i].GetYaxis().SetRangeUser(-0.35,0.35)
    hBSA[i].SetMarkerStyle(20)
    hBSA[i].Draw("AP")
    if i % tbinSize == (tbinSize-1):
        count +=1
    