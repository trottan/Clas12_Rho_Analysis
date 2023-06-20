tVals = [0.4971369440598597, 0.6569550834832922, 0.8170997962207208, 0.5030596180913797, 0.6596071029685012, 0.8182945024243566, 0.5325006497950542, 0.6784170258792163, 0.8247163870744437, 0.6830860790677337, 0.8015944237187014, 0.9209668879406028, 0.8951726733776048, 1.0407006380693573, 1.186574493472095]

hWSPhi = []
for i in range(0,(6)):
    hWSPhi.append(ROOT.TGraphErrors())
    
    
for i in range(0,5):
    count = 0
    for x in parsFit[i]:
        y =  tVals[i][count]
        yerr = parsFitError[i][count]
        hWSPhi[i].SetPoint(hWSPhi[i].GetN(),y,x)
        print(x,y,count)
        hWSPhi[i].SetPointError(hWSPhi[i].GetN()-1, 0, yerr)
        count = count + 1
        t1= ROOT.TF1("f1","pol2(0)",0,360)
        t1.SetParameter(0,1e5)
        #hWSPhi.Fit(t1,"QR")
c4 = ROOT.TCanvas("c1","c1", 1600,800)
c4.Divide(2,3,0,0)
c4.Draw()

ROOT.gStyle.SetTitleFontSize(0.05)

hWSPhi[5] = hWSPhi[4].Clone()

legend = ROOT.TLegend(0.05, 0.35, 0.95, 0.2)
legend.SetNColumns(4)
legend.AddEntry(hWSPhi[0],"Combined (w/o low rho fits)","P")



for i in range(0,6):
    c4.cd(i+1)
    hWSPhi[i].SetTitle("A^{sin(#phi)}_{LU} vs -t (Q^{2} bin" + str(i+1)+")")

    hWSPhi[i].GetXaxis().SetNdivisions(5)
    hWSPhi[i].GetYaxis().SetNdivisions(5)
    #hWSPhi[i].SetMarkerSize(2)
    
    hWSPhi[i].GetXaxis().SetLimits(0.2,1.2)
    hWSPhi[i].GetXaxis().SetTitle("-t")


    hWSPhi[i].GetXaxis().SetLabelSize(0.05)
    hWSPhi[i].GetXaxis().SetTitleSize(0.05)
    hWSPhi[i].GetXaxis().SetTitleOffset(0.9)

    hWSPhi[i].GetYaxis().SetRangeUser(-0.1,0.2)
    hWSPhi[i].GetYaxis().SetTitle("A^{sin(#phi)}_{LU}")
    hWSPhi[i].GetYaxis().SetLabelSize(0.05)
    hWSPhi[i].GetYaxis().SetTitleSize(0.05)
    hWSPhi[i].GetYaxis().SetTitleOffset(0.8)

    
    hWSPhi[i].SetMarkerStyle(8)
    hWSPhi[i].SetMarkerSize(1.5)
    hWSPhi[i].Draw("AP")


c4.cd(6)
hWSPhi[5].Draw("AF")
hWSPhi[i].SetTitle()
legend.Draw()