canvas = ROOT.TCanvas("canvas")
pdf_file = ROOT.TPDF("Q2xB_IMplotsCombine.pdf", 0)


canvas.SaveAs(pdf_file.GetName()+"[")
#for i in h:
#    i.Draw()
#    canvas.SaveAs(pdf_file.GetName()+"")


for i in range(0,len(h)):
    h[i].Draw()
    canvas.SaveAs(pdf_file.GetName()+"")


canvas.SaveAs(pdf_file.GetName()+"]")

pdf_file.Close()