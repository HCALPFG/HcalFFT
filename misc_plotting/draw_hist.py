import ROOT as r
in_file = r.TFile("hists.root")
graph = in_file.Get("hist")

canvas = r.TCanvas()
canvas.cd()

graph.Draw()
graph.GetYaxis().SetTitleOffset(1.5)
graph.GetYaxis().SetTitle("Histogram mean [Linear ADC]")
graph.GetXaxis().SetTitle("Reel [mm]")

canvas.SaveAs("sourcing_zoomed_hist.png")
