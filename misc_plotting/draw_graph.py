import ROOT as r
in_file = r.TFile("HFP_QuestionChannels.root");
graph = in_file.Get("HFP13_ETA38_PHI25_T10_SRCTUBE_Ieta38_Iphi25_Depth2 Run 221509reelPosition");

canvas = r.TCanvas()
canvas.cd()

graph.Draw("AP")
graph.GetHistogram().GetYaxis().SetTitleOffset(1.5)
graph.GetHistogram().GetYaxis().SetRangeUser(-0.005, 0.015)
graph.GetHistogram().GetYaxis().SetTitle("Histogram mean [Linear ADC]")
graph.GetHistogram().GetXaxis().SetTitle("Reel [mm]")
graph.GetHistogram().GetXaxis().SetRangeUser(5800, 6800)
graph.Draw("AP")

canvas.SaveAs("sourcing_zoomed_plot.png")
