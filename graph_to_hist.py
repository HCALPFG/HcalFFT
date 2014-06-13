import ROOT as r

in_file = r.TFile("HFP_QuestionChannels.root");
out_file = r.TFile("hists.root", "RECREATE");
graph = in_file.Get("HFP13_ETA38_PHI25_T10_SRCTUBE_Ieta38_Iphi25_Depth2 Run 221509reelPosition");

d_x_y = {}

n = graph.GetN()
x_array = graph.GetX()
y_array = graph.GetY()

hist = r.TH1F("hist", "hist", 1001, 5799.5, 6800.5)
x_min = hist.GetXaxis().GetXmin()
x_max = hist.GetXaxis().GetXmax()
y_min = hist.GetBinContent(hist.GetMinimumBin())
y_max = hist.GetBinContent(hist.GetMaximumBin())

x_values = [] 

for i in range(1,n):
    x = x_array[i]
    y = y_array[i]
    if x not in d_x_y.keys():
        d_x_y[x] = []
        x_values.append ( x ) 
    d_x_y[x].append ( y )

for x_value in x_values:
    if x_value < x_min : continue
    if x_value > x_max : continue
    y_values = d_x_y[x_value]
    y_value = 0.
    for tmp_y_value in y_values:
        y_value = y_value + tmp_y_value
    y_value = y_value / float(len(y_values))
    bin = hist.FindBin(x_value)
    hist.SetBinContent(bin,y_value)
    
out_file.cd()
hist.Write()
    
