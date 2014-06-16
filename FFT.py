import ROOT as r
import numpy as n
import array as a
import sys

r.gStyle.SetOptStat(0)

#------------------------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------------------------

def getRealFrequency(bin_number, sample_frequency, nFFT):
    if bin_number <= (n_FFT / 2): 
        real_frequency = float(bin_number) * sample_frequency / float(nFFT);
    else:
        real_frequency = -1.0 * float(nFFT - bin_number) * sample_frequency / float(nFFT);
    return real_frequency;

def getFrequencySpacing(sample_frequency, nFFT):
    spacing = sample_frequency / float(nFFT)
    return spacing

#------------------------------------------------------------------------------------------------
# User declares only these values
#------------------------------------------------------------------------------------------------

in_file = r.TFile("HFP_QuestionChannels.root");
graph = in_file.Get("HFP13_ETA38_PHI25_T10_SRCTUBE_Ieta38_Iphi25_Depth2 Run 221509reelPosition");
pdf_name = "fft.pdf"

#------------------------------------------------------------------------------------------------
# Make the graph into a histogram
#------------------------------------------------------------------------------------------------

n_points = graph.GetN()
x_array = graph.GetX()
y_array = graph.GetY()

# Declare the histogram

hist = r.TH1F("original_hist", "original_hist", 1001, 5799.5, 6800.5)
hist.GetXaxis().SetTitle("Reel [mm]");
hist.GetYaxis().SetTitle("Histogram mean [Linear ADC]");

# Define x axis limits
x_min = hist.GetXaxis().GetXmin()
x_max = hist.GetXaxis().GetXmax()

# Get the x/y values in the graph

d_x_y = {}
x_values = [] 
for i in range(1,n_points):
    x = x_array[i]
    y = y_array[i]
    if x not in d_x_y.keys():
        d_x_y[x] = []
        x_values.append ( x ) 
    d_x_y[x].append ( y )

#------------------------------------------------------------------------------------------------
# Fill the x/y values into the histogram
#------------------------------------------------------------------------------------------------

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

#------------------------------------------------------------------------------------------------
# Store some useful info about the histogram
#------------------------------------------------------------------------------------------------

n_bins           = hist.GetNbinsX()
x_min            = hist.GetXaxis().GetXmin()
x_max            = hist.GetXaxis().GetXmax()
n_FFT            = 2 * ((n_bins + 1) / 2 + 1);
n_dimension      = 1;
sample_period    = (x_max - x_min) / float (n_bins);
sample_frequency = 1.0 / sample_period;

#------------------------------------------------------------------------------------------------
# Store the input points
#------------------------------------------------------------------------------------------------

input_points = n.array([])
input_points.resize(n_FFT);
for bin in range (1, n_bins + 1):
    input_points[bin - 1] = hist.GetBinContent(bin)

#------------------------------------------------------------------------------------------------
# Do the FFT!
#------------------------------------------------------------------------------------------------

fft = r.TVirtualFFT.FFT(n_dimension, a.array("i",[n_FFT]), "R2C");
fft.SetPoints(input_points);
fft.Transform();

#------------------------------------------------------------------------------------------------
# Store the information about the DC:
#------------------------------------------------------------------------------------------------

dc_real = r.Double()
dc_imag = r.Double()
fft.GetPointComplex(0, dc_real, dc_imag);

#------------------------------------------------------------------------------------------------
# Make the histograms corresponding to the FFT output
#------------------------------------------------------------------------------------------------

mid_bin = n_FFT / 2;

my_map = list(range(1,mid_bin))
my_map_reverse = list(my_map)
my_map_reverse.reverse()
my_map_reverse = [i * -1 for i in my_map_reverse]
my_map = my_map_reverse + my_map
my_map = [i + mid_bin for i in my_map]

middle_bin = n_FFT / 2;

spacing = getFrequencySpacing(sample_frequency, n_FFT)

min_true_frequency = (getRealFrequency(middle_bin, sample_frequency, n_FFT) + (0.5 * getFrequencySpacing(sample_frequency, n_FFT))) * -1.0;
max_true_frequency = min_true_frequency * -1.0;

real_output_hist  = r.TH1F("real_output" ,"real_output" , n_FFT + 1, min_true_frequency, max_true_frequency);
imag_output_hist  = r.TH1F("imag_output" ,"imag_output" , n_FFT + 1, min_true_frequency, max_true_frequency);
mag_output_hist   = r.TH1F("mag_output"  ,"mag_output"  , n_FFT + 1, min_true_frequency, max_true_frequency);
phase_output_hist = r.TH1F("phase_output","phase_output", n_FFT + 1, min_true_frequency, max_true_frequency);

print "bin width =", real_output_hist.GetBinWidth(1)
print "plot range =", min_true_frequency, max_true_frequency
print "middle bin =", n_FFT / 2

#------------------------------------------------------------------------------------------------
# Fill the FFT histograms
#------------------------------------------------------------------------------------------------

real = r.Double()
imag = r.Double()

used_bins = []

for i in range (1, n_FFT):
    
    fft.GetPointComplex(i, real, imag);
    mag   = r.TMath.Sqrt((real * real) + (imag * imag));
    phase = r.TMath.ATan(imag / real);
    true_frequency = getRealFrequency(i, sample_frequency, n_FFT);
    bin = real_output_hist.FindBin(true_frequency);

    bin_low_edge  = real_output_hist.GetBinLowEdge(bin)
    bin_high_edge = real_output_hist.GetBinLowEdge(bin+1)
    up_spacing   = bin_high_edge - true_frequency
    down_spacing = true_frequency - bin_low_edge

    if i == middle_bin:
        print "middle bin:", real, imag

    if i == 0:
        print "DC:", real, imag
    
    real_output_hist .SetBinContent ( bin, real  );
    imag_output_hist .SetBinContent ( bin, imag  );
    mag_output_hist  .SetBinContent ( bin, mag   );
    phase_output_hist.SetBinContent ( bin, phase );

    if bin not in used_bins:
        used_bins.append ( bin ) 
    else:
        print "Used this bin already:", bin

#------------------------------------------------------------------------------------------------
# Style the FFT histograms
#------------------------------------------------------------------------------------------------

real_output_hist .Scale ( 1.0 / float(n_bins) );
imag_output_hist .Scale ( 1.0 / float(n_bins) );
mag_output_hist  .Scale ( 1.0 / float(n_bins) );
phase_output_hist.Scale ( 1.0 / float(n_bins) );

real_output_hist .GetXaxis().SetTitle("Reel frequency [1/mm]");
imag_output_hist .GetXaxis().SetTitle("Reel frequency [1/mm]");
mag_output_hist  .GetXaxis().SetTitle("Reel frequency [1/mm]");
phase_output_hist.GetXaxis().SetTitle("Reel frequency [1/mm]");

real_output_hist .GetYaxis().SetTitle("Real component [A.U.]");
imag_output_hist .GetYaxis().SetTitle("Imaginary component [A.U.]");
mag_output_hist  .GetYaxis().SetTitle("Magnitude [A.U.]");
phase_output_hist.GetYaxis().SetTitle("Phase [A.U.]");

hist             .GetYaxis().SetTitleOffset(1.4)
real_output_hist .GetYaxis().SetTitleOffset(1.4)
imag_output_hist .GetYaxis().SetTitleOffset(1.4)
mag_output_hist  .GetYaxis().SetTitleOffset(1.4)
phase_output_hist.GetYaxis().SetTitleOffset(1.4)

#------------------------------------------------------------------------------------------------
# Draw the FFT histograms
#------------------------------------------------------------------------------------------------

canvas = r.TCanvas()
canvas.cd()
canvas.SetTopMargin(0.1)

canvas.Print(pdf_name + "[")

hist.Draw()
canvas.SaveAs("png/original_histogram.png")
canvas.Print(pdf_name)

real_output_hist.Draw()
canvas.SaveAs("png/FFT_real.png")
canvas.Print(pdf_name)

imag_output_hist.Draw()
canvas.SaveAs("png/FFT_imag.png")
canvas.Print(pdf_name)

mag_output_hist.Draw()
canvas.SaveAs("png/FFT_magnitude.png")
canvas.Print(pdf_name)

phase_output_hist.Draw()
canvas.SaveAs("png/FFT_phase.png")
canvas.Print(pdf_name)

canvas.Print(pdf_name + "]")

#------------------------------------------------------------------------------------------------
# Done!
#------------------------------------------------------------------------------------------------
