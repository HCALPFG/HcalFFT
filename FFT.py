import sys
import ROOT as r
import numpy as n
import array as a

# r.gStyle.SetOptStat(0)

#------------------------------------------------------------------------------------------------
# User declares only these values
#------------------------------------------------------------------------------------------------

in_file = r.TFile("HFP_QuestionChannels.root");
graph = in_file.Get("HFP13_ETA38_PHI25_T10_SRCTUBE_Ieta38_Iphi25_Depth2 Run 221509reelPosition");
out_file = r.TFile ("fft_out_file.root","RECREATE")

#------------------------------------------------------------------------------------------------
# Make the graph into a histogram
#------------------------------------------------------------------------------------------------

d_x_y = {}

n_points = graph.GetN()
x_array = graph.GetX()
y_array = graph.GetY()

# Declare the histogram
original_histogram = r.TH1F("hist", "hist", 1001, 5799.5, 6800.5)

# Define x axis limits
x_min = original_histogram.GetXaxis().GetXmin()
x_max = original_histogram.GetXaxis().GetXmax()

# Get the x/y values in the graph

x_values = [] 
for i in range(1,n_points):
    x = x_array[i]
    y = y_array[i]
    if x not in d_x_y.keys():
        d_x_y[x] = []
        x_values.append ( x ) 
    d_x_y[x].append ( y )

# Fill the x/y values into the histogram

for x_value in x_values:
    if x_value < x_min : continue
    if x_value > x_max : continue
    y_values = d_x_y[x_value]
    y_value = 0.
    for tmp_y_value in y_values:
        y_value = y_value + tmp_y_value
    y_value = y_value / float(len(y_values))
    bin = original_histogram.FindBin(x_value)
    original_histogram.SetBinContent(bin,y_value)

# Store some useful info about the histogram

original_histogram_nbins = original_histogram.GetNbinsX()
original_function_xmin  = original_histogram.GetXaxis().GetXmin()
original_function_xmax  = original_histogram.GetXaxis().GetXmax()

#------------------------------------------------------------------------------------------------
# Declare canvases
#------------------------------------------------------------------------------------------------

pdf_name = "fft.pdf"

original_canvas       = r.TCanvas()
fft_mag_canvas        = r.TCanvas()
fft_mag_zeroFirstBin_canvas = r.TCanvas()
fft_phase_canvas      = r.TCanvas()
reconstruction_canvas = r.TCanvas()
fft_mag_zeroFirstBin_upperHalf_canvas = r.TCanvas()

#------------------------------------------------------------------------------------------------
# Draw original histogram
#------------------------------------------------------------------------------------------------

original_canvas.cd()

original_histogram.GetXaxis().SetTitle("Reel [mm]")
original_histogram.GetYaxis().SetTitle("Counts")
original_histogram.GetYaxis().SetRangeUser(-0.005, 0.015)

# Draw histogram
original_histogram.Draw()

# Save canvas
original_canvas.SaveAs("png/original_histogram.png");
original_canvas.Print(pdf_name + "[")
original_canvas.Print(pdf_name)

#------------------------------------------------------------------------------------------------
# Now do the magnitude transform
#------------------------------------------------------------------------------------------------

fft_mag_canvas.cd()

# Initialize the Virtual FFT
r.TVirtualFFT.SetTransform(0);

# Declare the FFT histogram
fft_mag_histogram_raw = r.NULL

# Do FFT
fft_mag_histogram_raw = original_histogram.FFT(fft_mag_histogram_raw,"MAG")

# Rescale x-axis on magnitude histogram

fft_mag_histogram_xmin  = 0
fft_mag_histogram_xmax  = fft_mag_histogram_raw.GetXaxis().GetXmax() / original_histogram.GetXaxis().GetXmax();

fft_mag_histogram = r.TH1D("fft_magnitude","fft_magnitude", 
                           fft_mag_histogram_raw.GetNbinsX(), 
                           fft_mag_histogram_xmin, 
                           fft_mag_histogram_xmax)

for bin in range(1, fft_mag_histogram.GetNbinsX() + 1):
    fft_mag_histogram.SetBinContent(bin, fft_mag_histogram_raw.GetBinContent(bin))

# fft_mag_histogram.GetXaxis().SetRangeUser(0., fft_mag_histogram_xmax / 2.0)
fft_mag_histogram.Scale(2.0 / float(original_histogram_nbins))

# Label axes
fft_mag_histogram.GetXaxis().SetTitle("Frequency [1/mm]");
fft_mag_histogram.GetYaxis().SetTitle("FFT magnitude");
fft_mag_histogram.GetYaxis().SetTitleOffset(1.3)

# Draw histogram
fft_mag_histogram.Draw()

# Save canvas
fft_mag_canvas.SaveAs("png/FFT_magnitude.png");
fft_mag_canvas.Print(pdf_name);

#------------------------------------------------------------------------------------------------
# Zero the first bin
#------------------------------------------------------------------------------------------------

fft_mag_zeroFirstBin_canvas.cd()

fft_mag_histogram_zeroFirstBin = fft_mag_histogram.Clone()
fft_mag_histogram_zeroFirstBin.SetBinContent(1,0)

# Draw histogram
fft_mag_histogram_zeroFirstBin.Draw()

# Save canvas
fft_mag_zeroFirstBin_canvas.SaveAs("png/FFT_magnitude_zeroFirstBin.png");
fft_mag_zeroFirstBin_canvas.Print(pdf_name);


#------------------------------------------------------------------------------------------------
# Look at only half of the plot
#------------------------------------------------------------------------------------------------

fft_mag_zeroFirstBin_upperHalf_canvas.cd()

fft_mag_histogram_zeroFirstBin_upperHalf = fft_mag_histogram_zeroFirstBin.Clone()
fft_mag_histogram_zeroFirstBin_upperHalf_xmax = fft_mag_histogram_zeroFirstBin.GetXaxis().GetXmax()
fft_mag_histogram_zeroFirstBin_upperHalf_xmin = fft_mag_histogram_zeroFirstBin_upperHalf_xmax / 2.0
fft_mag_histogram_zeroFirstBin_upperHalf.GetXaxis().SetRangeUser(fft_mag_histogram_zeroFirstBin_upperHalf_xmin,
                                                                 fft_mag_histogram_zeroFirstBin_upperHalf_xmax)

# Draw histogram
fft_mag_histogram_zeroFirstBin_upperHalf.Draw()

# Save canvas
fft_mag_zeroFirstBin_upperHalf_canvas.SaveAs("png/FFT_magnitude_zeroFirstBin_upperHalf.png");
fft_mag_zeroFirstBin_upperHalf_canvas.Print(pdf_name);

out_file.cd()
fft_mag_histogram_zeroFirstBin_upperHalf.Write("FFT_magnitude")

#------------------------------------------------------------------------------------------------
# Now do the phase transform
#------------------------------------------------------------------------------------------------

fft_phase_canvas.cd()

# What does this do?
r.TVirtualFFT.SetTransform(0);

# Declare the FFT histogram
fft_phase_histogram_raw = r.NULL

# Do FFT
fft_phase_histogram_raw = original_histogram.FFT(fft_phase_histogram_raw,"PH")

# Rescale x-axis on phase histogram

fft_phase_histogram_xmin  = 0
fft_phase_histogram_xmax  = fft_phase_histogram_raw.GetXaxis().GetXmax() / original_histogram.GetXaxis().GetXmax();

fft_phase_histogram = r.TH1D("fft_phase","fft_phase", 
                             fft_phase_histogram_raw.GetNbinsX(), 
                             fft_phase_histogram_xmin, 
                             fft_phase_histogram_xmax)

for bin in range(1, fft_phase_histogram.GetNbinsX() + 1):
    fft_phase_histogram.SetBinContent(bin, fft_phase_histogram_raw.GetBinContent(bin))

fft_phase_histogram.Scale(1.0 / r.TMath.Sqrt(float(original_histogram_nbins)))

# Label axes
fft_phase_histogram.GetXaxis().SetTitle("Frequency [1/mm]");
fft_phase_histogram.GetYaxis().SetTitle("FFT phase");

# Draw histogram
fft_phase_histogram.Draw()

# Save canvas
fft_phase_canvas.SaveAs("png/FFT_phase.png");
fft_phase_canvas.Print(pdf_name);

#------------------------------------------------------------------------------------------------
# Revert back to original (inverse FFT)
#------------------------------------------------------------------------------------------------

# Go to the reconstruction canvas

reconstruction_canvas.cd()

# First get the original FFT

fft = r.TVirtualFFT.GetCurrentTransform()

# Get full real * i*imaginary contents of the original FFT

fft_real_array = n.array ([])
fft_imag_array = n.array ([])

fft_real_array.resize(original_histogram_nbins)
fft_imag_array.resize(original_histogram_nbins)

fft.GetPointsComplex(fft_real_array, fft_imag_array)

# Set up the inverse FFT and run it

fft_inverse = r.TVirtualFFT.FFT(1, a.array("i",[original_histogram_nbins]) , "C2R M K")
fft_inverse.SetPointsComplex(fft_real_array, fft_imag_array)
fft_inverse.Transform();

# Reconstruct the original histogram

original_histogram_reconstruction_raw = r.NULL
original_histogram_reconstruction_raw = r.TH1D.TransformHisto(fft_inverse, original_histogram_reconstruction_raw, "Re")
original_histogram_reconstruction = r.TH1D("original_histogram_reconstruction","original_histogram_reconstruction",
                                           original_histogram_nbins, original_function_xmin, original_function_xmax)

for bin in range (1, original_histogram_nbins + 1):
    original_histogram_reconstruction.SetBinContent(bin, original_histogram_reconstruction_raw.GetBinContent(bin))
original_histogram_reconstruction.Scale (1.0 / float(original_histogram_nbins))

original_histogram_reconstruction.GetXaxis().SetTitle("Reel [mm]")
original_histogram_reconstruction.GetYaxis().SetTitle("Counts")

original_histogram_reconstruction.Draw();

reconstruction_canvas.SaveAs("png/reconstruction.png")
reconstruction_canvas.Print(pdf_name)

#------------------------------------------------------------------------------------------------
# Close pdf canvas
#------------------------------------------------------------------------------------------------

original_canvas.Print(pdf_name + "]")
