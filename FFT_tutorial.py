import ROOT as r
import numpy as n
import array as a

r.gStyle.SetOptStat(0)

#------------------------------------------------------------------------------------------------
# User declares only these values
#------------------------------------------------------------------------------------------------

# frequencies = [0.5, 1., 2.]
# phases      = [ 0., 0., 0.] 
# amps        = [ 1., 2., 0.5]


frequencies = [ 0.5 ]
phases      = [ 0.0 ] 
amps        = [ 1.0 ]

original_histogram_nbins = 200
original_function_xmin  = 0.;
original_function_xmax  = 10;

#------------------------------------------------------------------------------------------------
# Declare canvases
#------------------------------------------------------------------------------------------------

pdf_name = "fft_example.pdf"

original_canvas       = r.TCanvas()
fft_mag_canvas        = r.TCanvas()
fft_phase_canvas      = r.TCanvas()
reconstruction_canvas = r.TCanvas()

#------------------------------------------------------------------------------------------------
# Draw original function and histogram
#------------------------------------------------------------------------------------------------

original_canvas.cd()

# Get the size of the histogram
minimum_frequency = min(frequencies)
maximum_frequency = max(frequencies)
original_function_width = original_function_xmax - original_function_xmin 

# Declare function

function_string = ""
for i_frequency, frequency in enumerate(frequencies):
    function_string = function_string + str(amps[i_frequency])+"*sin(2.0*TMath::Pi()*(x*"+str(frequency)+"+"+str(phases[i_frequency])+"))+"
function_string = function_string[:-1]
print function_string

original_function = r.TF1("original_function", function_string, original_function_xmin, original_function_xmax)
# original_function = r.TF1("original_function", "sin(2.*1.5*TMath::Pi()*x)", original_function_xmin, original_function_xmax)

# Declare histogram
original_histogram = r.TH1D("original_histogram", "original_histogram", original_histogram_nbins, original_function_xmin, original_function_xmax)
original_histogram.GetXaxis().SetTitle("Time [s]");
original_histogram.GetYaxis().SetTitle("A.U.");

# Fill histogram
for ibin in range(0, original_histogram_nbins):
    bin_low_edge = original_function_xmin + float(ibin) * (original_function_width / float(original_histogram_nbins))
    bin          = original_histogram.FindBin( bin_low_edge )
    bin_center   = original_histogram.GetBinCenter( bin ) 
    y_value      = original_function.Eval( bin_center ) 
    original_histogram.SetBinContent(bin, y_value )
    
# Draw histogram and function
original_histogram.Draw()
# original_function.Draw("SAME")

# Save canvas
original_canvas.SaveAs("png/tutorial_original_histogram.png");
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
fft_mag_histogram.GetXaxis().SetTitle("Frequency [Hz]");
fft_mag_histogram.GetYaxis().SetTitle("FFT magnitude");
# fft_mag_histogram.GetXaxis().SetRangeUser(0, 3.0)

# Draw histogram
fft_mag_histogram.Draw()

# Save canvas
fft_mag_canvas.SaveAs("png/tutorial_FFT_magnitude.png");
fft_mag_canvas.Print(pdf_name);

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
fft_phase_histogram.GetXaxis().SetTitle("Frequency [Hz]");
fft_phase_histogram.GetYaxis().SetTitle("FFT phase");

# Draw histogram
fft_phase_histogram.Draw()

# Save canvas
fft_phase_canvas.SaveAs("png/tutorial_FFT_phase.png");
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

original_histogram_reconstruction.GetXaxis().SetTitle("Time [s]")
original_histogram_reconstruction.GetYaxis().SetTitle("A.U.")

original_histogram_reconstruction.Draw();

reconstruction_canvas.SaveAs("png/tutorial_reconstruction.png")
reconstruction_canvas.Print(pdf_name)

#------------------------------------------------------------------------------------------------
# Close pdf canvas
#------------------------------------------------------------------------------------------------

original_canvas.Print(pdf_name + "]")
