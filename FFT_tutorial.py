import ROOT as r
import numpy as n
import array as a
import sys

r.gStyle.SetOptStat(0)

#------------------------------------------------------------------------------------------------
# Define histograms
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

frequencies = [0.5, 1., 2.]
phases      = [ 0., 0., 0.] 
amps        = [ 1., 2., 0.5]

original_histogram_nbins = 200;
original_function_xmin  = 0.;
original_function_xmax  = 10;

pdf_name = "fft_tutorial.pdf"

#------------------------------------------------------------------------------------------------
# Draw original function and histogram
#------------------------------------------------------------------------------------------------

# Get the size of the histogram
original_function_width = original_function_xmax - original_function_xmin 

# Declare function

function_string = ""
for i_frequency, frequency in enumerate(frequencies):
    function_string = function_string + str(amps[i_frequency])+"*sin(2.0*TMath::Pi()*x*"+str(frequency)+")+"
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

hist = original_histogram

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

middle_bin = n_FFT / 2;

spacing = getFrequencySpacing(sample_frequency, n_FFT)

min_true_frequency = (getRealFrequency(middle_bin, sample_frequency, n_FFT) + (0.5 * getFrequencySpacing(sample_frequency, n_FFT))) * -1.0;
max_true_frequency = min_true_frequency * -1.0;

real_output_hist  = r.TH1F("real_output" ,"real_output" , n_FFT + 1, min_true_frequency, max_true_frequency);
imag_output_hist  = r.TH1F("imag_output" ,"imag_output" , n_FFT + 1, min_true_frequency, max_true_frequency);
mag_output_hist   = r.TH1F("mag_output"  ,"mag_output"  , n_FFT + 1, min_true_frequency, max_true_frequency);
phase_output_hist = r.TH1F("phase_output","phase_output", n_FFT + 1, min_true_frequency, max_true_frequency);

#------------------------------------------------------------------------------------------------
# Fill the FFT histograms
#------------------------------------------------------------------------------------------------

real = r.Double()
imag = r.Double()

used_bins = []

for i in range (0, n_FFT):
    
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

real_output_hist .GetXaxis().SetTitle("Frequency [Hz]");
imag_output_hist .GetXaxis().SetTitle("Frequency [Hz]");
mag_output_hist  .GetXaxis().SetTitle("Frequency [Hz]");
phase_output_hist.GetXaxis().SetTitle("Frequency [Hz]");

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
canvas.SaveAs("png/tutorial_original_histogram.png")
canvas.Print(pdf_name)

real_output_hist.Draw()
canvas.SaveAs("png/tutorial_FFT_real.png")
canvas.Print(pdf_name)

imag_output_hist.Draw()
canvas.SaveAs("png/tutorial_FFT_imag.png")
canvas.Print(pdf_name)

mag_output_hist.Draw()
canvas.SaveAs("png/tutorial_FFT_magnitude.png")
canvas.Print(pdf_name)

phase_output_hist.Draw()
canvas.SaveAs("png/tutorial_FFT_phase.png")
canvas.Print(pdf_name)

canvas.Print(pdf_name + "]")

#------------------------------------------------------------------------------------------------
# Done!
#------------------------------------------------------------------------------------------------
