#For generating a comparison rates plot for through-going stuff in the FORMOSA detector, not as a stack
import ROOT

# Open ROOT file
root_file = ROOT.TFile.Open("/eos/experiment/formosa/milliqanOffline/Run3Detector/scripts/outputs/exploringWindowTiming_run737_highHV_caseComparison.root")

# Get histograms from the file
hist1 = root_file.Get("rates_throughGoing_withCut")
hist2 = root_file.Get("rates_withCut")
hist3 = root_file.Get("rates_withCut_withPanels")
hist4 = root_file.Get("rates_withCut_withPanels_threshold")
hist5 = root_file.Get("rates_leftover")
hist6 = root_file.Get("rates_withCut_onePanel")
hist7 = root_file.Get("rates_withCut_onePanel_threshold")

# Create a canvas to draw on
canvas = ROOT.TCanvas("canvas", "Histogram Plot", 800, 600)

# Create THStack to stack histograms
stack = ROOT.THStack("stack", "Stacked Histograms")

# Draw histograms
hist1.SetLineColor(1)
hist2.SetFillColor(4)
hist2.SetLineColor(4)
hist3.SetFillColor(2)
hist3.SetLineColor(2)
hist4.SetFillColor(6)
hist4.SetLineColor(6)
hist5.SetFillColor(7)
hist5.SetLineColor(7)
hist6.SetFillColor(8)
hist6.SetLineColor(8)
hist7.SetFillColor(ROOT.kOrange)
hist7.SetLineColor(ROOT.kOrange)

stack.Add(hist3)
stack.Add(hist4)
stack.Add(hist2)
stack.Add(hist5)
stack.Add(hist6)
stack.Add(hist7)

#Draw Stack
hist1.SetTitle("Comparing Through-Going Populations, Run709, Loose Timing Cut")
hist1.Draw("same")
hist2.Draw("same")
hist3.Draw("same")
hist4.Draw("same")
hist6.Draw("same")
hist7.Draw("same")
hist5.Draw("same")

# Create legend
legend = ROOT.TLegend(0.5, 0.6, 0.9, 0.9)  # Adjust coordinates as needed
legend.AddEntry(hist1, "All through-going activity", "l")
legend.AddEntry(hist2, "Through-going activity, no panels", "l")
legend.AddEntry(hist3, "Through-going activity, with panels, below panel threshold", "l")
legend.AddEntry(hist4, "Through-going activity, with panels, above panel threshold", "l")
legend.AddEntry(hist6, "Through-going, with one panels", "l")
legend.AddEntry(hist7, "Through-going, with two panels (only one above threshold)", "l")
legend.AddEntry(hist5, "Leftover Through-Going Events", "l")
legend.Draw()

ROOT.gStyle.SetOptStat(0)

# Update canvas
canvas.Update()

# Keep the program running until the user closes the canvas
canvas.SaveAs("test.root")
