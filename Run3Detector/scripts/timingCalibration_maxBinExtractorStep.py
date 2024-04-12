#This is STEP TWO in the timing calibration gen process.
#It extracts the max bins from the calibrated files output by STEP ONE
import ROOT

# Open the ROOT file
root_file = ROOT.TFile.Open("outputs/throughgoingBarDifferences.root")

# List to store bin centers
bin_centers = {}

# Specify the string pattern
histo_name_pattern = "barDifferencesCH"

# Iterate over all histograms in the file
for key in root_file.GetListOfKeys():
    obj = key.ReadObj()
    # Check if the object is a histogram and its name contains the specified pattern
    if isinstance(obj, ROOT.TH1) and histo_name_pattern in obj.GetName():
        if obj.GetName().replace(histo_name_pattern,"") == str(16) or obj.GetName().replace(histo_name_pattern,"") == str(17):
            continue
        # Get the bin with the maximum content
        max_bin = obj.GetMaximumBin()
        # Get the bin center of the maximum bin
        bin_center = obj.GetBinCenter(max_bin)
        bin_centers[int(obj.GetName().replace(histo_name_pattern,""))] = bin_center

# Print the list of bin centers
print(bin_centers)

# Close the ROOT file
root_file.Close()
