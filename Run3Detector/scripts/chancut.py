import awkward as ak
import uproot
import ROOT
from array import array

def side(files, side, throughgoing=False):
    var = ["area", "height", "chan", "layer", "event_time_fromTDC"]  # Variables of interest
    sidehist = ROOT.TH1D(f"{side}hist_thru={throughgoing}", f"{side}hist_thru={throughgoing}", 20, 0-0.5, 20-0.5)
    
    for data in uproot.iterate(files, var, library='ak', step_size=1000):  # Looping over the data
        
        # Define the initial mask based on the side parameter
        if side == "t":  # Selects top
            mask = ak.any(
                (data["chan"] == 0) |
                (data["chan"] == 1) |
                (data["chan"] == 4) |
                (data["chan"] == 5) |
                (data["chan"] == 8) |
                (data["chan"] == 9) |
                (data["chan"] == 12) |
                (data["chan"] == 13),
                axis=-1
            )

        elif side == "b":  # Selects bottom
            mask = ak.any(
                (data["chan"] == 2) |
                (data["chan"] == 3) |
                (data["chan"] == 6) |
                (data["chan"] == 7) |
                (data["chan"] == 10) |
                (data["chan"] == 11) |
                (data["chan"] == 14) |
                (data["chan"] == 15),
                axis=-1
            )
        
        elif side == "r":
            mask = ak.any(
                (data["chan"] == 13) |
                (data["chan"] == 15),
                axis=-1
            )
        
        elif side == "l":
            mask = ak.any(
                (data["chan"] == 0) |
                (data["chan"] == 2),
                axis=-1
            )
        
        else:
            raise ValueError("Please enter a valid side")
        
        # Apply the throughgoing mask if throughgoing is True
        if throughgoing:
            throughgoing_mask = (
                ak.any((data["layer"] == 0), axis=-1) & 
                ak.any((data["layer"] == 1), axis=-1) & 
                ak.any((data["layer"] == 2), axis=-1) & 
                ak.any((data["layer"] == 3), axis=-1)
            )
            mask = mask & throughgoing_mask
        
        # Filter data using mask
        sidedata = data[mask]
        flat = ak.flatten(sidedata['chan'], axis=None)
        events = len(flat)
        weight = ak.ones_like(flat)
        if events > 0:
            sidehist.FillN(events, array('d', flat), array('d', weight), 1)
    
    return sidehist

files = ["MilliQan_run729_incomplete.root"]

if __name__ == "__main__":
    outfile = ROOT.TFile("side_histograms.root", "RECREATE")
    
    for side_param in ["t", "b", "l", "r"]:
        hist = side(files=files, side=side_param)
        hist.Write()
        hist_throughgoing = side(files=files, side=side_param, throughgoing=True)
        hist_throughgoing.Write()

    outfile.Close()
