#For looking at the muons for comparison to the weird background that shows up in the begining of our timing windows
import uproot
import awkward
import sys
import os
from array import array

original_directory = os.getcwd()
sys.path.append('/eos/experiment/formosa/steenis-general-functions')
from general_functions import *
sys.path.append(original_directory)
#------------------------------------------------------------------------------------------------#
filelist = ["/eos/experiment/formosa/commissioning/data/hadded_outputs/MilliQan_Run709_v35_throughGoingPanelsSkim.root:t"]
variables = ["timeFit", "area", "chan", "layer", "timeFit_module_calibrated", "height"]

outfile = ROOT.TFile("outputs/exploringWindowTiming_muons.root", "RECREATE")
timeFitDiff = ROOT.TH1F("timeFitDifferences", "timeFitDifferences", 200,-100,100)
shapeHist = ROOT.TH2F("shapeDist", "shapeDist", 200, 0, 1300, 100, 0, 200e3)
wfTimeHist = ROOT.TH2F("wfTime", "wfTime", 200, 0, 1300, 100, 0, 200e3)
wfTime_layer = ROOT.TH2F("wfTime_layer", "wfTime_layer", 200, 0, 1300, 4, 0, 4)

panelThresholds = {16:7000, 18:4000}
barThreshold = 80000

counter = 0
for data in uproot.iterate(filelist, variables, library='ak', step_size=100000): 
    if counter%100==0:
        print(counter, " out of ", len(filelist), " runs")
    counter+=1

    print("NEW STEP")

    barThresholdMask = data['area']>barThreshold
    layer0 = data['layer']==0
    layer1 = data['layer']==1
    layer2 = data['layer']==2
    layer3 = data['layer']==3
    chan16 = data['chan']==16
    notChan16 = data['chan']!=16
    chan18 = data['chan']==18
    notChan18 = data['chan']!=18
    panelThresholdMaskFront = data['area']>panelThresholds[16]
    panelThresholdMaskBack = data['area']>panelThresholds[18]

    throughGoingMask = ak.any(((layer0)), axis=-1) & ak.any(((layer1)), axis=-1) & ak.any(((layer2)), axis=-1) & ak.any(((layer3)), axis=-1)
    locationMask = throughGoingMask & ak.any(chan16, axis=-1) & ak.any(chan18, axis=-1) 
    data = data[locationMask] #We only want to consider throughGoingNonPanel events (event-level cut) 
    
    maxPerLayerMask = ((data['area']==ak.max(data['area'][data['layer']==0], axis=-1))) | ((data['area']==ak.max(data['area'][data['layer']==1], axis=-1))) | ((data['area']==ak.max(data['area'][data['layer']==2], axis=-1))) | ((data['area']==ak.max(data['area'][data['layer']==3], axis=-1)))
    data = data[maxPerLayerMask] #Only the maxes of each layer (pulse-level cut)
    
    onePerLayerCut = ak.num(data['area'], axis=-1)==4
    data = data[onePerLayerCut] #Gets rid of the cases where a max-pulse is reapeated in area
    
    tempDiffs = data['timeFit_module_calibrated'][data['layer']==0] - data['timeFit_module_calibrated'][data['layer']==3]
    tempDiffs = ak.pad_none(tempDiffs, 1, axis=-1)
    tempDiffs = ak.flatten(tempDiffs, axis=None)
    tempDiffs = ak.broadcast_arrays(tempDiffs, data['area'])[0]
    tempDiffs = ak.drop_none(tempDiffs)
    data['timeDiffs'] = tempDiffs
    '''timeDiffMask = (ak.any(data['timeDiffs'] < -40, axis=-1)) & (ak.any(data['timeDiffs'] > -65, axis=-1))
    data = data[timeDiffMask] #Removing events that don't have a layer0-layer3 timing difference in the negative peak.'''
    
    '''muonSeperatorMask = data['area'] > (20.31*data['height'])+4593.75
    nonSaturationMask = data['height'] < 1230
    nonMuonMask = muonSeperatorMask & nonSaturationMask
    data = data[nonMuonMask] #Restrict to stuff outside the muon band (pulse-level cut)'''
 
    timeDiffs03 = ak.flatten(data['timeDiffs'], axis=None)
    areas = ak.flatten(data['area'], axis=None)
    heights = ak.flatten(data['height'], axis=None)
    timeFit_module_calibrated = ak.flatten(data['timeFit_module_calibrated'], axis=None)
    layers = ak.flatten(data['layer'], axis=None)

    if len(timeDiffs03) != len(areas) or len(timeDiffs03) != len(heights) or len(timeDiffs03) != len(timeFit_module_calibrated):
        print("RED ALERT! Different Lengths of Data Arrays!")
    if len(timeDiffs03)==0:
        continue
    
    n_events = len(timeDiffs03)
    weights = ak.ones_like(timeDiffs03)
    timeFitDiff.FillN(n_events, array('d', timeDiffs03), array('d', weights), 1)
    shapeHist.FillN(n_events, array('d', heights), array('d', areas), array('d', weights), 1)
    wfTimeHist.FillN(n_events, array('d', timeFit_module_calibrated), array('d', areas), array('d', weights), 1)
    wfTime_layer.FillN(n_events, array('d', timeFit_module_calibrated), array('d', layers), array('d', weights), 1)
            
timeFitDiff.SetTitle("Cut Check, Layer0-Layer3")
timeFitDiff.GetXaxis().SetTitle("Time Difference [ns]")
timeFitDiff.Write()

shapeHist.SetTitle("Cut Check, Pulse Shape")
shapeHist.GetXaxis().SetTitle("Height [mV]")
shapeHist.GetYaxis().SetTitle("Area [pVs]")
shapeHist.Write()

wfTimeHist.SetTitle("Timing within Readout Window, Through-Going w Panels")
wfTimeHist.GetXaxis().SetTitle("timeFit_module_calibrated [ns]")
wfTimeHist.GetYaxis().SetTitle("Area [pVs]")
wfTimeHist.Write()

wfTime_layer.SetTitle("Timing within Readout Window, Through-Going w Panels")
wfTime_layer.GetXaxis().SetTitle("timeFit_module_calibrated [ns]")
wfTime_layer.GetYaxis().SetTitle("Layer")
wfTime_layer.Write()

outfile.Close()
