#For comparing Matthew's method of calculating rates with one including a trigger-able time cut for through-going With-Panel activity in the detctor.
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
filelist = ["/eos/experiment/formosa/commissioning/data/hadded_outputs/MilliQan_Run737_incomplete_muonSkim.root:t"]
#filelist = ["/eos/experiment/formosa/commissioning/data/hadded_outputs/MilliQan_Run709_v35_throughGoingPanelsSkim.root:t"]
variables = ["timeFit", "area", "chan", "layer", "timeFit_module_calibrated", "height", "event_time_fromTDC"]

timing_range = [1714.04e6, 1714.11e6] #uproot_range_finder(filelist[0], "event_time_fromTDC")
outfile = ROOT.TFile("outputs/exploringWindowTiming_comparingRates_timeCut_withPanels_run737_incomplete.root", "RECREATE")

times = ROOT.TH1F("times", "time", 500,timing_range[0],timing_range[1])
times_withCut = ROOT.TH1F("times_withCut", "times_withCut", 500,timing_range[0],timing_range[1])
rates = ROOT.TH1F("rates", "rates", 500,timing_range[0],timing_range[1])
rates_withCut = ROOT.TH1F("rates_withCut", "rates_withCut", 500,timing_range[0],timing_range[1])
windowTimes = ROOT.TH2F("windowTimes", "windowTimes", 200,0,1300, 4, 0, 4)
windowTimes_withCut = ROOT.TH2F("windowTimes_withCut", "windowTimes_withCut", 200, 0, 1300, 4, 0, 4)

panelThresholds = {16:7000, 18:4000}
barThreshold = 80000

counter = 0
counter_limit = None
for data in uproot.iterate(filelist, variables, library='ak', step_size=100000): 
    if counter%100==0:
        print(counter, " out of ", len(filelist), " runs")  
          
    if counter_limit is not None and counter==counter_limit: 
        continue
    counter+=1

    print("NEW STEP")
    
    #Broadcast event_time_fromTDC to proper shape
    data['event_time_fromTDC'] = ak.broadcast_arrays(ak.flatten(ak.pad_none(data['event_time_fromTDC'], 1, axis=-1), axis=None), data['area'])[0]

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
    timingCut = (data['timeFit_module_calibrated']>500) & (data['timeFit_module_calibrated']<700)

    throughGoingMask = ak.any(((layer0)), axis=-1) & ak.any(((layer1)), axis=-1) & ak.any(((layer2)), axis=-1) & ak.any(((layer3)), axis=-1)
    throughGoingMask_withCut = ak.any(((layer0) & timingCut), axis=-1) & ak.any(((layer1) & timingCut), axis=-1) & ak.any(((layer2) & timingCut), axis=-1) & ak.any(((layer3) & timingCut), axis=-1) 
    locationMask = throughGoingMask & ak.any(chan16, axis=-1) & ak.any(chan18, axis=-1) 
    locationMask_withCut = throughGoingMask_withCut & ak.any(chan16, axis=-1) & ak.any(chan18, axis=-1) 
    data_copy = data
    data = data[locationMask] #We only want to consider throughGoingNonPanel events (event-level cut) 
    data_withCut = data_copy[locationMask_withCut] #To compare for the case where we look for through-going stuff that is triggerable in time
 
    maxes = data['event_time_fromTDC'][data['area']==ak.max(data['area'], axis=-1)] #Only want one pulse from each event
    event_time_fromTDC_noCut = ak.flatten(maxes, axis=None)
    n_events_noCut = len(event_time_fromTDC_noCut)
    weights_noCut = ak.ones_like(event_time_fromTDC_noCut)
    times.FillN(n_events_noCut, array('d', event_time_fromTDC_noCut), array('d', weights_noCut), 1)
    
    maxes_withCut = data_withCut['event_time_fromTDC'][data_withCut['area']==ak.max(data_withCut['area'], axis=-1)]
    event_time_fromTDC_cut = ak.flatten(maxes_withCut, axis=None)
    n_events_cut = len(event_time_fromTDC_cut)
    weights_cut = ak.ones_like(event_time_fromTDC_cut)
    times_withCut.FillN(n_events_cut, array('d', event_time_fromTDC_cut), array('d', weights_cut), 1)
    
    #For checking the timing cut
    windowTimes_data = ak.flatten(data['timeFit_module_calibrated'], axis=None)
    layers_data = ak.flatten(data['layer'], axis=None)
    windowTimes.FillN(n_events_noCut, array('d', windowTimes_data), array('d', layers_data), array('d', weights_noCut), 1)
    windowTimes_withCut_data = ak.flatten(data_withCut['timeFit_module_calibrated'], axis=None)
    layers_withCut_data = ak.flatten(data_withCut['layer'], axis=None)
    windowTimes_withCut.FillN(n_events_cut, array('d', windowTimes_withCut_data), array('d', layers_withCut_data), array('d', weights_cut), 1)
    
times.SetTitle("event_time_fromTDC for Through-Going With-Panel Activity")
times.GetXaxis().SetTitle("Event Time [ns]")
times.Write()

times_withCut.SetTitle("event_time_fromTDC for Through-Going With-Panel Activity, Triggerable Time Cut")
times_withCut.GetXaxis().SetTitle("Event Time [ns]")
times_withCut.Write()

rates.Add(times)
rates.GetXaxis().SetTimeDisplay(1)
rates.GetXaxis().SetLabelSize(0.02)
rates.GetXaxis().SetLabelOffset(0.01)
rates.GetXaxis().SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 02:00:00")

rates_withCut.Add(times_withCut)
rates_withCut.GetXaxis().SetTimeDisplay(1)
rates_withCut.GetXaxis().SetLabelSize(0.02)
rates_withCut.GetXaxis().SetLabelOffset(0.01)
rates_withCut.GetXaxis().SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 02:00:00")

width = rates.GetBinWidth(30)
inv_width = 1/width
rates.Scale(inv_width)

if rates_withCut.GetBinWidth(30) != width:
    raise ValueError("The two rates bin widths are different!")

else:
    rates_withCut.Scale(inv_width)

rates.SetTitle("Rates for Through-Going With-Panel Activity")
rates.GetXaxis().SetTitle("Event Time [ns]")
rates.GetYaxis().SetTitle("Rate [Hz]")
rates.Write()

rates_withCut.SetTitle("Rates for Through-Going With-Panel Activity, Triggerable Time Cut")
rates_withCut.GetXaxis().SetTitle("Event Time [ns]")
rates_withCut.GetYaxis().SetTitle("Rate [Hz]")
rates_withCut.Write()

windowTimes.SetTitle("Checking Timing Cut, No Cut")
windowTimes.GetXaxis().SetTitle("timeFit_module_calibrated [ns]")
windowTimes.GetYaxis().SetTitle("Layer")
windowTimes.Write()

windowTimes_withCut.SetTitle("Checking Timing Cut, With Cut")
windowTimes_withCut.GetXaxis().SetTitle("timeFit_module_calibrated [ns]")
windowTimes_withCut.GetYaxis().SetTitle("Layer")
windowTimes_withCut.Write()

outfile.Close()
