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
filelist = ["/eos/experiment/formosa/commissioning/data/hadded_outputs/MilliQan_Run737_incomplete.root:t"]
#filelist = ["/eos/experiment/formosa/commissioning/data/hadded_outputs/MilliQan_Run709_v35_matchedANDprocessed.root:t"]
#filelist = ["/eos/experiment/formosa/commissioning/data/hadded_outputs/MilliQan_Run740_v35_matchedANDprocessed.root:t"]
variables = ["timeFit", "area", "chan", "layer", "timeFit_module_calibrated", "height", "event_time_fromTDC", "tTrigger"]

timing_range = [1714.04e6, 1714.07e6]#uproot_range_finder(filelist[0], "event_time_fromTDC") #[1712.48e6, 1712.51e6]#
outfile = ROOT.TFile("outputs/exploringWindowTiming_run737_highHV_caseComparison.root", "RECREATE")

HV_setting = "high" #"high"

if HV_setting == "low":
    panelThresholds = {16:7000, 18:4000}
elif HV_setting == "high":
    panelThresholds = {16:30000, 18:18000}
else:
    raise ValueError("Wrong HV Setting Selected!")

dataSets = {}
counter = 0
counter_limit = None
for data in uproot.iterate(filelist, variables, library='ak', step_size=100000): 
    if counter%100==0:
        print(counter, " out of ", len(filelist), " runs")  
          
    if counter_limit is not None and counter==counter_limit: 
        continue
    counter+=1

    print("NEW STEP")
    
    #Broadcast event_tibame_fromTDC to proper shape
    data['event_time_fromTDC'] = ak.broadcast_arrays(ak.flatten(ak.pad_none(data['event_time_fromTDC'], 1, axis=-1), axis=None), data['area'])[0]
    data['tTrigger'] = ak.broadcast_arrays(ak.flatten(ak.pad_none(data['tTrigger'], 1, axis=-1), axis=None), data['area'])[0]

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
    locationMask = throughGoingMask & ak.all(notChan16, axis=-1) & ak.all(notChan18, axis=-1) 
    locationMask_withCut = throughGoingMask_withCut & ak.all(notChan16, axis=-1) & ak.all(notChan18, axis=-1) 
    locationMask_withPanels = throughGoingMask & ak.any(chan16, axis=-1) & ak.any(chan18, axis=-1) &  ak.all((chan16 & (~panelThresholdMaskFront)) | (chan18 & (~panelThresholdMaskBack)) | (notChan16 & notChan18), axis=-1)
    locationMask_withCut_withPanels = throughGoingMask_withCut & ak.any(chan16, axis=-1) & ak.any(chan18, axis=-1) &  ak.all((chan16 & (~panelThresholdMaskFront)) | (chan18 & (~panelThresholdMaskBack)) | (notChan16 & notChan18), axis=-1)
    locationMask_withPanels_threshold = throughGoingMask & ak.any(chan16 & panelThresholdMaskFront, axis=-1) & ak.any(chan18 & panelThresholdMaskBack, axis=-1) 
    locationMask_withCut_withPanels_threshold = throughGoingMask_withCut & ak.any(chan16 & panelThresholdMaskFront, axis=-1) & ak.any(chan18 & panelThresholdMaskBack, axis=-1) 
    
    onePanelHitOnlyMask_withCut = throughGoingMask_withCut & ~((ak.any(chan16, axis=-1) & ak.any(chan18, axis=-1))) & (ak.any(chan16, axis=-1) | ak.any(chan18, axis=-1)) #To look for events where there's only a pulse in one of the panels
    onePanelAboveThreshold_withCut = throughGoingMask_withCut & (ak.any(chan16, axis=-1) & ak.any(chan18, axis=-1)) & (ak.any(chan16 & panelThresholdMaskFront, axis=-1) | ak.any(chan18 & panelThresholdMaskBack, axis=-1)) & ~(ak.any(chan16 & panelThresholdMaskFront, axis=-1) & ak.any(chan18 & panelThresholdMaskBack, axis=-1)) #To look for events where there's only a pulse above threshold in one panel, but both panels have hits

    #When adding a new histo to generate, add it here and in the else section below
    if counter == 1:
        dataCopy_throughGoing = data[throughGoingMask] #For all the through-going stuff, no timing cut
        dataSets['dataCopy_throughGoing'] = [dataCopy_throughGoing, "All through-going events, no timing cut", "throughGoing"]
        
        dataCopy_throughGoing_withCut = data[throughGoingMask_withCut] #For all the through-going stuff, loose timing cut
        dataSets['dataCopy_throughGoing_withCut'] = [dataCopy_throughGoing_withCut, "All through-going events, loose timing cut", "throughGoing_withCut"]
        
        dataCopy = data[locationMask] #Through-going, no panels, no time cut
        dataSets['dataCopy'] = [dataCopy, "Through-going, no panels, no time cut", ""]
        
        dataCopy_withCut = data[locationMask_withCut] #Through-going, no panels, loose time cut
        dataSets['dataCopy_withCut'] = [dataCopy_withCut, "Through-going, no panels, loose time cut", "withCut"]
        
        dataCopy_withPanels = data[locationMask_withPanels] #Through-going, with panels, no time cut, all below panel threshold
        dataSets['dataCopy_withPanels'] = [dataCopy_withPanels, "Through-going, with panels, no time cut, all below panel threshold", "withPanels"]
        
        dataCopy_withCut_withPanels = data[locationMask_withCut_withPanels] #Through-going, with panels, loose time cut, all below panel threshold
        dataSets['dataCopy_withCut_withPanels'] = [dataCopy_withCut_withPanels, "Through-going, with panels, loose time cut, all below panel threshold", "withCut_withPanels"]
        
        dataCopy_withPanels_threshold = data[locationMask_withPanels_threshold] #Through-going, with panels, no time cut, one per panel above threshold
        dataSets['dataCopy_withPanels_threshold'] = [dataCopy_withPanels_threshold, "Through-going, with panels, no time cut, one per panel above threshold", "withPanels_threshold"]
        
        dataCopy_withCut_withPanels_threshold = data[locationMask_withCut_withPanels_threshold] #Through-going, with panels, loose time cut, one per panel above threshold
        dataSets['dataCopy_withCut_withPanels_threshold'] = [dataCopy_withCut_withPanels_threshold, "Through-going, with panels, loose time cut, one per panel above threshold", "withCut_withPanels_threshold"]
        
        dataCopy_withCut_onePanel = data[onePanelHitOnlyMask_withCut] #Through-going, with one panels, loose time cut
        dataSets['dataCopy_withCut_onePanel'] = [dataCopy_withCut_onePanel, "Through-going, with one panel, loose time cut", "withCut_onePanel"]
        
        dataCopy_withCut_onePanel_threshold = data[onePanelAboveThreshold_withCut] #Through-going, with two panels (only one above threshold), loose time cut
        dataSets['dataCopy_withCut_onePanel_threshold'] = [dataCopy_withCut_onePanel_threshold, "Through-going, with two panels (only one above threshold), loose time cut", "withCut_onePanel_threshold"]
        
        leftoverData_withCut = data[throughGoingMask_withCut & (~locationMask_withCut) & (~locationMask_withCut_withPanels) & (~locationMask_withCut_withPanels_threshold) & (~onePanelHitOnlyMask_withCut) & (~onePanelAboveThreshold_withCut)] #Tells us what through-going stuff is missed by with the panel selections
        dataSets['leftoverData_withCut'] = [leftoverData_withCut, "Through-going stuff is missed by with the panel selections", "leftover"]
        
    else:
        dataSets['dataCopy_throughGoing'][0] = data[throughGoingMask]
        dataSets['dataCopy_throughGoing_withCut'][0] = data[throughGoingMask_withCut]
        dataSets['dataCopy'][0] = data[locationMask]
        dataSets['dataCopy_withCut'][0] = data[locationMask_withCut]
        dataSets['dataCopy_withPanels'][0] = data[locationMask_withPanels]
        dataSets['dataCopy_withCut_withPanels'][0] = data[locationMask_withCut_withPanels]
        dataSets['dataCopy_withPanels_threshold'][0] = data[locationMask_withPanels_threshold]
        dataSets['dataCopy_withCut_withPanels_threshold'][0] = data[locationMask_withCut_withPanels_threshold]
        dataSets['dataCopy_withCut_onePanel'][0] = data[onePanelHitOnlyMask_withCut]
        dataSets['dataCopy_withCut_onePanel_threshold'][0] = data[onePanelAboveThreshold_withCut]
        dataSets['leftoverData_withCut'][0] = data[throughGoingMask_withCut & (~locationMask_withCut) & (~locationMask_withCut_withPanels) & (~locationMask_withCut_withPanels_threshold) & (~onePanelHitOnlyMask_withCut) & (~onePanelAboveThreshold_withCut)]
        
    
    for key in dataSets.keys():
        iDataset = dataSets[key][0]
        maxes = iDataset['area']==ak.max(iDataset['area'], axis=-1) #We only want one pulse per event
        maxedTimes = ak.flatten(iDataset['event_time_fromTDC'][maxes], axis=None) #Get one event time per event
        maxedTriggers = ak.flatten(iDataset['tTrigger'][maxes], axis=None)
        n_events = len(maxedTimes)
        weights = ak.ones_like(maxedTimes)
        
        if len(maxedTimes) != len(maxedTriggers):
            raise ValueError("Triggers and Times have two different lengths :(")
        
        if counter==1:
            #Initializing the histo objects
            dataSets[key].append(ROOT.TH1D(f"times_{dataSets[key][2]}", f"times_{dataSets[key][2]}", 500, timing_range[0], timing_range[1])) #event_time_fromTDC
            dataSets[key][3].SetDirectory(0)
            dataSets[key].append(ROOT.TH1D(f"triggerBits_{dataSets[key][2]}", f"triggerBits_{dataSets[key][2]}", 100, 0, 100)) #Trigger Bits fired
            dataSets[key][4].SetDirectory(0)
            dataSets[key].append(ROOT.TH1D(f"rates_{dataSets[key][2]}", f"rates_{dataSets[key][2]}", 500, timing_range[0], timing_range[1])) #calculated rates
            dataSets[key][5].SetDirectory(0)
            
        if n_events==0:
            print(f"{str(key)} has no entries for counter {counter}")
            continue

        #print(n_events, array('d', maxedTimes)[0:10], array('d', weights)[0:10], 1)
        #print(n_events, array('d', maxedTriggers)[0:10], array('d', weights)[0:10], 1)
        
        dataSets[key][3].FillN(n_events, array('d', maxedTimes), array('d', weights), 1)
        dataSets[key][4].FillN(n_events, array('d', maxedTriggers), array('d', weights), 1)
        

for filled_key in dataSets.keys():
    print(f"Writing {str(filled_key)}")
    filled_array = dataSets[filled_key]
    times = filled_array[3]
    triggers = filled_array[4]
    rates = filled_array[5]
    
    times.SetTitle(filled_array[1])
    times.GetXaxis().SetTitle("Event Time [ns]")
    times.Write()

    rates.Add(times)
    rates.GetXaxis().SetTimeDisplay(1)
    rates.GetXaxis().SetLabelSize(0.02)
    rates.GetXaxis().SetLabelOffset(0.01)
    rates.GetXaxis().SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 02:00:00")

    width = rates.GetBinWidth(30)
    inv_width = 1/width
    rates.Scale(inv_width)

    rates.SetTitle(filled_array[1])
    rates.GetXaxis().SetTitle("Event Time [ns]")
    rates.GetYaxis().SetTitle("Rate [Hz]")
    rates.Write()

    triggers.SetTitle(filled_array[1])
    triggers.GetXaxis().SetTitle("Trigger Bits (Integer-Binned)")
    triggers.GetYaxis().SetTitle("Counts")
    triggers.Write()

outfile.Close()
