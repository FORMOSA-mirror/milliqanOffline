#This is STEP ONE in the timing calibration gen process
#It generates an output file with calibration histograms for each channel relative to the front panel
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
#filelist = files_from_directory("/eos/experiment/formosa/commissioning/data/offline/v34/700/0000", "Run709")
#filelist = [file + ":t" for file in filelist]
filelist = ["/eos/experiment/formosa/commissioning/data/hadded_outputs/MilliQan_Run709_v34_throughGoingSkim.root:t"]
variables = ["chan", "layer", "timeFit", "area"]

panelDifference = -34 #nanoseconds
panelThresholds = {16:7000, 18:4000}
barThreshold = 95000

outfile = ROOT.TFile("outputs/throughgoingBarDifferences.root", "RECREATE")
barCheck = ROOT.TH1F("barCheck", "barCheck", 100,0,120e3) #For verifying the bar cut

areaTimeChecks = []
barTimes = [] 
for j in range(19):
    temp_histo = ROOT.TH1F(f"barDifferencesCH{j}", f"barDifferencesCH{j}", 200,-50,50)
    temp_histo.SetDirectory(0)
    barTimes.append(temp_histo)
    
    temp_histo1 = ROOT.TH2F(f"areaVStimeCH{j}", f"areaVStimeCH{j}", 200,-50,50,100,0,120e3)
    temp_histo1.SetDirectory(0)
    areaTimeChecks.append(temp_histo1)

counter = 0
for data in uproot.iterate(filelist, variables, library='ak'): 
    if counter%100==0:
        print(counter, " out of ", len(filelist), " runs")
    counter+=1
    
    barThresholdMask = data['area']>barThreshold
    layer0 = data['layer']==0
    layer1 = data['layer']==1
    layer2 = data['layer']==2
    layer3 = data['layer']==3
    chan16 = data['chan']==16
    chan18 = data['chan']==18
    panelThresholdMaskFront = data['area']>panelThresholds[16]
    panelThresholdMaskBack = data['area']>panelThresholds[18]
    
    throughGoingMask = ak.any(((layer0) & (barThresholdMask)), axis=-1) & ak.any(((layer1) & (barThresholdMask)), axis=-1) & ak.any(((layer2) & (barThresholdMask)), axis=-1) & ak.any(((layer3) & (barThresholdMask)), axis=-1)
    panelsMask = ak.any(((chan16) & (panelThresholdMaskFront)), axis=-1) & ak.any(((chan18) & (panelThresholdMaskBack)), axis=-1)
    
    for i in range(19):
        if i==17: #No channel here
            continue
        elif i==16: #We don't want timing differences with the front panel itself...
            continue
        
        specificBarMask = data['chan']==i
        barMask = ak.any((specificBarMask)&(barThresholdMask), axis=-1)
        eventMask = throughGoingMask & panelsMask & barMask
        
        specificBar = data['chan'] == i
        frontPanel = data['chan'] == 16
        
        barData = ak.where(data['area'][eventMask & specificBar] == ak.max(data['area'][eventMask & specificBar], axis=-1), data['timeFit'][eventMask & specificBar], -1)
        barData = barData[barData>=0]
        frontData = ak.where(data['area'][eventMask & frontPanel] == ak.max(data['area'][eventMask & frontPanel], axis=-1), data['timeFit'][eventMask & frontPanel], -1)
        frontData = frontData[frontData>=0]
        
        x_data = ak.flatten(frontData - barData, axis=None)
        if(len(x_data)==0):
            continue
        
        n_events = len(x_data)
        weights = ak.ones_like(x_data)

        barTimes[i].FillN(n_events, array('d', x_data), array('d', weights), 1)
        
        #Exploring the relationship between area and these time differences
        barData_area = ak.where(data['area'][eventMask & specificBar] == ak.max(data['area'][eventMask & specificBar], axis=-1), data['area'][eventMask & specificBar], -1)
        barData_area = barData_area[barData_area>=0]
        y_data = ak.flatten(barData_area, axis=None)
        areaTimeChecks[i].FillN(n_events, array('d', x_data), array('d', y_data), array('d', weights), 1) 
        
        #Verifying the cut
        if i==0:
            barCheckData = data['area'][barThresholdMask]
            barCheckArray = ak.flatten(barCheckData, axis=None)
            barCheck.FillN(len(barCheckArray), array('d', barCheckArray), array('d', ak.ones_like(barCheckArray)), 1)
            
    
for k, bar in enumerate(barTimes):  
    bar.SetTitle(f"Front - CH{k} Times")
    bar.GetXaxis().SetTitle("Time Differences [ns]")
    bar.Write()
    
for l, check in enumerate(areaTimeChecks):  
    check.SetTitle(f"Front - CH{l} Times versus Area")
    check.GetXaxis().SetTitle("Time Differences [ns]")
    check.GetYaxis().SetTitle("Area [pVs]")
    check.Write()   

barCheck.SetTitle("Checking Area Cut on Bars")
barCheck.GetXaxis().SetTitle("Area [pVs]")
barCheck.Write()

outfile.Close()
