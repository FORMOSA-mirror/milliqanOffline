#This gets run after you have an estimate for the time difference between the two panels on the detector (step 2)
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
#filelist = files_from_directory("/eos/experiment/formosa/commissioning/data/offline/v34/700/0000", "Run709")[0:10]
#filelist = [file + ":t" for file in filelist]
filelist = ["/eos/experiment/formosa/commissioning/data/hadded_outputs/MilliQan_Run709_v34_throughGoingSkim.root:t"]
variables = ["timeFit", "area", "chan", "layer"]

outfile = ROOT.TFile("outputs/throughgoingBarDifferences_timeFitCalibrated.root", "RECREATE")
corrections = {0: 23.25, 1: 23.25, 2: 24.25, 3: 24.25, 4: 15.25, 5: 14.75, 6: 15.25, 7: 15.75, 8: 7.75, 9: 6.75, 10: 7.75, 11: 8.25, 12: -1.25, 13: -1.75, 14: 0.25, 15: -0.25, 16: 0, 17: 0, 18: -36.75}
timeFitCalib = ROOT.TH1F("timeFitCalibDifferences", "timeFitCalibDifferences", 200,-1e3,1e3)

panelThresholds = {16:7000, 18:4000}
barThreshold = 80000

counter = 0
for data in uproot.iterate(filelist, variables, library='ak'): 
    if counter%100==0:
        print(counter, " out of ", len(filelist), " runs")
    counter+=1
    
    data['timeFitCalib'] = ak.full_like(data['timeFit'], np.nan)
    for i in range(19):
        data['timeFitCalib'] = ak.where(data['chan']==i, data['timeFit'] + corrections[i], data['timeFitCalib'])

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
    
    throughGoingMask = ak.any(((layer0) & (barThresholdMask)), axis=-1) & ak.any(((layer1) & (barThresholdMask)), axis=-1) & ak.any(((layer2) & (barThresholdMask)), axis=-1) & ak.any(((layer3) & (barThresholdMask)), axis=-1)
    panelsMask = ak.any(((chan16) & (panelThresholdMaskFront)), axis=-1) & ak.any(((chan18) & (panelThresholdMaskBack)), axis=-1)
    #panelsMask = ak.any(((notChan16)), axis=-1) & ak.any(((notChan18)), axis=-1)

    pre_data0 = ak.drop_none(ak.where(data['area'][layer0 & throughGoingMask & panelsMask]>barThreshold, data['timeFitCalib'][layer0 & throughGoingMask & panelsMask], np.nan))
    pre_data3 = ak.drop_none(ak.where(data['area'][layer3 & throughGoingMask & panelsMask]>barThreshold, data['timeFitCalib'][layer3 & throughGoingMask & panelsMask], np.nan))

    layer0_data = ak.max(pre_data0, axis=-1)
    layer3_data = ak.max(pre_data3, axis=-1)

    x_data = ak.flatten(layer0_data - layer3_data, axis=None)

    if(len(x_data)==0):
        continue
    
    n_events = len(x_data)
    weights = ak.ones_like(x_data)
    timeFitCalib.FillN(n_events, array('d', x_data), array('d', weights), 1)
            
timeFitCalib.SetTitle("Through-Going Activity in FORMOSA After Corrections, Layer0-Layer3")
timeFitCalib.GetXaxis().SetTitle("Time Difference [ns]")
timeFitCalib.Write()
outfile.Close()