#This script creates timing calibrations for the FORMOSA demonstrator.
#It is the first step to be run since it gives the time difference between the front and back panels
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
filelist = files_from_directory("/eos/experiment/formosa/commissioning/data/offline/v34/700/0000", "Run709")
filelist = [file + ":t" for file in filelist]
variables = ["chan", "layer", "timeFit", "area"]

panelThresholds = {16:7000, 18:4000}

outfile = ROOT.TFile("outputs/throughgoingSlabDifferences.root", "RECREATE")
slabTimes = ROOT.TH1F("slabDifferences", "slabDifferences", 200,-50,50)
areaCheckFront = ROOT.TH1F("areaCheckFront", "areaCheckFront", 100,0,25e3)
areaCheckBack = ROOT.TH1F("areaCheckBack", "areaCheckBack", 100,0,25e3)

counter = 0
for data in uproot.iterate(filelist, variables, library='ak'): 
    if counter%100==0:
        print(counter, " out of ", len(filelist), " runs")
    panelMask = ak.any(data['chan']==16, axis=-1) & ak.any(data['chan']==18, axis=-1)
    frontPanelMask = (data['area']>=panelThresholds[16]) & (data['chan']==16)
    backPanelMask = (data['area']>=panelThresholds[18]) & (data['chan']==18)

    frontMask = panelMask & frontPanelMask
    backMask = panelMask & backPanelMask
    
    frontData = ak.where(data['area'][frontMask] == ak.max(data['area'][frontMask], axis=-1), data['timeFit'][frontMask], -1)
    frontData = frontData[frontData>=0]
    backData = ak.where(data['area'][backMask] == ak.max(data['area'][backMask], axis=-1), data['timeFit'][backMask], -1)
    backData = backData[backData>=0]
    
    #frontData = ak.max(data['timeFit'][frontMask], axis=-1)
    #backData = ak.max(data['timeFit'][backMask], axis=-1)
    
    x_data = ak.flatten(frontData - backData, axis=None)
    if(len(x_data)==0):
        continue
    
    n_events = len(x_data)
    weights = ak.ones_like(x_data)

    slabTimes.FillN(n_events, array('d', x_data), array('d', weights), 1)
    
    #Verifying the cut
    frontCheck = data['area'][frontMask]
    backCheck = data['area'][backMask]
    
    frontCheckArray = ak.flatten(frontCheck, axis=None)
    backCheckArray = ak.flatten(backCheck, axis=None)
    
    areaCheckFront.FillN(len(frontCheckArray), array('d', frontCheckArray), array('d', ak.ones_like(frontCheckArray)), 1)
    areaCheckBack.FillN(len(backCheckArray), array('d', backCheckArray), array('d', ak.ones_like(backCheckArray)), 1)
    
    counter+=1
    
slabTimes.SetTitle("Front - Back Panel Times")
slabTimes.GetXaxis().SetTitle("Time Differences [ns]")
slabTimes.Write()

areaCheckFront.SetTitle("Checking Area Cut Front Panel")
areaCheckFront.GetXaxis().SetTitle("Area [pVs]")
areaCheckFront.Write()

areaCheckBack.SetTitle("Checking Area Cut Back Panel")
areaCheckBack.GetXaxis().SetTitle("Area [pVs]")
areaCheckBack.Write()

outfile.Close()