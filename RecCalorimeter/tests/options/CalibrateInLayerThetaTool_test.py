#
# File: RecCalorimeter/tests/options/CalibrateInLayerThetaTool_test.py
# Purpose: Test for CalibrateInLayerThetaTool
#

import os

import Configurables as C
import ROOT
from RecCalorimeter.k4RecCalorimeterPluginsConf import CalibrateInLayerThetaTool
from RecCalorimeter.k4RecCalorimeterTestsConf import (
    k4__recCalo__CalibrateInLayerThetaToolTestAlg,
)


sf_file = "calibrate_layer_theta_sf_test.root"
sf_hist_name = "ecal_sf_layer_theta"

root_file = ROOT.TFile(sf_file, "RECREATE")
sf_hist = ROOT.TH2D(sf_hist_name, sf_hist_name, 11, -0.5, 10.5, 1024, -0.5, 1023.5)
sf_hist.SetBinContent(sf_hist.GetXaxis().FindBin(0), sf_hist.GetYaxis().FindBin(100), 0.5)
sf_hist.SetBinContent(sf_hist.GetXaxis().FindBin(1), sf_hist.GetYaxis().FindBin(200), 0.25)
sf_hist.SetBinContent(sf_hist.GetXaxis().FindBin(3), sf_hist.GetYaxis().FindBin(400), 0.3)
sf_hist.Write()
root_file.Close()

compact_file = "ALLEGRO_o1_v03.xml"
path_to_detector = (
    os.environ.get("K4GEO", "")
    + "/FCCee/ALLEGRO/compact/"
    + os.path.splitext(compact_file)[0]
)

geoSvc = C.GeoSvc("GeoSvc", detectors=[os.path.join(path_to_detector, compact_file)])

calibrateTool = CalibrateInLayerThetaTool(
    "CalibrateECalBarrelLayerThetaTest",
    readoutName="ECalBarrelModuleThetaMerged",
    layerFieldName="layer",
    thetaFieldName="theta",
    sf2DFileName=sf_file,
    sf2DHistName=sf_hist_name,
    fallbackSamplingFraction=[0.1, 0.2, 0.2, 0.3],
    allowLayerFallback=True,
    failOnMissingBin=False,
    maxMissingBinWarnings=5,
)

appmgr = C.ApplicationMgr(
    TopAlg=[
        k4__recCalo__CalibrateInLayerThetaToolTestAlg(CalibrateTool=calibrateTool)
    ],
    ExtSvc=[geoSvc],
    EvtMax=1,
)
