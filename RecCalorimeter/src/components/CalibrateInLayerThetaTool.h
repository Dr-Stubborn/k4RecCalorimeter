#ifndef RECCALORIMETER_CALIBRATEINLAYERTHETATOOL_H
#define RECCALORIMETER_CALIBRATEINLAYERTHETATOOL_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"

#include "DDSegmentation/BitFieldCoder.h"

// Interfaces
#include "RecCaloCommon/ICalibrateCaloHitsTool.h"
class IGeoSvc;
class TH2D;

#include <atomic>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

/** @class CalibrateInLayerThetaTool
 *
 * Tool for energy calibration to the electromagnetic scale using a two-dimensional
 * sampling-fraction map indexed by calorimeter layer and theta cell id.
 */
class CalibrateInLayerThetaTool : public extends<AlgTool, k4::recCalo::ICalibrateCaloHitsTool> {
public:
  using base_class::base_class;
  ~CalibrateInLayerThetaTool() = default;

  virtual StatusCode initialize() override final;
  virtual StatusCode finalize() override final;

  virtual void calibrate(std::unordered_map<CellID, double>& aHits) const override final;
  virtual void calibrate(std::vector<std::pair<CellID, double>>& aHits) const override final;

private:
  enum class SfSource { Map, LayerFallback, DefaultFallback, None };

  void calibrateCell(CellID cID, double& energy) const;
  double samplingFraction(CellID cID, SfSource& source) const;
  double fallbackSamplingFraction(int layerId, SfSource& source) const;
  void reportMissingBin(CellID cID, int layerId, int thetaId, const char* reason) const;
  [[noreturn]] void failMissingBin(CellID cID, int layerId, int thetaId, const char* reason) const;

  dd4hep::DDSegmentation::BitFieldCoder* m_decoder = nullptr;

  ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc"};
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "", "Name of the detector readout"};
  Gaudi::Property<std::string> m_layerFieldName{this, "layerFieldName", "layer", "Identifier of layers"};
  Gaudi::Property<std::string> m_thetaFieldName{this, "thetaFieldName", "theta", "Identifier of theta cells"};
  Gaudi::Property<int> m_firstLayerId{this, "firstLayerId", 0, "ID of first layer for fallbackSamplingFraction"};
  Gaudi::Property<std::string> m_sf2DFileName{this, "sf2DFileName", "", "ROOT file containing the layer-theta SF TH2D"};
  Gaudi::Property<std::string> m_sf2DHistName{this, "sf2DHistName", "ecal_sf_layer_theta", "Name of the layer-theta SF TH2D"};
  Gaudi::Property<std::vector<double>> m_fallbackSamplingFraction{
      this, "fallbackSamplingFraction", {}, "Fallback sampling fraction per layer"};
  Gaudi::Property<bool> m_allowLayerFallback{
      this, "allowLayerFallback", true, "Use fallbackSamplingFraction when the 2D SF map has no valid bin"};
  Gaudi::Property<bool> m_failOnMissingBin{
      this, "failOnMissingBin", false, "Fail the event processing when no valid SF is available"};
  Gaudi::Property<double> m_defaultSamplingFraction{
      this, "defaultSamplingFraction", 1.0, "Last-resort sampling fraction if layer fallback is unavailable"};
  Gaudi::Property<unsigned int> m_maxMissingBinWarnings{
      this, "maxMissingBinWarnings", 20, "Maximum number of per-cell missing-bin warnings"};

  std::unique_ptr<TH2D> m_sf2D;
  int m_layerIndex = -1;
  int m_thetaIndex = -1;

  mutable std::atomic<unsigned long long> m_mapCalibrations{0};
  mutable std::atomic<unsigned long long> m_layerFallbackCalibrations{0};
  mutable std::atomic<unsigned long long> m_defaultFallbackCalibrations{0};
  mutable std::atomic<unsigned long long> m_missingBinWarnings{0};
};

#endif /* RECCALORIMETER_CALIBRATEINLAYERTHETATOOL_H */
