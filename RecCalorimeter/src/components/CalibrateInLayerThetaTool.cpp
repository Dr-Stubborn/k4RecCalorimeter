#include "CalibrateInLayerThetaTool.h"

// k4FWCore
#include "RecCaloCommon/k4RecCalorimeter_check.h"
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"

// Gaudi
#include "GaudiKernel/GaudiException.h"

// ROOT
#include "TAxis.h"
#include "TFile.h"
#include "TH2D.h"

#include <cmath>

DECLARE_COMPONENT(CalibrateInLayerThetaTool)

StatusCode CalibrateInLayerThetaTool::initialize() {
  K4RECCALORIMETER_CHECK(AlgTool::initialize());

  if (m_geoSvc->getDetector()->readouts().find(m_readoutName) == m_geoSvc->getDetector()->readouts().end()) {
    error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  m_layerIndex = m_decoder->index(m_layerFieldName);
  m_thetaIndex = m_decoder->index(m_thetaFieldName);
  if (m_layerIndex < 0 || m_thetaIndex < 0) {
    error() << "Readout <<" << m_readoutName << ">> does not contain required fields: " << m_layerFieldName << ", "
            << m_thetaFieldName << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_sf2DFileName.empty()) {
    error() << "Property sf2DFileName is empty." << endmsg;
    return StatusCode::FAILURE;
  }
  auto inputFile = std::unique_ptr<TFile>(TFile::Open(m_sf2DFileName.value().c_str(), "READ"));
  if (!inputFile || inputFile->IsZombie()) {
    error() << "Cannot open layer-theta SF file: " << m_sf2DFileName << endmsg;
    return StatusCode::FAILURE;
  }
  auto* object = inputFile->Get(m_sf2DHistName.value().c_str());
  auto* sfHist = dynamic_cast<TH2D*>(object);
  if (sfHist == nullptr) {
    error() << "Cannot find TH2D <<" << m_sf2DHistName << ">> in file " << m_sf2DFileName << endmsg;
    return StatusCode::FAILURE;
  }
  m_sf2D.reset(static_cast<TH2D*>(sfHist->Clone((name() + "_sf2D").c_str())));
  m_sf2D->SetDirectory(nullptr);

  if (m_allowLayerFallback && m_fallbackSamplingFraction.empty()) {
    warning() << "Layer fallback is enabled but fallbackSamplingFraction is empty. "
              << "Missing 2D SF bins will use defaultSamplingFraction=" << m_defaultSamplingFraction << endmsg;
  }
  if (m_defaultSamplingFraction.value() <= 0.0 || !std::isfinite(m_defaultSamplingFraction.value())) {
    error() << "defaultSamplingFraction must be finite and positive." << endmsg;
    return StatusCode::FAILURE;
  }

  info() << "Using layer-theta SF map " << m_sf2DFileName << ":" << m_sf2DHistName << " for readout "
         << m_readoutName << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode CalibrateInLayerThetaTool::finalize() {
  info() << "Layer-theta SF calibration summary: map=" << m_mapCalibrations
         << ", layerFallback=" << m_layerFallbackCalibrations
         << ", defaultFallback=" << m_defaultFallbackCalibrations << endmsg;
  return AlgTool::finalize();
}

void CalibrateInLayerThetaTool::calibrate(std::unordered_map<CellID, double>& aHits) const {
  for (auto& p : aHits) {
    calibrateCell(p.first, p.second);
  }
}

void CalibrateInLayerThetaTool::calibrate(std::vector<std::pair<CellID, double>>& aHits) const {
  for (auto& p : aHits) {
    calibrateCell(p.first, p.second);
  }
}

void CalibrateInLayerThetaTool::calibrateCell(CellID cID, double& energy) const {
  SfSource source = SfSource::None;
  const double sf = samplingFraction(cID, source);
  if (sf > 0.0 && std::isfinite(sf)) {
    energy /= sf;
  }
  if (source == SfSource::Map) {
    ++m_mapCalibrations;
  } else if (source == SfSource::LayerFallback) {
    ++m_layerFallbackCalibrations;
  } else if (source == SfSource::DefaultFallback) {
    ++m_defaultFallbackCalibrations;
  }
}

double CalibrateInLayerThetaTool::samplingFraction(CellID cID, SfSource& source) const {
  const int layerId = m_decoder->get(cID, m_layerIndex);
  const int thetaId = m_decoder->get(cID, m_thetaIndex);

  const int layerBin = m_sf2D->GetXaxis()->FindFixBin(layerId);
  const int thetaBin = m_sf2D->GetYaxis()->FindFixBin(thetaId);
  const bool inRange = layerBin >= 1 && layerBin <= m_sf2D->GetNbinsX() && thetaBin >= 1 &&
                       thetaBin <= m_sf2D->GetNbinsY();
  if (inRange) {
    const double sf = m_sf2D->GetBinContent(layerBin, thetaBin);
    if (sf > 0.0 && std::isfinite(sf)) {
      source = SfSource::Map;
      return sf;
    }
    reportMissingBin(cID, layerId, thetaId, "invalid or non-positive SF");
    return fallbackSamplingFraction(layerId, source);
  }

  reportMissingBin(cID, layerId, thetaId, "outside SF map axes");
  return fallbackSamplingFraction(layerId, source);
}

double CalibrateInLayerThetaTool::fallbackSamplingFraction(int layerId, SfSource& source) const {
  if (m_allowLayerFallback) {
    const int layerIndex = layerId - m_firstLayerId;
    if (layerIndex >= 0 && static_cast<size_t>(layerIndex) < m_fallbackSamplingFraction.size()) {
      const double sf = m_fallbackSamplingFraction[layerIndex];
      if (sf > 0.0 && std::isfinite(sf)) {
        source = SfSource::LayerFallback;
        return sf;
      }
    }
  }

  if (m_failOnMissingBin) {
    failMissingBin(0, layerId, -1, "no valid 2D SF or fallback SF");
  }
  source = SfSource::DefaultFallback;
  return m_defaultSamplingFraction;
}

void CalibrateInLayerThetaTool::reportMissingBin(CellID cID, int layerId, int thetaId, const char* reason) const {
  if (m_failOnMissingBin) {
    failMissingBin(cID, layerId, thetaId, reason);
  }
  const auto warningIndex = ++m_missingBinWarnings;
  if (warningIndex <= m_maxMissingBinWarnings.value()) {
    warning() << "Missing layer-theta SF bin: layer=" << layerId << ", theta=" << thetaId << ", cellID=" << cID
              << ", reason=" << reason << ". Applying configured fallback." << endmsg;
  } else if (warningIndex == static_cast<unsigned long long>(m_maxMissingBinWarnings.value()) + 1) {
    warning() << "Further missing layer-theta SF bin warnings are suppressed." << endmsg;
  }
}

void CalibrateInLayerThetaTool::failMissingBin(CellID cID, int layerId, int thetaId, const char* reason) const {
  throw GaudiException("Missing layer-theta SF bin for cellID=" + std::to_string(cID) +
                           ", layer=" + std::to_string(layerId) + ", theta=" + std::to_string(thetaId) +
                           ", reason=" + reason,
                       name(), StatusCode::FAILURE);
}
