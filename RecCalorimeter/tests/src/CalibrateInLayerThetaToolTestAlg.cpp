/**
 * @file RecCalorimeter/tests/src/CalibrateInLayerThetaToolTestAlg.cpp
 * @brief Test for CalibrateInLayerThetaTool
 */

#undef NDEBUG

#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
#include "RecCaloCommon/ICalibrateCaloHitsTool.h"
#include "k4FWCore/GaudiChecks.h"
#include "k4Interface/IGeoSvc.h"

#include <cmath>
#include <unordered_map>
#include <vector>

namespace k4::recCalo {

class CalibrateInLayerThetaToolTestAlg : public Algorithm {
public:
  using Algorithm::Algorithm;

  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;

private:
  using CellID = dd4hep::DDSegmentation::CellID;

  CellID makeCellID(int layer, int theta) const;
  StatusCode checkClose(double value, double expected, const char* label) const;

  ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", ""};
  ToolHandle<k4::recCalo::ICalibrateCaloHitsTool> m_tool{this, "CalibrateTool", "CalibrateInLayerThetaTool", ""};
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder = nullptr;
};

DECLARE_COMPONENT(k4::recCalo::CalibrateInLayerThetaToolTestAlg);

StatusCode CalibrateInLayerThetaToolTestAlg::initialize() {
  K4_GAUDI_CHECK(m_geoSvc.retrieve());
  K4_GAUDI_CHECK(m_tool.retrieve());
  m_decoder = m_geoSvc->getDetector()->readout("ECalBarrelModuleThetaMerged").idSpec().decoder();
  return StatusCode::SUCCESS;
}

StatusCode CalibrateInLayerThetaToolTestAlg::execute() {
  std::vector<std::pair<CellID, double>> cells;
  cells.emplace_back(makeCellID(0, 100), 10.0);
  cells.emplace_back(makeCellID(1, 200), 8.0);
  cells.emplace_back(makeCellID(2, 300), 4.0);

  m_tool->calibrate(cells);
  K4_GAUDI_CHECK(checkClose(cells[0].second, 20.0, "2D SF layer 0 theta 100"));
  K4_GAUDI_CHECK(checkClose(cells[1].second, 32.0, "2D SF layer 1 theta 200"));
  K4_GAUDI_CHECK(checkClose(cells[2].second, 20.0, "layer fallback layer 2 theta 300"));

  std::unordered_map<CellID, double> cellMap;
  const CellID id = makeCellID(3, 400);
  cellMap[id] = 6.0;
  m_tool->calibrate(cellMap);
  K4_GAUDI_CHECK(checkClose(cellMap[id], 20.0, "unordered_map 2D SF layer 3 theta 400"));

  return StatusCode::SUCCESS;
}

CalibrateInLayerThetaToolTestAlg::CellID CalibrateInLayerThetaToolTestAlg::makeCellID(int layer, int theta) const {
  CellID id = 0;
  m_decoder->set(id, "system", 4);
  m_decoder->set(id, "cryo", 0);
  m_decoder->set(id, "type", 0);
  m_decoder->set(id, "subtype", 0);
  m_decoder->set(id, "layer", layer);
  m_decoder->set(id, "module", 0);
  m_decoder->set(id, "theta", theta);
  return id;
}

StatusCode CalibrateInLayerThetaToolTestAlg::checkClose(double value, double expected, const char* label) const {
  if (std::abs(value - expected) > 1e-12) {
    error() << label << ": expected " << expected << ", got " << value << endmsg;
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

} // namespace k4::recCalo
