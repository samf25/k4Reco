/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "DDPlanarDigi.h"

#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

#include <GaudiKernel/ITHistSvc.h>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDSegmentation/BitFieldCoder.h"

#include "TMath.h"

#include <fmt/format.h>
#include <cmath>

DDPlanarDigi::DDPlanarDigi(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {
                           KeyValues("SimTrackerHitCollectionName", {"SimTrackerHits"}),
                           KeyValues("HeaderName", {"EventHeader"}),
                       },
                       {KeyValues("TrackerHitCollectionName", {"VTXTrackerHits"}),
                        KeyValues("SimTrkHitRelCollection", {"VTXTrackerHitRelations"})}) {
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
  }

  if (m_resULayer.size() != m_resVLayer.size()) {
    error() << "DDPlanarDigi - Inconsistent number of resolutions given for U and V coordinate: "
            << "ResolutionU  :" << m_resULayer.size() << " != ResolutionV : " << m_resVLayer.size();

    throw std::runtime_error("DDPlanarDigi: Inconsistent number of resolutions given for U and V coordinate");
  }
}

StatusCode DDPlanarDigi::initialize() {
  m_geoSvc = serviceLocator()->service(m_geoSvcName);
  if (!m_geoSvc) {
    error() << "Unable to retrieve the GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  const auto detector = m_geoSvc->getDetector();

  const auto         surfMan = detector->extension<dd4hep::rec::SurfaceManager>();
  dd4hep::DetElement det     = detector->detector(m_subDetName.value());
  surfaceMap                 = surfMan->map(m_subDetName.value());

  if (!surfaceMap) {
    throw std::runtime_error(fmt::format("Could not find surface map for detector: {} in SurfaceManager", det.name()));
  }

  // Get and store the name for a debug message later
  (void)this->getProperty("SimTrackerHitCollectionName", m_collName);

  SmartIF<ITHistSvc> histSvc;
	histSvc = serviceLocator()->service("THistSvc");

  m_histograms[hu] = new TH1F("hu", "smearing u", 50, -5., +5.);
  histSvc->regHist("/histos/" + m_collName + "/hu", m_histograms[hu]).ignore();
  m_histograms[hv] = new TH1F("hv", "smearing v", 50, -5., +5.);
  histSvc->regHist("/histos/" + m_collName + "/hv", m_histograms[hv]).ignore();
  m_histograms[hT] = new TH1F("hT", "smearing time", 50, -5., +5.);
  histSvc->regHist("/histos/" + m_collName + "/hT", m_histograms[hT]).ignore();

  m_histograms[diffu] = new TH1F("diffu", "diff u", 1000, -.1, +.1);
  histSvc->regHist("/histos/" + m_collName + "/diffu", m_histograms[diffu]).ignore();
  m_histograms[diffv] = new TH1F("diffv", "diff v", 1000, -.1, +.1);
  histSvc->regHist("/histos/" + m_collName + "/diffv", m_histograms[diffv]).ignore();
  m_histograms[diffT] = new TH1F("diffT", "diff time", 1000, -5., +5.);
  histSvc->regHist("/histos/" + m_collName + "/diffT", m_histograms[diffT]).ignore();

  m_histograms[hitE] = new TH1F("hitE", "hitEnergy in keV", 1000, 0, 200);
  histSvc->regHist("/histos/" + m_collName + "/hitE", m_histograms[hitE]).ignore();
  m_histograms[hitsAccepted] = new TH1F("hitsAccepted", "Fraction of accepted hits [%]", 201, 0, 100.5);
  histSvc->regHist("/histos/" + m_collName + "/hitsAccepted", m_histograms[hitsAccepted]).ignore();

  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> DDPlanarDigi::operator()(
    const edm4hep::SimTrackerHitCollection& simTrackerHits, const edm4hep::EventHeaderCollection& headers) const {
  auto seed = m_uidSvc->getUniqueID(headers[0].getEventNumber(), headers[0].getRunNumber(), this->name());
  debug() << "Using seed " << seed << " for event " << headers[0].getEventNumber() << " and run "
          << headers[0].getRunNumber() << endmsg;
  m_engine.SetSeed(seed);

  int nCreatedHits   = 0;
  int nDismissedHits = 0;

  auto trkhitVec = edm4hep::TrackerHitPlaneCollection();
  auto thsthcol  = edm4hep::TrackerHitSimTrackerHitLinkCollection();

  std::string cellIDEncodingString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(cellIDEncodingString);

  int nSimHits = simTrackerHits.size();
  debug() << "Processing collection " << m_collName << " with " << simTrackerHits.size() << " hits ... " << endmsg;

  for (const auto& hit : simTrackerHits) {
    m_histograms[hitE]->Fill(hit.getEDep() * (dd4hep::GeV / dd4hep::keV));

    if (hit.getEDep() < m_minEnergy) {
      debug() << "Hit with insufficient energy " << hit.getEDep() * (dd4hep::GeV / dd4hep::keV) << " keV" << endmsg;
      continue;
    }

    const int cellID0 = hit.getCellID();

    // get the measurement surface for this hit using the CellID
    dd4hep::rec::SurfaceMap::const_iterator sI = surfaceMap->find(cellID0);

    if (sI == surfaceMap->end()) {
      throw std::runtime_error(fmt::format("DDPlanarDigi::processEvent(): no surface found for cellID : {}", cellID0));
    }

    const dd4hep::rec::ISurface* surf  = sI->second;
    int                          layer = bitFieldCoder.get(cellID0, "layer");

    dd4hep::rec::Vector3D oldPos(hit.getPosition()[0], hit.getPosition()[1], hit.getPosition()[2]);
    dd4hep::rec::Vector3D newPos;

    //  Check if Hit is inside sensitive
    if (!surf->insideBounds(dd4hep::mm * oldPos)) {
      // debug() << "  hit at " << oldPos
      //         << " " << cellid_decoder( hit).valueString()
      //         << " is not on surface "
      //         << *surf
      //         << " distance: " << surf->distance(  dd4hep::mm * oldPos )
      //         << endmsg;

      if (m_forceHitsOntoSurface) {
        dd4hep::rec::Vector2D lv           = surf->globalToLocal(dd4hep::mm * oldPos);
        dd4hep::rec::Vector3D oldPosOnSurf = (1. / dd4hep::mm) * surf->localToGlobal(lv);

        debug() << " moved to " << oldPosOnSurf << " distance " << (oldPosOnSurf - oldPos).r() << endmsg;

        oldPos = oldPosOnSurf;

      } else {
        ++nDismissedHits;
        continue;
      }
    }

    // Smear time of the hit and apply the time window cut if needed
    double hitT = hit.getTime();

    if (m_resTLayer.size() && m_resTLayer[0] > 0) {
      float resT = m_resTLayer.size() > 1 ? m_resTLayer[layer] : m_resTLayer[0];

      double tSmear = resT > 0 ? m_engine.Gaus(0, resT) : 0;
      m_histograms[hT]->Fill(resT > 0 ? tSmear / resT : 0);
      m_histograms[diffT]->Fill(tSmear);

      hitT += tSmear;
      debug() << "smeared hit at T: " << hit.getTime() << " ns to T: " << hitT
              << " ns according to resolution: " << resT << " ns" << endmsg;
    }

    // Correcting for the propagation time
    if (m_correctTimesForPropagation) {
      double dt = oldPos.r() / (TMath::C() / 1e6);
      hitT -= dt;
      debug() << "corrected hit at R: " << oldPos.r() << " mm by propagation time: " << dt << " ns to T: " << hitT
              << " ns" << endmsg;
    }

    // Skip the hit if its time is outside the acceptance time window
    if (m_useTimeWindow) {
      // TODO: Check that the length of the time window is OK
      float timeWindow_min = m_timeWindowMin.size() > 1 ? m_timeWindowMin[layer] : m_timeWindowMin[0];
      float timeWindow_max = m_timeWindowMax.size() > 1 ? m_timeWindowMax[layer] : m_timeWindowMax[0];
      if (hitT < timeWindow_min || hitT > timeWindow_max) {
        debug() << "hit at T: " << hit.getTime() << " smeared to: " << hitT
                << " is outside the time window: hit dropped" << endmsg;
        ++nDismissedHits;
        continue;
      }
    }

    // Try to smear the hit position but ensure the hit is inside the sensitive region
    dd4hep::rec::Vector3D u = surf->u();
    dd4hep::rec::Vector3D v = surf->v();

    // Get local coordinates on surface
    dd4hep::rec::Vector2D lv = surf->globalToLocal(dd4hep::mm * oldPos);
    double                uL = lv[0] / dd4hep::mm;
    double                vL = lv[1] / dd4hep::mm;

    bool acceptHit = false;
    int  tries     = 0;

    // TODO: check lengths
    float resU = m_resULayer.size() > 1 ? m_resULayer[layer] : m_resULayer[0];
    float resV = m_resVLayer.size() > 1 ? m_resVLayer[layer] : m_resVLayer[0];

    while (tries < m_maxTries) {
      // if( tries > 0 ) debug() << "retry smearing for " <<  cellid_decoder( hit ).valueString() << " : retries " << tries << endmsg;

      double uSmear = m_engine.Gaus(0, resU);
      double vSmear = m_engine.Gaus(0, resV);

      dd4hep::rec::Vector3D newPosTmp;
      if (m_isStrip) {
        if (m_subDetName == "SET") {
          double xStripPos, yStripPos, zStripPos;
          //Find intersection of the strip with the z=centerOfSensor plane to set it as the center of the SET strip
          dd4hep::rec::Vector3D simHitPosSmeared =
              (1. / dd4hep::mm) * (surf->localToGlobal(dd4hep::rec::Vector2D((uL + uSmear) * dd4hep::mm, 0.)));
          zStripPos        = surf->origin()[2] / dd4hep::mm;
          double lineParam = (zStripPos - simHitPosSmeared[2]) / v[2];
          xStripPos        = simHitPosSmeared[0] + lineParam * v[0];
          yStripPos        = simHitPosSmeared[1] + lineParam * v[1];
          newPosTmp        = dd4hep::rec::Vector3D(xStripPos, yStripPos, zStripPos);
        } else {
          newPosTmp = (1. / dd4hep::mm) * (surf->localToGlobal(dd4hep::rec::Vector2D((uL + uSmear) * dd4hep::mm, 0.)));
        }
      } else {
        newPosTmp =
            (1. / dd4hep::mm) *
            (surf->localToGlobal(dd4hep::rec::Vector2D((uL + uSmear) * dd4hep::mm, (vL + vSmear) * dd4hep::mm)));
      }

      debug() << " hit at    : " << oldPos << " smeared to: " << newPosTmp << " uL: " << uL << " vL: " << vL
              << " uSmear: " << uSmear << " vSmear: " << vSmear << endmsg;

      if (surf->insideBounds(dd4hep::mm * newPosTmp)) {
        acceptHit = true;
        newPos    = newPosTmp;

        m_histograms[hu]->Fill(uSmear / resU);
        m_histograms[hv]->Fill(vSmear / resV);

        m_histograms[diffu]->Fill(uSmear);
        m_histograms[diffv]->Fill(vSmear);

        break;
      }

      // debug() << "  hit at " << newPosTmp
      //         << " " << cellid_decoder( hit).valueString()
      //         << " is not on surface "
      //         << " distance: " << surf->distance( dd4hep::mm * newPosTmp )
      //         << endmsg;

      ++tries;
    }

    if (!acceptHit) {
      debug() << "hit could not be smeared within ladder after " << m_maxTries << "  tries: hit dropped" << endmsg;
      ++nDismissedHits;
      continue;
    }

    auto trkHit = trkhitVec.create();

    trkHit.setCellID(cellID0);

    trkHit.setPosition(newPos.const_array());
    trkHit.setTime(hitT);
    trkHit.setEDep(hit.getEDep());

    float u_direction[2];
    u_direction[0] = u.theta();
    u_direction[1] = u.phi();

    float v_direction[2];
    v_direction[0] = v.theta();
    v_direction[1] = v.phi();

    debug() << " U[0] = " << u_direction[0] << " U[1] = " << u_direction[1] << " V[0] = " << v_direction[0]
            << " V[1] = " << v_direction[1] << endmsg;

    trkHit.setU(u_direction);
    trkHit.setV(v_direction);
    trkHit.setDu(resU);

    if (m_isStrip) {
      // store the resolution from the length of the wafer - in case a fitter might want to treat this as 2d hit ....
      double stripRes = surf->length_along_v() / dd4hep::mm / std::sqrt(12);
      trkHit.setDv(stripRes);
      // TODO: Set type?
      // trkHit.setType( UTIL::set_bit( trkHit.getType() ,  UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ) );

    } else {
      trkHit.setDv(resV);
    }

    auto association = thsthcol.create();
    association.setTo(hit);
    association.setFrom(trkHit);

    ++nCreatedHits;
  }

  // Filling the fraction of accepted hits in the event
  float accFraction = nSimHits > 0 ? float(nCreatedHits) / float(nSimHits) * 100.0 : 0.0;
  m_histograms[hitsAccepted]->Fill(accFraction);

  debug() << "Created " << nCreatedHits << " hits, " << nDismissedHits << " hits  dismissed" << endmsg;

  return std::make_tuple(std::move(trkhitVec), std::move(thsthcol));
}
