# v00-01-00

* 2024-09-30 jmcarcell ([PR#9](https://github.com/key4hep/k4Reco/pull/9))
  - Fix the steering file for DDPlanarDigi; remove an unused property

* 2024-09-24 jmcarcell ([PR#8](https://github.com/key4hep/k4Reco/pull/8))
  - Fix histograms after Gaudi v39. RootHistogram has been renamed to StaticHistogram, see https://github.com/key4hep/k4FWCore/pull/239

* 2024-09-19 jmcarcell ([PR#2](https://github.com/key4hep/k4Reco/pull/2))
  - Add an Overlay algorithm, ported from https://github.com/iLCSoft/Overlay/blob/master/include/OverlayTiming.h

* 2024-09-09 jmcarcell ([PR#7](https://github.com/key4hep/k4Reco/pull/7))
  - Use the Key4hepConfig flag to set the standard, compiler flags and rpath magic.

* 2024-08-09 jmcarcell ([PR#6](https://github.com/key4hep/k4Reco/pull/6))
  - Change Associations to Links

* 2024-07-12 jmcarcell ([PR#5](https://github.com/key4hep/k4Reco/pull/5))
  - Don't link to GaudiAlgLib

* 2024-07-07 jmcarcell ([PR#4](https://github.com/key4hep/k4Reco/pull/4))
  - Change the default value for the time resolution to -1, which will mean no time smearing by default.
  - Use the histogram service to save histograms instead of saving our own files (which meant that if the algorithm runs multiple time there will be multiple files, or a single one if the name is the name for all of them with the content of the last one that was saved).
  - Other minor changes like removing an unnecessary header or changing `and` to `&&`.

* 2024-06-27 tmadlener ([PR#3](https://github.com/key4hep/k4Reco/pull/3))
  - Switch from TrackerHit Plane asociation to TrackerHit association, because the former has been removed from EDM4hep in [EDM4hep#331](https://github.com/key4hep/EDM4hep/pull/331)

* 2024-06-24 jmcarcell ([PR#1](https://github.com/key4hep/k4Reco/pull/1))
  - Add the DDPlanarDigi processor as a Gaudi algorithm. It should have exactly the same logic as the original DDPlanarDigiProcessor. There has been a bit of clean up, some parameters have been changed and some variable names have been changed to be in the style of the rest of the key4hep repositories. Ported from https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Digitisers/include/DDPlanarDigiProcessor.h

# v00-01

* This file is also automatically populated by the tagging script