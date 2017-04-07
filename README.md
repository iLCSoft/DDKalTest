# DDKalTest
[![Build Status](https://travis-ci.org/iLCSoft/DDKalTest.svg?branch=master)](https://travis-ci.org/iLCSoft/DDKalTest)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/12349/badge.svg)](https://scan.coverity.com/projects/ilcsoft-ddkaltest)

interface between KalTest fitter and DD4hep based geometry

DDKalTest is distributed under the [GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.en.html)

[![License](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Description

Re-implentation of some of the code in KalDet, now using the DDRec:Surface provided by a DD4hep based tracking geometry as input for the measurement surfaces needed in KalTest. Intersection calculation is currently done in aidaTT. Material effects use averaged material from the DDRec:Surface.

Main classes:

* cylindrical measurement layers (1D,2D):
  * DDCylinderMeasLayer
  * DDCylinderHit

* planar measurement layers (1D,2D, along and orthogonal to z):
  * DDParallelPlanarMeasLayer
  * DDPlanarMeasLayer
  * DDDiskMeasLayer ( same as DDPlanarMeasLayer)
  * DDPlanarHit

* generic detector set up from a DD4hep based tracking model
  * DDKalDetector

* measurement layer base classes: handles material effects (E-loss and MSq)
  * DDVMeasLayer

## License and Copyright
Copyright (C), DDKalTest Authors

DDKalTest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License long with this program.  If not, see <http://www.gnu.org/licenses/>.
