# v01-01
A.Sailer
* DDKalTest::CellIDEncoding: implement set_encoding_string and add access protection.

F.Gaede
* fix treatment of ROOT in CMakeLists.txt

# v01-00-01
A.Sailer
* DDParallelPlanarMeasLayer: protect against basically infinite loop if phi value is very very large
* DDKaltest: Ignore warnings from external header files
 
# v01-00
F.Gaede
* fixed the sorting policy for zdisks -> major improvement in forward track fitting
* crosscheck of energy loss from aidaTT commented out
* made compatible with c++11 and ROOT6
* added DDSurfaces::ISurface* surface() and some debug output
* remove hard coded muon mass for track fits

# v00-03
F.Gaede
* use new (faster) standalone helix intersection methods
* workaround for hits that are smeared off the DDRec::surface in the old digitizer
* added DDConeMeasLayer
* fixed bug in energy loss and multiple scattering for non-planar surfaces ( compute normal at x'ing point)
* distinguish between zplanes, disks and generic planes for sorting policy
* switch back to use SurfaceHelper class in order to get all insensitive surfaces ( lost in SurfaceMap )
* fixed some compiler warnings
 
# v01-00

F.Gaede
* added DDSurfaces::ISurface* surface() and some debug output 
* remove hard coded muon mass for track fits 
* fixed the sorting policy for zdisks -> major improvement in forward track fitting
* crosscheck of energy loss from aidaTT
* commented out 
* made compatible with c++11 and ROOT6


# v00-03
 F.Gaede
* adjusted some loge levels for debug ...
* use new (faster) standalone helix intersection methods
* workaround for hits that are smeared off the DDRec::surface in the old digitizer
* added DDConeMeasLayer
* fixed bug in energy loss and multiple scattering for non-planar surfaces ( compute normal at x'ing point)
* distinguish between zplanes, disks and generic planes for sorting policy
* replaced WARNING w/ DEBUG for debug printout
* switch back to use SurfaceHelper class in order to get all insensitive suirfaces ( lost in SurfaceMap )
* fixed some compiler warnings

# v00-02
* changed to just use abstract ISurface and ICylinder

# v00-01
* first release with basic functionality
