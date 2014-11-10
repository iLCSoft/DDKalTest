




. /data/ilcsoft/HEAD/init_ilcsoft.sh
  

( if DD4hep not in ilcsoft: . ~/DD4hep/trunk/bin/thisdd4hep.sh )


cmake -C /data/ilcsoft/HEAD/ILCSoft.cmake -D DD4hep_DIR=~/DD4hep/trunk ..




-----------------------------------------------------------------------------
.  ~/DDSim/trunk/bin/thisDDSim.sh
. ../bin/thisDDKalTest.sh
