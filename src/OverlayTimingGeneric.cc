#include "OverlayTimingGeneric.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

#include <marlin/Global.h>
#include <marlin/ProcessorEventSeeder.h>

OverlayTimingGeneric aOverlayTimingGeneric;

OverlayTimingGeneric::OverlayTimingGeneric(): OverlayTiming("OverlayTimingGeneric")
{
  // modify processor description
  _description = "Processor to overlay events from the background taking the timing of the subdetectors into account";

  StringVec files;

  registerProcessorParameter("Delta_t",
                             "Time difference between bunches in the bunch train in ns",
                             _T_diff,
                             float(0.5));

  registerProcessorParameter("NBunchtrain",
                             "Number of bunches in a bunch train",
                             _nBunchTrain,
                             int(1));

  registerProcessorParameter("BackgroundFileNames",
                             "Name of the lcio input file(s) with background - assume one file per bunch crossing.",
                             _inputFileNames,
                             files);

  registerProcessorParameter("PhysicsBX",
                             "Number of the Bunch crossing of the physics event",
                             _BX_phys,
                             int(1));

  registerProcessorParameter("TPCDriftvelocity",
                             "[mm/ns] (float) - default 5.0e-2 (5cm/us)",
                             _tpcVdrift_mm_ns,
                             float(5.0e-2) );

  registerProcessorParameter("RandomBx",
                             "Place the physics event at an random position in the train: overrides PhysicsBX",
                             _randomBX,
                             bool(false) );

  registerProcessorParameter("NumberBackground",
                             "Number of Background events to overlay - either fixed or Poisson mean",
                             _NOverlay,
                             float(1) );

  registerProcessorParameter("Poisson_random_NOverlay",
                             "Draw random number of Events to overlay from Poisson distribution with  mean value NumberBackground",
                             _Poisson,
                             bool(false) );

  registerProcessorParameter("MergeMCParticles", 
                             "Merge the MC Particle collections",
                             _mergeMCParticles,
                             bool(true) );

  registerProcessorParameter("MCParticleCollectionName",
                             "The MC Particle Collection Name",
                             _mcParticleCollectionName,
                             std::string("MCParticle"));

  registerProcessorParameter("MCPhysicsParticleCollectionName",
                             "The output MC Particle Collection Name for the physics event" ,
                             _mcPhysicsParticleCollectionName,
                             std::string("MCPhysicsParticles"));

  //Collections with Integration Times
  registerProcessorParameter("Collection_IntegrationTimes",
                             "Integration times for the Collections",
                             _collectionTimesVec,
                             _collectionTimesVec);

  registerProcessorParameter("IntegrationTimeMin", 
                             "Lower border of the integration time window for all collections",
                             _integrationTimeMin,
                             float(-0.25) );

  registerProcessorParameter("AllowReusingBackgroundFiles",
                             "If true the same background file can be used for the same event",
                             m_allowReusingBackgroundFiles,
                             m_allowReusingBackgroundFiles);

  registerOptionalParameter("StartBackgroundFileIndex",
                            "Which background file to startWith",
                            m_startWithBackgroundFile,
                            m_startWithBackgroundFile);

  registerOptionalParameter("StartBackgroundEventIndex",
                            "Which background event to startWith",
                            m_startWithBackgroundEvent,
                            m_startWithBackgroundEvent);

    registerProcessorParameter("Start_Integration_Time",
                               "Starting integration time.  Should be shortly before the BX, but may need to be shifted earlier if the vertex is smeared in time.",
                               _DefaultStart_int,
                               float(-0.25));
}

void OverlayTimingGeneric::init()
{

  printParameters();

  overlay_Eventfile_reader = LCFactory::getInstance()->createLCReader();

  streamlog_out(WARNING) << "Attention! There are " << _inputFileNames.size()
                         << " files in the list of background files to overlay. Make sure that the total number of background events is sufficiently large for your needs!!"
                         << std::endl;

  marlin::Global::EVENTSEEDER->registerProcessor(this);

  _nRun = 0;
  _nEvt = 0;


  // parse the collectionTimesVec vector to get the collections and integration times
  std::string key;
  float value, low, high;
  size_t keyIndex = 0;
  for (size_t i = 0; i < _collectionTimesVec.size(); i++ ) {
    std::string str = _collectionTimesVec.at(i);
    if (str.find_first_not_of("-+.0123456789") != std::string::npos) {
      // Extracting the collection name
      key = str;
      keyIndex = i;
    } 
    else {
      // Extracting the 1 or 2 time values
      value = std::atof( str.c_str() );
      if (i - keyIndex == 1) {
        low = _integrationTimeMin;
        high = value;
      } else 
      if (i - keyIndex == 2) {
        low = high;
        high = value;
      }
      // Storing the time values for the last collection name to the map
      _collectionIntegrationTimes[key] = std::pair<float, float>(low, high);
    } 
  }

  streamlog_out(MESSAGE) << "Collection integration times:" << std::endl;
  for (auto const& entry : _collectionIntegrationTimes) {
    streamlog_out(MESSAGE) << "  " << entry.first << ": " << entry.second.first << " -> " << entry.second.second << std::endl;
  }


}

//------------------------------------------------------------------------------------------------------------------------------------------

void OverlayTimingGeneric::define_time_windows( std::string const& collectionName ) {

  this_start = _integrationTimeMin;
  this_stop= 0.0;
  TPC_hits = false;

  auto iter = _collectionIntegrationTimes.find( collectionName );
  if ( iter != _collectionIntegrationTimes.end() ) {
    this_start = iter->second.first;
    this_stop = iter->second.second;
  } else {
    throw std::runtime_error( "Cannot find integration time for collection " + collectionName );
  }

}
