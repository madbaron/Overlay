#include "OverlayTimingRandomMix.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/SimTrackerHit.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/MCParticleImpl.h>
#include <IMPL/SimCalorimeterHitImpl.h>
#include <IMPL/SimTrackerHitImpl.h>

#include <marlin/Exceptions.h>
#include <marlin/Global.h>
#include <marlin/ProcessorEventSeeder.h>

#include <algorithm>
#include <limits>
#include <random>
#include <set>
#include <iterator>
#include <filesystem>

using namespace lcio;
using namespace marlin;
namespace fs = std::filesystem;

void listFilesInFolder(const std::string &folderPath, StringVec &files)
{
  try
  {
    for (const auto &entry : fs::directory_iterator(folderPath))
    {
      if (fs::is_regular_file(entry.path()))
      {
        files.push_back(entry.path().string());
      }
    }
  }
  catch (const std::exception &ex)
  {
    std::cerr << "Error: " << ex.what() << std::endl;
  }
}

namespace overlay
{

  OverlayTimingRandomMix aOverlayTimingRandomMix;

  OverlayTimingRandomMix::OverlayTimingRandomMix(std::string const &procName) : Processor(procName)
  {
  }

  OverlayTimingRandomMix::OverlayTimingRandomMix() : Processor("OverlayTimingRandomMix")
  {
    // modify processor description
    _description = "Processor to overlay events from the background taking the timing of the subdetectors into account";

    registerProcessorParameter("PathToMuPlus",
                               "Path to the fluka lcio input file(s)",
                               _pathToMuPlus,
                               std::string("/data/BIB10TeV/sim_mp"));

    registerProcessorParameter("PathToMuMinus",
                               "Path to the fluka lcio input file(s)",
                               _pathToMuMinus,
                               std::string("/data/BIB10TeV/sim_mp"));

    registerProcessorParameter("NumberBackground",
                               "Number of Background events to overlay",
                               _NOverlay,
                               float(1));

    registerProcessorParameter("MergeMCParticles",
                               "Merge the MC Particle collections",
                               _mergeMCParticles,
                               bool(true));

    registerProcessorParameter("MCParticleCollectionName",
                               "The MC Particle Collection Name",
                               _mcParticleCollectionName,
                               std::string("MCParticle"));

    registerProcessorParameter("Collection_IntegrationTimes",
                               "Integration times for the Collections",
                               _collectionTimesVec,
                               _collectionTimesVec);

    registerProcessorParameter("IntegrationTimeMin",
                               "Lower border of the integration time window for all collections",
                               _integrationTimeMin,
                               float(-0.5));
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void OverlayTimingRandomMix::init()
  {
    streamlog_out(DEBUG) << " init called  " << std::endl;
    printParameters();

    overlay_Eventfile_reader = LCFactory::getInstance()->createLCReader();

    marlin::Global::EVENTSEEDER->registerProcessor(this);

    // parse the collectionTimesVec vector to get the collections and integration times
    std::string key;
    float value, low, high;
    size_t keyIndex = 0;
    for (size_t i = 0; i < _collectionTimesVec.size(); i++)
    {
      std::string str = _collectionTimesVec.at(i);
      if (str.find_first_not_of("-+.0123456789") != std::string::npos)
      {
        // Extracting the collection name
        key = str;
        keyIndex = i;
      }
      else
      {
        // Extracting the 1 or 2 time values
        value = std::atof(str.c_str());
        if (i - keyIndex == 1)
        {
          low = _integrationTimeMin;
          high = value;
        }
        else if (i - keyIndex == 2)
        {
          low = high;
          high = value;
        }
        // Storing the time values for the last collection name to the map
        _collectionIntegrationTimes[key] = std::pair<float, float>(low, high);
      }
    }

    // Get the list of files from the paths
    listFilesInFolder(_pathToMuMinus, _inputFileNamesMuMinus);
    listFilesInFolder(_pathToMuPlus, _inputFileNamesMuPlus);

    if ((_NOverlay > _inputFileNamesMuPlus.size()) || (_NOverlay > _inputFileNamesMuMinus.size()))
    {
      streamlog_out(WARNING) << "Attention! There are " << _inputFileNamesMuPlus.size() << " (" << _inputFileNamesMuMinus.size()
                             << ") files in the list of mu plus (minus) background files to overlay. Make sure that each background file list is equal or greater to NumberBackground!!"
                             << std::endl;
    }

    streamlog_out(MESSAGE) << "Collection integration times:" << std::endl;
    for (auto const &entry : _collectionIntegrationTimes)
    {
      streamlog_out(MESSAGE) << "  " << entry.first << ": " << entry.second.first << " -> " << entry.second.second << std::endl;
    }

    _nRun = 0;
    _nEvt = 0;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void OverlayTimingRandomMix::processRunHeader(EVENT::LCRunHeader *)
  {
    _nRun++;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void OverlayTimingRandomMix::modifyEvent(EVENT::LCEvent *evt)
  {

    std::random_device rd;
    std::mt19937 g(rd());

    // We have the physics event in evt. Now we merge the new overlay events with it.
    // First cut the collections in the physics event to the defined time windows
    const std::vector<std::string> *collection_names_in_Evt = evt->getCollectionNames();

    for (unsigned int j = 0, nCollections = collection_names_in_Evt->size(); j < nCollections; ++j)
    {
      const std::string Collection_name = collection_names_in_Evt->at(j);
      LCCollection *Collection_in_Physics_Evt = evt->getCollection(Collection_name);
      currentDest = Collection_name;
      if ((Collection_in_Physics_Evt->getTypeName() == LCIO::SIMCALORIMETERHIT) || (Collection_in_Physics_Evt->getTypeName() == LCIO::SIMTRACKERHIT))
      {
        try
        {
          define_time_windows(Collection_name);
        }
        catch (std::runtime_error &e)
        {
          streamlog_out(DEBUG) << "Skipping collection without integration times: " << Collection_name << std::endl;
          continue;
        }
        streamlog_out(DEBUG) << "Cropping collection: " << Collection_name << std::endl;
        crop_collection(Collection_in_Physics_Evt);
      }
    }

    // Make sure we have filenames to open and that we really want to overlay something
    if ((_NOverlay > 0.) && (overlay_Evt == nullptr) && (_inputFileNamesMuPlus.size() > 0) && (_inputFileNamesMuMinus.size() > 0))
    {

      std::vector<int> v_file_indices_mupl(_inputFileNamesMuPlus.size());           // vector of indices
      std::iota(std::begin(v_file_indices_mupl), std::end(v_file_indices_mupl), 0); // Fill with 0, 1, ...
      std::shuffle(v_file_indices_mupl.begin(), v_file_indices_mupl.end(), g);

      for (int k = 0; k < _NOverlay; ++k)
      {
        overlay_Eventfile_reader->open(_inputFileNamesMuPlus.at(v_file_indices_mupl[k]));
        streamlog_out(MESSAGE) << "Open mu plus background file: " << _inputFileNamesMuPlus.at(v_file_indices_mupl[k]) << std::endl;

        overlay_Evt = overlay_Eventfile_reader->readNextEvent(LCIO::UPDATE);

        // the overlay_Event is now open, start to merge its collections with the ones of the accumulated overlay events collections
        // all the preparatory work has been done now....
        // first, let's see which collections are in the event

        if (_mergeMCParticles)
        {
          // first include the MCParticles into the physics event
          try
          {
            // Do Not Need DestMap, because this is only MCParticles
            currentDest = _mcParticleCollectionName;
            streamlog_out(DEBUG) << "Merging MCParticles " << std::endl;
            merge_collections(overlay_Evt->getCollection(_mcParticleCollectionName), evt->getCollection(_mcParticleCollectionName), 0.); // last argument is time offset
          }
          catch (DataNotAvailableException &e)
          {
            streamlog_out(ERROR) << "Failed to extract MCParticle collection: " << e.what() << std::endl;
            throw e;
          }
        }

        collection_names_in_Evt = overlay_Evt->getCollectionNames();

        for (unsigned int j = 0; j < collection_names_in_Evt->size(); ++j)
        {
          const std::string Collection_name = collection_names_in_Evt->at(j);

          LCCollection *Collection_in_overlay_Evt = overlay_Evt->getCollection(Collection_name);
          LCCollection *Collection_in_Physics_Evt = 0;

          // Skip the MCParticle collection
          if (Collection_name == _mcParticleCollectionName)
          {
            continue;
          }

          try
          {
            define_time_windows(Collection_name);
          }
          catch (std::runtime_error &e)
          {
            continue;
          }

          // the event can only make contributions to the readout, if the bx does not happen after the integration time stopped.
          // and we are only interested in Calorimeter or Trackerhits.

          if ((this_stop > 0.) &&
              ((Collection_in_overlay_Evt->getTypeName() == LCIO::SIMCALORIMETERHIT) || (Collection_in_overlay_Evt->getTypeName() == LCIO::SIMTRACKERHIT)))
          {
            // Open the same collection in the physics event
            try
            {
              Collection_in_Physics_Evt = evt->getCollection(Collection_name);
            }
            catch (DataNotAvailableException &e)
            {
              streamlog_out(DEBUG) << "Add new Collection" << Collection_in_overlay_Evt->getTypeName() << " with name " << Collection_name << std::endl;
              LCCollectionVec *new_collection = new LCCollectionVec(Collection_in_overlay_Evt->getTypeName());

              StringVec stringKeys;
              Collection_in_overlay_Evt->getParameters().getStringKeys(stringKeys);
              for (unsigned i = 0, nStringKeys = stringKeys.size(); i < nStringKeys; ++i)
              {
                StringVec vals;
                Collection_in_overlay_Evt->getParameters().getStringVals(stringKeys[i], vals);
                new_collection->parameters().setValues(stringKeys[i], vals);
              }
              StringVec intKeys;
              Collection_in_overlay_Evt->getParameters().getIntKeys(intKeys);
              for (unsigned i = 0, nIntKeys = intKeys.size(); i < nIntKeys; ++i)
              {
                IntVec vals;
                Collection_in_overlay_Evt->getParameters().getIntVals(intKeys[i], vals);
                new_collection->parameters().setValues(intKeys[i], vals);
              }
              StringVec floatKeys;
              Collection_in_overlay_Evt->getParameters().getFloatKeys(floatKeys);
              for (unsigned i = 0, nFloatKeys = floatKeys.size(); i < nFloatKeys; ++i)
              {
                FloatVec vals;
                Collection_in_overlay_Evt->getParameters().getFloatVals(floatKeys[i], vals);
                new_collection->parameters().setValues(floatKeys[i], vals);
              }

              evt->addCollection(new_collection, Collection_name);
              Collection_in_Physics_Evt = evt->getCollection(Collection_name);
            }

            // Set DestMap back to the one for the Collection Name...
            currentDest = Collection_name;
            streamlog_out(DEBUG) << "Now overlaying collection " << Collection_name
                                 << " And we have " << collDestMap[currentDest].size() << " Hits in destMap"
                                 << std::endl;
            // Now we merge the collections
            merge_collections(Collection_in_overlay_Evt, Collection_in_Physics_Evt, 0.);
          }
        }
        // overlay_Eventfile_reader->close();
      }

      // Do the same for mu minus
      std::vector<int> v_file_indices_mumi(_inputFileNamesMuMinus.size());          // vector of indices
      std::iota(std::begin(v_file_indices_mumi), std::end(v_file_indices_mumi), 0); // Fill with 0, 1, ...
      std::shuffle(v_file_indices_mumi.begin(), v_file_indices_mumi.end(), g);

      for (int k = 0; k < _NOverlay; ++k)
      {
        overlay_Eventfile_reader->open(_inputFileNamesMuMinus.at(v_file_indices_mumi[k]));
        streamlog_out(MESSAGE) << "Open mu minus background file: " << _inputFileNamesMuMinus.at(v_file_indices_mumi[k]) << std::endl;

        overlay_Evt = overlay_Eventfile_reader->readNextEvent(LCIO::UPDATE);

        // the overlay_Event is now open, start to merge its collections with the ones of the accumulated overlay events collections
        // all the preparatory work has been done now....
        // first, let's see which collections are in the event

        if (_mergeMCParticles)
        {
          // first include the MCParticles into the physics event
          try
          {
            // Do Not Need DestMap, because this is only MCParticles
            currentDest = _mcParticleCollectionName;
            streamlog_out(DEBUG) << "Merging MCParticles " << std::endl;
            merge_collections(overlay_Evt->getCollection(_mcParticleCollectionName), evt->getCollection(_mcParticleCollectionName), 0.); // last argument is time offset
          }
          catch (DataNotAvailableException &e)
          {
            streamlog_out(ERROR) << "Failed to extract MCParticle collection: " << e.what() << std::endl;
            throw e;
          }
        }

        collection_names_in_Evt = overlay_Evt->getCollectionNames();

        for (unsigned int j = 0; j < collection_names_in_Evt->size(); ++j)
        {
          const std::string Collection_name = collection_names_in_Evt->at(j);

          LCCollection *Collection_in_overlay_Evt = overlay_Evt->getCollection(Collection_name);
          LCCollection *Collection_in_Physics_Evt = 0;

          // Skip the MCParticle collection
          if (Collection_name == _mcParticleCollectionName)
          {
            continue;
          }

          try
          {
            define_time_windows(Collection_name);
          }
          catch (std::runtime_error &e)
          {
            continue;
          }

          // the event can only make contributions to the readout, if the bx does not happen after the integration time stopped.
          // and we are only interested in Calorimeter or Trackerhits.

          if ((this_stop > 0.) &&
              ((Collection_in_overlay_Evt->getTypeName() == LCIO::SIMCALORIMETERHIT) || (Collection_in_overlay_Evt->getTypeName() == LCIO::SIMTRACKERHIT)))
          {
            // Open the same collection in the physics event
            Collection_in_Physics_Evt = evt->getCollection(Collection_name);

            // Set DestMap back to the one for the Collection Name...
            currentDest = Collection_name;
            streamlog_out(DEBUG) << "Now overlaying collection " << Collection_name
                                 << " And we have " << collDestMap[currentDest].size() << " Hits in destMap"
                                 << std::endl;
            // Now we merge the collections
            merge_collections(Collection_in_overlay_Evt, Collection_in_Physics_Evt, 0.);
          }
        }
        // overlay_Eventfile_reader->close();
      }
    } // If we have any files, and more than 0 events to overlay end
    else
    {
      streamlog_out(DEBUG) << "Nothing to overlay " << std::endl;
    }

    ++_nEvt;
    // we clear the map of calorimeter hits for the next event
    collDestMap.clear();
    const std::vector<std::string> *collection_names_in_evt = evt->getCollectionNames();

    for (unsigned int i = 0; i < collection_names_in_evt->size(); ++i)
    {
      streamlog_out(DEBUG) << "Collection " << collection_names_in_evt->at(i) << " has now " << evt->getCollection(collection_names_in_evt->at(i))->getNumberOfElements() << " elements" << std::endl;
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void OverlayTimingRandomMix::define_time_windows(std::string const &collectionName)
  {

    this_start = _integrationTimeMin;
    this_stop = 0.0;

    auto iter = _collectionIntegrationTimes.find(collectionName);
    if (iter != _collectionIntegrationTimes.end())
    {
      this_start = iter->second.first;
      this_stop = iter->second.second;
    }
    else
    {
      throw std::runtime_error("Cannot find integration time for collection " + collectionName);
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void OverlayTimingRandomMix::crop_collection(EVENT::LCCollection *collection)
  {
    const int number_of_elements = collection->getNumberOfElements();

    if (number_of_elements > 0)
    {
      if (collection->getTypeName() == LCIO::SIMTRACKERHIT)
      {
        for (int k = number_of_elements - 1; k >= 0; --k)
        {
          SimTrackerHit *TrackerHit = static_cast<SimTrackerHit *>(collection->getElementAt(k));
          const float _time_of_flight = time_of_flight(TrackerHit->getPosition()[0], TrackerHit->getPosition()[1], TrackerHit->getPosition()[2]);
          if (!((TrackerHit->getTime() > (this_start + _time_of_flight)) && (TrackerHit->getTime() < (this_stop + _time_of_flight))))
          {
            collection->removeElementAt(k);
            delete TrackerHit;
          }
        }
      }
      else if (collection->getTypeName() == LCIO::SIMCALORIMETERHIT)
      {
        // we count from top to bottom, in order not to get confused when removing and adding elements!
        for (int i = number_of_elements - 1; i >= 0; --i)
        {
          SimCalorimeterHit *CalorimeterHit = static_cast<SimCalorimeterHit *>(collection->getElementAt(i));
          int not_within_time_window = 0;

          // check whether all entries are within the time window
          const float _time_of_flight = time_of_flight(CalorimeterHit->getPosition()[0], CalorimeterHit->getPosition()[1], CalorimeterHit->getPosition()[2]);

          for (int j = 0; j < CalorimeterHit->getNMCContributions(); ++j)
          {
            // we need to shift the time window to account for the time of flight of the particle...
            if (!((CalorimeterHit->getTimeCont(j) > (this_start + _time_of_flight)) && (CalorimeterHit->getTimeCont(j) < (this_stop + _time_of_flight))))
            {
              ++not_within_time_window;
              // std::cout << " calo hit : " << j << " Time : " << CalorimeterHit->getTimeCont(j) << " ?> " << this_start + _time_of_flight << " ?< " << this_stop + _time_of_flight << std::endl;
            }
          }

          // if one and not all MC contribution is not within the time window....
          if (not_within_time_window == 0)
          {
            collDestMap[currentDest].insert(DestMap::value_type(cellID2long(CalorimeterHit->getCellID0(), CalorimeterHit->getCellID1()), CalorimeterHit));
          }
          else if ((not_within_time_window > 0) && (not_within_time_window < CalorimeterHit->getNMCContributions()))
          {
            SimCalorimeterHitImpl *newCalorimeterHit = new SimCalorimeterHitImpl();

            for (int j = 0; j < CalorimeterHit->getNMCContributions(); ++j)
            {
              if ((CalorimeterHit->getTimeCont(j) > (this_start + _time_of_flight)) && (CalorimeterHit->getTimeCont(j) < (this_stop + _time_of_flight)))
              {
                newCalorimeterHit->addMCParticleContribution(CalorimeterHit->getParticleCont(j), CalorimeterHit->getEnergyCont(j), CalorimeterHit->getTimeCont(j));
              }
            }

            newCalorimeterHit->setCellID0(CalorimeterHit->getCellID0());
            newCalorimeterHit->setCellID1(CalorimeterHit->getCellID1());
            float ort[3] = {CalorimeterHit->getPosition()[0], CalorimeterHit->getPosition()[1], CalorimeterHit->getPosition()[2]};
            newCalorimeterHit->setPosition(ort);

            collection->removeElementAt(i);
            delete CalorimeterHit;

            collection->addElement(newCalorimeterHit);
            collDestMap[currentDest].insert(DestMap::value_type(cellID2long(newCalorimeterHit->getCellID0(), newCalorimeterHit->getCellID1()), newCalorimeterHit));
          }
          else if (not_within_time_window == CalorimeterHit->getNMCContributions())
          {
            collection->removeElementAt(i);
            delete CalorimeterHit;
          }
        }
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void OverlayTimingRandomMix::merge_collections(EVENT::LCCollection *source_collection, EVENT::LCCollection *dest_collection, float time_offset)
  {
    // first, calculate the integration time, depending on the subdetector
    // time offset is the time of the physics event, after the start of the bunch train
    // adding the time offset shall move the background event relative to the physics event...
    const int number_of_elements = source_collection->getNumberOfElements();
    int mergedN = 0;
    streamlog_out(DEBUG) << "We are starting the merge with " << dest_collection->getNumberOfElements() << std::endl;
    if (number_of_elements > 0)
    {
      if (source_collection->getTypeName() == LCIO::MCPARTICLE)
      {
        for (int i = number_of_elements - 1; i >= 0; --i)
        {
          MCParticleImpl *MC_Part = static_cast<MCParticleImpl *>(source_collection->getElementAt(i));
          MC_Part->setTime(MC_Part->getTime() + time_offset);
          dest_collection->addElement(MC_Part);
          source_collection->removeElementAt(i);
        }
      }
      else if (source_collection->getTypeName() == LCIO::SIMTRACKERHIT)
      {
        // If truth is removed in the overlay, we need to remove the pointer since it will become invalid
        if (not _mergeMCParticles)
        {
          for (int k = 0; k < number_of_elements; ++k)
          {
            SimTrackerHitImpl *TrackerHit = static_cast<SimTrackerHitImpl *>(source_collection->getElementAt(k));
            EVENT::MCParticle *mcParticle = TrackerHit->getMCParticle();
            // if a valida pointer exists, keep momentum information
            if (mcParticle)
            {
              const double *preserveMomentum = mcParticle->getMomentum();
              TrackerHit->setMomentum(preserveMomentum[0], preserveMomentum[1], preserveMomentum[2]);
            }
            TrackerHit->setMCParticle(nullptr);
          }
        }
        // Adjust timing information on hits depending on BX and set overlay flag
        for (int k = number_of_elements - 1; k >= 0; --k)
        {
          SimTrackerHitImpl *TrackerHit = static_cast<SimTrackerHitImpl *>(source_collection->getElementAt(k));

          TrackerHit->setOverlay(true);

          const float _time_of_flight = time_of_flight(TrackerHit->getPosition()[0], TrackerHit->getPosition()[1], TrackerHit->getPosition()[2]);

          if (((TrackerHit->getTime() + time_offset) > (this_start + _time_of_flight)) && ((TrackerHit->getTime() + time_offset) < (this_stop + _time_of_flight)))
          {
            TrackerHit->setTime(TrackerHit->getTime() + time_offset);
            dest_collection->addElement(TrackerHit);
            source_collection->removeElementAt(k);
          }
        }
      }
      else if (source_collection->getTypeName() == LCIO::SIMCALORIMETERHIT)
      {
        // create a map of dest Collection
        for (int k = number_of_elements - 1; k >= 0; --k)
        {
          SimCalorimeterHit *CalorimeterHit = static_cast<SimCalorimeterHit *>(source_collection->getElementAt(k));
          const float _time_of_flight = time_of_flight(CalorimeterHit->getPosition()[0], CalorimeterHit->getPosition()[1], CalorimeterHit->getPosition()[2]);

          // check whether there is already a hit at this position
          const unsigned long long lookfor = cellID2long(CalorimeterHit->getCellID0(), CalorimeterHit->getCellID1());
          DestMap::const_iterator destMapIt = collDestMap[currentDest].find(lookfor);
          if (destMapIt == collDestMap[currentDest].end())
          {
            // There is no Hit at this position -- the new hit can be added, if it is not outside the window
            SimCalorimeterHitImpl *newCalorimeterHit = new SimCalorimeterHitImpl();
            bool add_Hit = false;

            for (int j = 0; j < CalorimeterHit->getNMCContributions(); ++j)
            {
              if (((CalorimeterHit->getTimeCont(j) + time_offset) > (this_start + _time_of_flight)) && ((CalorimeterHit->getTimeCont(j) + time_offset) < (this_stop + _time_of_flight)))
              {
                add_Hit = true;
                newCalorimeterHit->addMCParticleContribution(CalorimeterHit->getParticleCont(j), CalorimeterHit->getEnergyCont(j), CalorimeterHit->getTimeCont(j) + time_offset);
              }
            }
            if (add_Hit)
            {
              newCalorimeterHit->setCellID0(CalorimeterHit->getCellID0());
              newCalorimeterHit->setCellID1(CalorimeterHit->getCellID1());
              float ort[3] = {CalorimeterHit->getPosition()[0], CalorimeterHit->getPosition()[1], CalorimeterHit->getPosition()[2]};
              newCalorimeterHit->setPosition(ort);
              dest_collection->addElement(newCalorimeterHit);
              collDestMap[currentDest].insert(DestMap::value_type(cellID2long(newCalorimeterHit->getCellID0(), newCalorimeterHit->getCellID1()), newCalorimeterHit));
            }
            else
            {
              delete newCalorimeterHit;
            }
          }
          else
          {
            // there is already a hit at this position....
            SimCalorimeterHitImpl *newCalorimeterHit = static_cast<SimCalorimeterHitImpl *>(destMapIt->second);
            ++mergedN;
            if ((newCalorimeterHit->getPosition()[0] - CalorimeterHit->getPosition()[0]) *
                        (newCalorimeterHit->getPosition()[0] - CalorimeterHit->getPosition()[0]) +
                    (newCalorimeterHit->getPosition()[1] - CalorimeterHit->getPosition()[1]) *
                        (newCalorimeterHit->getPosition()[1] - CalorimeterHit->getPosition()[1]) +
                    (newCalorimeterHit->getPosition()[2] - CalorimeterHit->getPosition()[2]) *
                        (newCalorimeterHit->getPosition()[2] - CalorimeterHit->getPosition()[2]) >
                10)
            {
              streamlog_out(ERROR) << "HITS DO NOT MATCH in " << currentDest << "!!!" << std::endl;
              streamlog_out(ERROR) << "X New  " << newCalorimeterHit->getPosition()[0]
                                   << "  Old  " << CalorimeterHit->getPosition()[0] << std::endl;
              streamlog_out(ERROR) << "Y New  " << newCalorimeterHit->getPosition()[1]
                                   << "  Old  " << CalorimeterHit->getPosition()[1] << std::endl;
              streamlog_out(ERROR) << "Z New  " << newCalorimeterHit->getPosition()[2]
                                   << "  Old  " << CalorimeterHit->getPosition()[2] << std::endl;
              streamlog_out(ERROR) << "ID0New  " << newCalorimeterHit->getCellID0()
                                   << "   Old  " << CalorimeterHit->getCellID0() << std::endl;
              streamlog_out(ERROR) << "ID1New  " << newCalorimeterHit->getCellID1()
                                   << "   Old  " << CalorimeterHit->getCellID1() << std::endl;

              //                      std::exit(1);
              streamlog_out(ERROR) << "ID1New  " << cellID2long(newCalorimeterHit->getCellID0(), newCalorimeterHit->getCellID1())
                                   << "   Old  " << cellID2long(CalorimeterHit->getCellID0(), CalorimeterHit->getCellID1())
                                   << std::endl;
            }
            for (int j = 0; j < CalorimeterHit->getNMCContributions(); ++j)
            {
              if (((CalorimeterHit->getTimeCont(j) + time_offset) > (this_start + _time_of_flight)) && ((CalorimeterHit->getTimeCont(j) + time_offset) < (this_stop + _time_of_flight)))
              {
                newCalorimeterHit->addMCParticleContribution(CalorimeterHit->getParticleCont(j), CalorimeterHit->getEnergyCont(j), CalorimeterHit->getTimeCont(j) + time_offset);
              }
            }
          }
        }
      }
    }
    streamlog_out(DEBUG) << "We are ending the merge with " << dest_collection->getNumberOfElements()
                         << " and we merged " << mergedN << "  others  "
                         << std::endl;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void OverlayTimingRandomMix::check(EVENT::LCEvent *)
  {
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void OverlayTimingRandomMix::end()
  {
    delete overlay_Eventfile_reader;
    overlay_Eventfile_reader = nullptr;
  }

} // namespace
