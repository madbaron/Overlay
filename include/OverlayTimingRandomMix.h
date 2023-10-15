#ifndef OverlayTimingRandomMix_h
#define OverlayTimingRandomMix_h 1

#include "marlin/Processor.h"
#include "marlin/EventModifier.h"

#include "lcio.h"

#include <cmath>
#include <limits>

namespace EVENT
{
  class SimCalorimeterHit;
  class LCRunHeader;
  class LCEvent;
  class LCCollection;
}

namespace overlay
{

  /** OverlayTimingRandomMix processor for overlaying background to each bunch crossing of a bunch train.
   *
   *  A physics event is placed at a random or fixed position of the bunch train.
   *  Then background events are overlayed to each bunch crossing of the train and merged into the physics event.
   *
   *  For the merging, integration times can be given for each collection or subdetector, respectively.
   *  Then only hits are added to the physics event's collections, which fall into the integration time window of
   *  the specific subdetector
   *
   *  @author F. Meloni DESY (based on Overlay processor by F. Gaede)
   *
   *  @param BackgroundFileNames (StringVec) The names (with absolute or relative pathes) of the files from
   *  which the background should be read.
   *  It is the users responsibility to provide sufficient statistics for the signal
   *  sample under study.
   *
   *  @param NumberBackground - default 1 -- Number of background events to overlay to one BX. This is either fixed (if Poisson_random_NOverlay = false) of the mean value if  (if Poisson_random_NOverlay = true)
   *
   *  @param Poisson_random_NOverlay - default false -- Overlay a random number of background events to each BX. If t, the number is pulled from a Poisson distribution.
   *
   *  @param RandomSeed (int) random seed - default 42
   *
   */
  class OverlayTimingRandomMix : public marlin::Processor, public marlin::EventModifier
  {
  public:
    virtual marlin::Processor *newProcessor();

    OverlayTimingRandomMix();
    OverlayTimingRandomMix(std::string const &name);
    OverlayTimingRandomMix(OverlayTimingRandomMix const &) = delete;
    OverlayTimingRandomMix &operator=(OverlayTimingRandomMix const &) = delete;

    virtual void init();

    virtual const std::string &name() const;

    virtual void processRunHeader(EVENT::LCRunHeader *run);

    virtual void modifyEvent(EVENT::LCEvent *evt);

    virtual void check(EVENT::LCEvent *evt);

    virtual void end();

  protected:
    float time_of_flight(float x, float y, float z) const;

    virtual void define_time_windows(const std::string &collectionName);

    void crop_collection(EVENT::LCCollection *collection);

    void merge_collections(EVENT::LCCollection *source_collection, EVENT::LCCollection *dest_collection, float time_offset);

    unsigned long long cellID2long(unsigned int id0, unsigned int id1) const;

    unsigned int _nRun = 0;
    unsigned int _nEvt = 0;
    StringVec _inputFileNamesMuPlus{};
    StringVec _inputFileNamesMuMinus{};

    std::vector<std::string> _collectionTimesVec{"BeamCalCollection", "10"};
    std::map<std::string, std::pair<float, float>> _collectionIntegrationTimes{};

    float _NOverlay = 1;
    float _integrationTimeMin = -0.5;

    IO::LCReader *overlay_Eventfile_reader = NULL;
    LCEvent *overlay_Evt = nullptr;
    int m_eventCounter = 0;
    int m_currentFileIndex = 0;

    float this_start = -0.25;
    float this_stop = std::numeric_limits<float>::max();
    bool TPC_hits = false;

    bool _mergeMCParticles = true;
    std::string _mcParticleCollectionName = "";
    std::string currentDest = "";

    typedef std::map<unsigned long long, EVENT::SimCalorimeterHit *> DestMap;
    typedef std::map<std::string, DestMap> CollDestMap;
    CollDestMap collDestMap{};
  };

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline marlin::Processor *OverlayTimingRandomMix::newProcessor()
  {
    return new OverlayTimingRandomMix;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline const std::string &OverlayTimingRandomMix::name() const
  {
    return marlin::Processor::name();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline float OverlayTimingRandomMix::time_of_flight(float x, float y, float z) const
  {
    // returns the time of flight to the radius in ns
    //  mm/m/s = 10^{-3}s = 10^6 ns d.h. 299 mm/ns
    return std::sqrt((x * x) + (y * y) + (z * z)) / 299.792458;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline unsigned long long OverlayTimingRandomMix::cellID2long(unsigned int id0, unsigned int id1) const
  {
    const unsigned long long newID = ((unsigned long long)(id0) << sizeof(unsigned int) * 8 | (unsigned long long)id1);
    return newID;
  }

} // namespace

#endif
