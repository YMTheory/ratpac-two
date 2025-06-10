#ifndef DsG4Cerenkov_h
#define DsG4Cerenkov_h 1

#include "G4DynamicParticle.hh"
#include "G4ForceCondition.hh"
#include "G4GPILSelection.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4VProcess.hh"
#include "globals.hh"

class G4Material;
class G4ParticleDefinition;
class G4PhysicsTable;
class G4Step;
class G4Track;
class G4VParticleChange;

class DsG4Cerenkov : public G4VProcess

{
 public:
  explicit DsG4Cerenkov(const G4String& processName = "DsG4Cerenkov", G4ProcessType type = fElectromagnetic);
  ~DsG4Cerenkov();

  explicit DsG4Cerenkov(const DsG4Cerenkov& right);

  DsG4Cerenkov& operator=(const DsG4Cerenkov& right) = delete;

  G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override;

  void BuildPhysicsTable(const G4ParticleDefinition& aParticleType) override;

  void PreparePhysicsTable(const G4ParticleDefinition& part) override;
  void Initialise();

  G4double GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition*);

  G4double PostStepGetPhysicalInteractionLength(const G4Track& aTrack, G4double, G4ForceCondition*) override;

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) override;

  virtual G4double AlongStepGetPhysicalInteractionLength(const G4Track&, G4double, G4double, G4double&,
                                                         G4GPILSelection*) override {
    return -1.0;
  };

  virtual G4double AtRestGetPhysicalInteractionLength(const G4Track&, G4ForceCondition*) override { return -1.0; };

  virtual G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) override { return nullptr; };

  virtual G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) override { return nullptr; };

  void SetTrackSecondariesFirst(const G4bool state);

  G4bool GetTrackSecondariesFirst() const;

  void SetMaxBetaChangePerStep(const G4double d);

  G4double GetMaxBetaChangePerStep() const;

  void SetMaxNumPhotonsPerStep(const G4int NumPhotons);

  G4int GetMaxNumPhotonsPerStep() const;

  void SetStackPhotons(const G4bool);

  G4bool GetStackPhotons() const;

  G4int GetNumPhotons() const;

  G4PhysicsTable* GetPhysicsTable() const;

  void DumpPhysicsTable() const;

  G4double GetAverageNumberOfPhotons(const G4double charge, const G4double beta, const G4Material* aMaterial,
                                     G4MaterialPropertyVector* Rindex) const;

  void DumpInfo() const override { ProcessDescription(G4cout); };
  void ProcessDescription(std::ostream& out) const override;

  void SetVerboseLevel(G4int);
  // sets verbosity

 protected:
  G4PhysicsTable* thePhysicsTable;

 private:
  G4bool fTrackSecondariesFirst;
  G4double fMaxBetaChange;
  G4int fMaxPhotons;

  G4bool fStackingFlag;

  G4int fNumPhotons;

  G4int secID = -1;
};

inline G4bool DsG4Cerenkov::GetTrackSecondariesFirst() const { return fTrackSecondariesFirst; }

inline G4double DsG4Cerenkov::GetMaxBetaChangePerStep() const { return fMaxBetaChange; }

inline G4int DsG4Cerenkov::GetMaxNumPhotonsPerStep() const { return fMaxPhotons; }

inline G4bool DsG4Cerenkov::GetStackPhotons() const { return fStackingFlag; }

inline G4int DsG4Cerenkov::GetNumPhotons() const { return fNumPhotons; }

inline G4PhysicsTable* DsG4Cerenkov::GetPhysicsTable() const { return thePhysicsTable; }

#endif
