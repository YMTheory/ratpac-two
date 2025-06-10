#include <RAT/DsG4Cerenkov.hh>

#include "G4LossTableManager.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleMomentum.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsModelCatalog.hh"
#include "G4Poisson.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DsG4Cerenkov::DsG4Cerenkov(const G4String& processName, G4ProcessType type)
    : G4VProcess(processName, type), fNumPhotons(0) {
  secID = G4PhysicsModelCatalog::GetModelID("model_Cerenkov");
  SetProcessSubType(fCerenkov);

  thePhysicsTable = nullptr;

  if (verboseLevel > 0) {
    G4cout << GetProcessName() << " is created." << G4endl;
  }
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DsG4Cerenkov::~DsG4Cerenkov() {
  if (thePhysicsTable != nullptr) {
    thePhysicsTable->clearAndDestroy();
    delete thePhysicsTable;
  }
}

void DsG4Cerenkov::ProcessDescription(std::ostream& out) const {
  out << "The Cerenkov effect simulates optical photons created by the\n";
  out << "passage of charged particles through matter. Materials need\n";
  out << "to have the property RINDEX (refractive index) defined.\n";
  G4VProcess::DumpInfo();

  G4OpticalParameters* params = G4OpticalParameters::Instance();
  out << "Maximum beta change per step: " << params->GetCerenkovMaxBetaChange();
  out << "Maximum photons per step: " << params->GetCerenkovMaxPhotonsPerStep();
  out << "Track secondaries first: " << params->GetCerenkovTrackSecondariesFirst();
  out << "Stack photons: " << params->GetCerenkovStackPhotons();
  out << "Verbose level: " << params->GetCerenkovVerboseLevel();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool DsG4Cerenkov::IsApplicable(const G4ParticleDefinition& aParticleType) {
  return (aParticleType.GetPDGCharge() != 0.0 && aParticleType.GetPDGMass() != 0.0 &&
          aParticleType.GetParticleName() != "chargedgeantino" && !aParticleType.IsShortLived())
             ? true
             : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DsG4Cerenkov::Initialise() {
  G4OpticalParameters* params = G4OpticalParameters::Instance();
  SetMaxBetaChangePerStep(params->GetCerenkovMaxBetaChange());
  SetMaxNumPhotonsPerStep(params->GetCerenkovMaxPhotonsPerStep());
  SetTrackSecondariesFirst(params->GetCerenkovTrackSecondariesFirst());
  SetStackPhotons(params->GetCerenkovStackPhotons());
  SetVerboseLevel(params->GetCerenkovVerboseLevel());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DsG4Cerenkov::BuildPhysicsTable(const G4ParticleDefinition&) {
  if (thePhysicsTable) return;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  std::size_t numOfMaterials = G4Material::GetNumberOfMaterials();

  thePhysicsTable = new G4PhysicsTable(numOfMaterials);

  // loop over materials
  for (std::size_t i = 0; i < numOfMaterials; ++i) {
    G4PhysicsFreeVector* cerenkovIntegral = nullptr;

    // Retrieve vector of refraction indices for the material
    // from the material's optical properties table
    G4Material* aMaterial = (*theMaterialTable)[i];
    G4MaterialPropertiesTable* MPT = aMaterial->GetMaterialPropertiesTable();

    if (MPT) {
      cerenkovIntegral = new G4PhysicsFreeVector();
      G4MaterialPropertyVector* refractiveIndex = MPT->GetProperty(kRINDEX);

      if (refractiveIndex) {
        // Retrieve the first refraction index in vector
        // of (photon energy, refraction index) pairs
        G4double currentRI = (*refractiveIndex)[0];
        if (currentRI > 1.0) {
          // Create first (photon energy, Cerenkov Integral) pair
          G4double currentPM = refractiveIndex->Energy(0);
          G4double currentCAI = 0.0;

          cerenkovIntegral->InsertValues(currentPM, currentCAI);

          // Set previous values to current ones prior to loop
          G4double prevPM = currentPM;
          G4double prevCAI = currentCAI;
          G4double prevRI = currentRI;

          // loop over all (photon energy, refraction index)
          // pairs stored for this material
          for (std::size_t ii = 1; ii < refractiveIndex->GetVectorLength(); ++ii) {
            currentRI = (*refractiveIndex)[ii];
            currentPM = refractiveIndex->Energy(ii);
            currentCAI =
                prevCAI + (currentPM - prevPM) * 0.5 * (1.0 / (prevRI * prevRI) + 1.0 / (currentRI * currentRI));

            cerenkovIntegral->InsertValues(currentPM, currentCAI);

            prevPM = currentPM;
            prevCAI = currentCAI;
            prevRI = currentRI;
          }
        }
      } else {
      }
    }

    // The Cerenkov integral for a given material will be inserted in
    // thePhysicsTable according to the position of the material in
    // the material table.
    thePhysicsTable->insertAt(i, cerenkovIntegral);
  }
}

G4VParticleChange* DsG4Cerenkov::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)

{
  G4cout << "DsG4Cerenkov Post Step Do It run here" << G4endl;
  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4Material* aMaterial = aTrack.GetMaterial();

  G4StepPoint* pPreStepPoint = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

  G4ThreeVector x0 = pPreStepPoint->GetPosition();
  G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
  G4double t0 = pPreStepPoint->GetGlobalTime();

  G4MaterialPropertiesTable* MPT = aMaterial->GetMaterialPropertiesTable();
  if (!MPT) return pParticleChange;

  G4MaterialPropertyVector* Rindex = MPT->GetProperty(kRINDEX);
  if (!Rindex) return pParticleChange;

  G4double charge = aParticle->GetDefinition()->GetPDGCharge();
  G4double beta = (pPreStepPoint->GetBeta() + pPostStepPoint->GetBeta()) * 0.5;

  G4double MeanNumberOfPhotons = GetAverageNumberOfPhotons(charge, beta, aMaterial, Rindex);

  if (MeanNumberOfPhotons <= 0.0) {
    // return unchanged particle and no secondaries
    aParticleChange.SetNumberOfSecondaries(0);
    return pParticleChange;
  }

  G4double step_length = aStep.GetStepLength();
  MeanNumberOfPhotons = MeanNumberOfPhotons * step_length;
  fNumPhotons = (G4int)G4Poisson(MeanNumberOfPhotons);

  if (fNumPhotons <= 0 || !fStackingFlag) {
    // return unchanged particle and no secondaries
    aParticleChange.SetNumberOfSecondaries(0);
    return pParticleChange;
  }
  ////////////////////////////////////////////////////////////////
  aParticleChange.SetNumberOfSecondaries(fNumPhotons);

  if (fTrackSecondariesFirst) {
    if (aTrack.GetTrackStatus() == fAlive) aParticleChange.ProposeTrackStatus(fSuspend);
  }

  ////////////////////////////////////////////////////////////////
  G4double Pmin = Rindex->Energy(0);
  G4double Pmax = Rindex->GetMaxEnergy();
  G4double dp = Pmax - Pmin;

  G4double nMax = Rindex->GetMaxValue();
  G4double BetaInverse = 1. / beta;

  G4double maxCos = BetaInverse / nMax;
  G4double maxSin2 = (1.0 - maxCos) * (1.0 + maxCos);

  G4double beta1 = pPreStepPoint->GetBeta();
  G4double beta2 = pPostStepPoint->GetBeta();

  G4double MeanNumberOfPhotons1 = GetAverageNumberOfPhotons(charge, beta1, aMaterial, Rindex);
  G4double MeanNumberOfPhotons2 = GetAverageNumberOfPhotons(charge, beta2, aMaterial, Rindex);

  for (G4int i = 0; i < fNumPhotons; ++i) {
    // Determine photon energy
    G4double rand;
    G4double sampledEnergy, sampledRI;
    G4double cosTheta, sin2Theta;

    // sample an energy
    do {
      rand = G4UniformRand();
      sampledEnergy = Pmin + rand * dp;
      sampledRI = Rindex->Value(sampledEnergy);
      cosTheta = BetaInverse / sampledRI;

      sin2Theta = (1.0 - cosTheta) * (1.0 + cosTheta);
      rand = G4UniformRand();

      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    } while (rand * maxSin2 > sin2Theta);

    // Create photon momentum direction vector. The momentum direction is still
    // with respect to the coordinate system where the primary particle
    // direction is aligned with the z axis
    rand = G4UniformRand();
    G4double phi = twopi * rand;
    G4double sinPhi = std::sin(phi);
    G4double cosPhi = std::cos(phi);
    G4double sinTheta = std::sqrt(sin2Theta);
    G4ParticleMomentum photonMomentum(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);

    // Rotate momentum direction back to global reference system
    photonMomentum.rotateUz(p0);

    // Determine polarization of new photon
    G4ThreeVector photonPolarization(cosTheta * cosPhi, cosTheta * sinPhi, -sinTheta);

    // Rotate back to original coord system
    photonPolarization.rotateUz(p0);

    // Generate a new photon:
    auto aCerenkovPhoton = new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), photonMomentum);

    aCerenkovPhoton->SetPolarization(photonPolarization);
    aCerenkovPhoton->SetKineticEnergy(sampledEnergy);

    G4double NumberOfPhotons, N;

    do {
      rand = G4UniformRand();
      NumberOfPhotons = MeanNumberOfPhotons1 - rand * (MeanNumberOfPhotons1 - MeanNumberOfPhotons2);
      N = G4UniformRand() * std::max(MeanNumberOfPhotons1, MeanNumberOfPhotons2);
      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    } while (N > NumberOfPhotons);

    G4double delta = rand * aStep.GetStepLength();
    G4double deltaTime = delta / (pPreStepPoint->GetVelocity() +
                                  rand * (pPostStepPoint->GetVelocity() - pPreStepPoint->GetVelocity()) * 0.5);

    G4double aSecondaryTime = t0 + deltaTime;
    G4ThreeVector aSecondaryPosition = x0 + rand * aStep.GetDeltaPosition();

    // Generate new G4Track object:
    G4Track* aSecondaryTrack = new G4Track(aCerenkovPhoton, aSecondaryTime, aSecondaryPosition);

    aSecondaryTrack->SetTouchableHandle(aStep.GetPreStepPoint()->GetTouchableHandle());
    aSecondaryTrack->SetParentID(aTrack.GetTrackID());
    aSecondaryTrack->SetCreatorModelID(secID);
    aParticleChange.AddSecondary(aSecondaryTrack);
  }

  if (verboseLevel > 1) {
    G4cout << "\n Exiting from DsG4Cerenkov::DoIt -- NumberOfSecondaries = " << aParticleChange.GetNumberOfSecondaries()
           << G4endl;
  }

  return pParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DsG4Cerenkov::PreparePhysicsTable(const G4ParticleDefinition&) { Initialise(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DsG4Cerenkov::GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*) { return 1.; }

G4double DsG4Cerenkov::PostStepGetPhysicalInteractionLength(const G4Track& aTrack, G4double,
                                                            G4ForceCondition* condition) {
  *condition = NotForced;
  G4double StepLimit = DBL_MAX;
  fNumPhotons = 0;

  const G4Material* aMaterial = aTrack.GetMaterial();
  std::size_t materialIndex = aMaterial->GetIndex();

  // If Physics Vector is not defined no Cerenkov photons
  if (!(*thePhysicsTable)[materialIndex]) {
    return StepLimit;
  }

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();

  G4double kineticEnergy = aParticle->GetKineticEnergy();
  const G4ParticleDefinition* particleType = aParticle->GetDefinition();
  G4double mass = particleType->GetPDGMass();

  G4double beta = aParticle->GetTotalMomentum() / aParticle->GetTotalEnergy();
  G4double gamma = aParticle->GetTotalEnergy() / mass;

  G4MaterialPropertiesTable* aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();

  G4MaterialPropertyVector* Rindex = nullptr;
  if (aMaterialPropertiesTable) Rindex = aMaterialPropertiesTable->GetProperty(kRINDEX);

  G4double nMax;
  if (Rindex) {
    nMax = Rindex->GetMaxValue();
  } else {
    return StepLimit;
  }

  G4double BetaMin = 1. / nMax;
  if (BetaMin >= 1.) return StepLimit;

  G4double GammaMin = 1. / std::sqrt(1. - BetaMin * BetaMin);
  if (gamma < GammaMin) return StepLimit;

  G4double kinEmin = mass * (GammaMin - 1.);
  G4double RangeMin = G4LossTableManager::Instance()->GetRange(particleType, kinEmin, couple);
  G4double Range = G4LossTableManager::Instance()->GetRange(particleType, kineticEnergy, couple);
  G4double Step = Range - RangeMin;

  // If the step is smaller than G4ThreeVector::getTolerance(), it may happen
  // that the particle does not move. See bug 1992.
  static const G4double minAllowedStep = G4ThreeVector::getTolerance();
  if (Step < minAllowedStep) return StepLimit;

  if (Step < StepLimit) StepLimit = Step;

  // If user has defined an average maximum number of photons to be generated in
  // a Step, then calculate the Step length for that number of photons.
  if (fMaxPhotons > 0) {
    const G4double charge = aParticle->GetDefinition()->GetPDGCharge();
    G4double MeanNumberOfPhotons = GetAverageNumberOfPhotons(charge, beta, aMaterial, Rindex);
    Step = 0.;
    if (MeanNumberOfPhotons > 0.0) Step = fMaxPhotons / MeanNumberOfPhotons;
    if (Step > 0. && Step < StepLimit) StepLimit = Step;
  }

  // If user has defined an maximum allowed change in beta per step
  if (fMaxBetaChange > 0.) {
    G4double dedx = G4LossTableManager::Instance()->GetDEDX(particleType, kineticEnergy, couple);
    G4double deltaGamma = gamma - 1. / std::sqrt(1. - beta * beta * (1. - fMaxBetaChange) * (1. - fMaxBetaChange));

    Step = mass * deltaGamma / dedx;
    if (Step > 0. && Step < StepLimit) StepLimit = Step;
  }

  *condition = StronglyForced;
  return StepLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DsG4Cerenkov::GetAverageNumberOfPhotons(const G4double charge, const G4double beta,
                                                 const G4Material* aMaterial, G4MaterialPropertyVector* Rindex) const
// This routine computes the number of Cerenkov photons produced per
// Geant4-unit (millimeter) in the current medium.
{
  constexpr G4double Rfact = 369.81 / (eV * cm);
  if (beta <= 0.0) return 0.0;
  G4double BetaInverse = 1. / beta;

  // Vectors used in computation of Cerenkov Angle Integral:
  //    - Refraction Indices for the current material
  //    - new G4PhysicsFreeVector allocated to hold CAI's
  std::size_t materialIndex = aMaterial->GetIndex();

  // Retrieve the Cerenkov Angle Integrals for this material
  G4PhysicsVector* CerenkovAngleIntegrals = ((*thePhysicsTable)(materialIndex));

  std::size_t length = CerenkovAngleIntegrals->GetVectorLength();
  if (0 == length) return 0.0;

  // Min and Max photon energies
  G4double Pmin = Rindex->Energy(0);
  G4double Pmax = Rindex->GetMaxEnergy();

  // Min and Max Refraction Indices
  G4double nMin = Rindex->GetMinValue();
  G4double nMax = Rindex->GetMaxValue();

  // Max Cerenkov Angle Integral
  G4double CAImax = (*CerenkovAngleIntegrals)[length - 1];

  G4double dp, ge;
  // If n(Pmax) < 1/Beta -- no photons generated
  if (nMax < BetaInverse) {
    dp = 0.0;
    ge = 0.0;
  }
  // otherwise if n(Pmin) >= 1/Beta -- photons generated
  else if (nMin > BetaInverse) {
    dp = Pmax - Pmin;
    ge = CAImax;
  }
  // If n(Pmin) < 1/Beta, and n(Pmax) >= 1/Beta, then we need to find a P such
  // that the value of n(P) == 1/Beta. Interpolation is performed by the
  // GetEnergy() and Value() methods of the G4MaterialPropertiesTable and
  // the Value() method of G4PhysicsVector.
  else {
    Pmin = Rindex->GetEnergy(BetaInverse);
    dp = Pmax - Pmin;

    G4double CAImin = CerenkovAngleIntegrals->Value(Pmin);
    ge = CAImax - CAImin;

    if (verboseLevel > 1) {
      G4cout << "CAImin = " << CAImin << G4endl << "ge = " << ge << G4endl;
    }
  }

  // Calculate number of photons
  G4double NumPhotons = Rfact * charge / eplus * charge / eplus * (dp - ge * BetaInverse * BetaInverse);

  return NumPhotons;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DsG4Cerenkov::SetTrackSecondariesFirst(const G4bool state) {
  fTrackSecondariesFirst = state;
  G4OpticalParameters::Instance()->SetCerenkovTrackSecondariesFirst(fTrackSecondariesFirst);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DsG4Cerenkov::SetMaxBetaChangePerStep(const G4double value) {
  fMaxBetaChange = value * CLHEP::perCent;
  G4OpticalParameters::Instance()->SetCerenkovMaxBetaChange(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DsG4Cerenkov::SetMaxNumPhotonsPerStep(const G4int NumPhotons) {
  fMaxPhotons = NumPhotons;
  G4OpticalParameters::Instance()->SetCerenkovMaxPhotonsPerStep(fMaxPhotons);
}

void DsG4Cerenkov::SetStackPhotons(const G4bool stackingFlag) {
  fStackingFlag = stackingFlag;
  G4OpticalParameters::Instance()->SetCerenkovStackPhotons(fStackingFlag);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DsG4Cerenkov::DumpPhysicsTable() const {
  G4cout << "Dump Physics Table!" << G4endl;
  for (std::size_t i = 0; i < thePhysicsTable->entries(); ++i) {
    (*thePhysicsTable)[i]->DumpValues();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DsG4Cerenkov::SetVerboseLevel(G4int verbose) {
  verboseLevel = verbose;
  G4OpticalParameters::Instance()->SetCerenkovVerboseLevel(verboseLevel);
}
