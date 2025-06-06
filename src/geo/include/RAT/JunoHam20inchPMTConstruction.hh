#ifndef __RAT_JunoHam20inchPMTConstruction__
#define __RAT_JunoHam20inchPMTConstruction__

#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4OpticalSurface.hh>
#include <G4PVPlacement.hh>
#include <G4VSensitiveDetector.hh>
#include <G4VSolid.hh>
#include <RAT/DB.hh>
#include <RAT/PMTConstruction.hh>
#include <string>
#include <vector>

namespace RAT {

struct JunoHam20inchPMTConstructionParams {
  JunoHam20inchPMTConstructionParams() {
    invisible = false;
    efficiencyCorrection = 1.0;  // default to 1.0 for no correction
    photocathode_MINrho = 0.0;
    photocathode_MAXrho = 0.0;
  };

  bool invisible;

  double dynodeRadius;         // mm
  double dynodeTop;            // mm
  double dynodeHeight;         // mm
  double photocathode_MINrho;  // mm
  double photocathode_MAXrho;  // mm

  std::vector<double> rInner, zInner, rEdge, zEdge;

  G4Material *exterior;
  G4Material *glass;
  G4Material *vacuum;
  G4Material *dynode;

  G4OpticalSurface *photocathode;
  G4OpticalSurface *mirror;
  G4OpticalSurface *dynode_surface;

  double efficiencyCorrection;
};

// Construction for PMTs based on G4Polycon
class JunoHam20inchPMTConstruction : public PMTConstruction {
 public:
  JunoHam20inchPMTConstruction(DBLinkPtr params, G4LogicalVolume *mother);
  virtual ~JunoHam20inchPMTConstruction() {}

  virtual G4LogicalVolume *BuildVolume(const std::string &prefix);
  virtual G4VSolid *BuildSolid(const std::string &prefix);
  virtual G4PVPlacement *PlacePMT(G4RotationMatrix *pmtrot, G4ThreeVector pmtpos, const std::string &name,
                                  G4LogicalVolume *logi_pmt, G4VPhysicalVolume *mother_phys, bool booleanSolid,
                                  int copyNo);

 protected:
  G4LogicalVolume *body_log;

  G4PVPlacement *inner1_phys;
  G4PVPlacement *inner2_phys;
  G4PVPlacement *central_gap_phys;
  G4PVPlacement *dynode_phys;

  JunoHam20inchPMTConstructionParams fParams;
};

}  // namespace RAT

#endif
