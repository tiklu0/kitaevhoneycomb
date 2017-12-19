//This has been done by Mr.Arnab Barman Roy & Mr. Basudev Sahoo of IITKGP.....
 

#ifndef ABSteppingAction_h
#define ABSteppingAction_h 1

#include "G4UserSteppingAction.hh"

class ABDetectorConstruction;
class ABEventAction;

/// Stepping action class.
///
/// In UserSteppingAction() there are collected the energy deposit and track 
/// lengths of charged particles in Absober and Gap layers and
/// updated in ABEventAction.

class ABSteppingAction : public G4UserSteppingAction
{
public:
  ABSteppingAction(const ABDetectorConstruction* detectorConstruction,
                    ABEventAction* eventAction);
  virtual ~ABSteppingAction();

  virtual void UserSteppingAction(const G4Step* step);
    
private:
  const ABDetectorConstruction* fDetConstruction;
  ABEventAction*  fEventAction;  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
