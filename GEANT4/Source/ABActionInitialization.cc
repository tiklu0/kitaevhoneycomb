//This has been done by Mr.Arnab Barman Roy & Mr. Basudev Sahoo of IITKGP.....

#include "ABActionInitialization.hh"
#include "ABPrimaryGeneratorAction.hh"
#include "ABRunAction.hh"
#include "ABEventAction.hh"
#include "ABSteppingAction.hh"
#include "ABDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ABActionInitialization::ABActionInitialization
                            (ABDetectorConstruction* detConstruction)
 : G4VUserActionInitialization(),
   fDetConstruction(detConstruction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ABActionInitialization::~ABActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABActionInitialization::BuildForMaster() const
{
  SetUserAction(new ABRunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABActionInitialization::Build() const
{
  SetUserAction(new ABPrimaryGeneratorAction);
  SetUserAction(new ABRunAction);
  auto eventAction = new ABEventAction;
  SetUserAction(eventAction);
  SetUserAction(new ABSteppingAction(fDetConstruction,eventAction));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
