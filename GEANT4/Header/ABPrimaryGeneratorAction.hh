#ifndef ABPrimaryGeneratorAction_hh
#define ABPrimaryGeneratorAction_hh
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
    
class Detectorgeometry;
class G4ParticleGun;
class G4Event;

class ABPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
ABPrimaryGeneratorAction(); //constructor
virtual ~ABPrimaryGeneratorAction(); //destructor
virtual void GeneratePrimaries(G4Event*);

private:
G4ParticleGun* fParticleGun;

};
#endif
