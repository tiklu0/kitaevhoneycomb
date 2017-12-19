#include "ABPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh" 
#include <math.h> 
#define PI 3.14159265

ABPrimaryGeneratorAction::ABPrimaryGeneratorAction(): G4VUserPrimaryGeneratorAction(),fParticleGun(0){

G4int n_particle = 1;
fParticleGun = new G4ParticleGun(n_particle);


// Source position

G4double position = -20.0*mm;                                //should be negative in mm, (distance from the source+0.5mm)*(-1)

G4double sinangle= 34.5/(sqrt((position+0.5)*(position+0.5)+(34.5)*(34.5)));


// source initials

G4double sinTheta = sinangle*G4UniformRand();
G4double phi = 2*PI * G4UniformRand();
G4double cosTheta = sqrt(1-sinTheta*sinTheta);



G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
fParticleGun->SetParticleDefinition(particle);
fParticleGun->SetParticleEnergy(1332*keV);
fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sinTheta * cos(phi) , sinTheta * sin(phi), cosTheta));



fParticleGun->SetParticlePosition(G4ThreeVector(0*mm, 0*mm, position));
 }
    
    
    
ABPrimaryGeneratorAction::~ABPrimaryGeneratorAction(){
      
delete fParticleGun;
}
    
    

    
void ABPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){ 
    
      
G4ThreeVector position = fParticleGun->GetParticlePosition();
      
      
fParticleGun->SetParticlePosition(position);
      
fParticleGun->GeneratePrimaryVertex(anEvent);
}





