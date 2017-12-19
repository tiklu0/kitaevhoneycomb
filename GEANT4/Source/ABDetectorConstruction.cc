//This has been done by Mr.Arnab Barman Roy & Mr. Basudev Sahoo of IITKGP.....

#include "ABDetectorConstruction.hh"
//#include "G4PhysicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4Tubs.hh"
#include "G4Torus.hh"
#include "G4Sphere.hh" 

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* ABDetectorConstruction::fMagFieldMessenger = nullptr; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ABDetectorConstruction::ABDetectorConstruction()
 : G4VUserDetectorConstruction(),
   fAbsorberPV(nullptr),
   fGapPV(nullptr),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ABDetectorConstruction::~ABDetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ABDetectorConstruction::Construct()
{
  // Define materials 
  //DefineMaterials();
  
  // Define volumes
  //return DefineVolumes();


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4NistManager* nist = G4NistManager::Instance();
G4Material* world_mat=nist->FindOrBuildMaterial("G4_AIR");
G4Material* cover_mat=nist->FindOrBuildMaterial("G4_Al");
G4Material* det_mat=nist->FindOrBuildMaterial("G4_Ge");




//the required parameters
G4double q=34.5*mm;
G4double dt=0.7*mm;



//world
G4Box* sworld =new G4Box("World",40*mm,40*mm,60*mm); 

G4LogicalVolume* lworld=new G4LogicalVolume(sworld,world_mat,"world");

G4VPhysicalVolume* pworld=new G4PVPlacement(0,G4ThreeVector(),lworld,"world",0,false,0);




//the solids

//part1

G4VSolid* part1=new G4Tubs("part1cylinder",0,q-1*mm,1*mm,0*degree,360*degree);

//part2

G4VSolid* part2=new G4Tubs("part2cylinder",q-1*mm,q,56*mm,0*degree,360*degree);

//part3

G4VSolid* part3=new G4Tubs("part3cylinder",0*mm,28.5*mm,0.03*mm,0*degree,360*degree);

//part4

G4VSolid* part4=new G4Tubs("part4cylinder",28.5*mm,29.3*mm,56*mm,0*degree,360*degree); //56->26

//part5

   //complexobject3

G4VSolid* cobj3toroid=new G4Torus("cobj3t1",8*mm-dt,8*mm,18.3*mm,0*degree,360*degree);

G4VSolid* cobj3cylinder1=new G4Tubs("cobj3c1",0*mm,30*mm,8*mm,0*degree,360*degree);


G4VSolid* cobj3cylinder2=new G4Tubs("cobj3c2",0*mm,18.3*mm,8*mm,0*degree,360*degree);


G4VSolid* cobj3det1= new G4SubtractionSolid("cobj3Detector1", cobj3toroid, cobj3cylinder1,0, G4ThreeVector(0,0,-8*mm));


G4VSolid* cobj3det2= new G4SubtractionSolid("cobj3Detector2", cobj3toroid, cobj3cylinder2,0,G4ThreeVector(0,0,0));

    //complexobject3

G4VSolid* part5cylinder1=new G4Tubs("part5c1",0*mm,18.3*mm,dt,0*degree,360*degree);

G4VSolid* part5cylinder2= new G4Tubs("part5c2",26.3*mm-dt, 26.3*mm, 21.9*mm,0*degree,360*degree);

G4VSolid* protopart5= new G4UnionSolid("protopart5m", cobj3det2, part5cylinder1,0, G4ThreeVector(0,0,8*mm));

G4VSolid* part5= new G4UnionSolid("part5m", protopart5, part5cylinder2,0, G4ThreeVector(0,0,-21.9*mm));


//part6

    //complexobject1

G4VSolid* cobj1sphere=new G4Sphere("cobj1shell",0*mm,5*mm,0*degree,360*degree,0,64.158*degree);

G4VSolid* cobj1cylinder=new G4Tubs("cobj1cylinder",0*mm,4.5*mm,17.74*mm,0*degree,360*degree);

G4VSolid* cobj1det2=new G4Tubs("cobj1det",0*mm,26.3*mm-dt,21.9*mm,0*degree,360*degree);

G4VSolid* cobj1union= new G4UnionSolid("cobj1Sphere+Cylinder",cobj1cylinder,cobj1sphere,0,G4ThreeVector(0,0,15.561*mm)); 

G4VSolid* cobj1det= new G4SubtractionSolid("cobj1Detector", cobj1det2, cobj1union,0, G4ThreeVector(0,0,-4.16*mm));

    

    //complexobject2

G4VSolid* cobj2toroid=new G4Torus("cobj2t1",0,8*mm-dt,18.3*mm,0*degree,360*degree);

G4VSolid* cobj2cylinder1=new G4Tubs("cobj2c1",0*mm,30*mm,8*mm,0*degree,360*degree);


G4VSolid* cobj2cylinder2=new G4Tubs("cobj2c2",0*mm,18.3*mm,8*mm,0*degree,360*degree);


G4VSolid* cobj2det1= new G4SubtractionSolid("cobj2Detector1", cobj2toroid, cobj2cylinder1,0, G4ThreeVector(0,0,-8*mm));


G4VSolid* cobj2det2= new G4SubtractionSolid("cobj2Detector2", cobj2toroid, cobj2cylinder2,0,G4ThreeVector(0,0,0));

     

G4VSolid* vcobj1cylinder=new G4Tubs("vcobj1cylinder",0*mm,18.3*mm,4*mm-dt*0.5,0*degree,360*degree);

G4VSolid* vcobj1det1= new G4UnionSolid("vcobj1Detector1", cobj2det2, vcobj1cylinder,0, G4ThreeVector(0,0,4*mm-0.5*dt));

G4VSolid* part6= new G4UnionSolid("part6", vcobj1det1, cobj1det,0,G4ThreeVector(0,0,-21.9*mm));




//The Logical Volumes

G4LogicalVolume* lpart1=new G4LogicalVolume(part1,cover_mat,"part1",0,0,0);

G4LogicalVolume* lpart2=new G4LogicalVolume(part2,cover_mat,"part2",0,0,0);

G4LogicalVolume* lpart3=new G4LogicalVolume(part3,cover_mat,"part3",0,0,0);

G4LogicalVolume* lpart4=new G4LogicalVolume(part4,cover_mat,"part4",0,0,0);

G4LogicalVolume* lpart5=new G4LogicalVolume(part5,det_mat,"part5",0,0,0);

G4LogicalVolume* lpart6=new G4LogicalVolume(part6,det_mat,"part6",0,0,0);



//The physical volumes


G4RotationMatrix Ra;
Ra.rotateY(0*deg);

//part1

G4ThreeVector uz1= G4ThreeVector(0,0,(55+3.515)*mm);
   
G4ThreeVector position = uz1;

G4Transform3D transform = G4Transform3D(Ra,position);

G4VPhysicalVolume* ppart1=new G4PVPlacement(transform,lpart1,"ppart1",lworld,false,0);



//part2

G4ThreeVector uz2 = G4ThreeVector(0,0,55.5*mm);
   
G4ThreeVector position2 = uz2;

G4Transform3D transform2 = G4Transform3D(Ra,position2);

G4VPhysicalVolume* ppart2=new G4PVPlacement(transform2,lpart2,"ppart2",lworld,false,0);



//part3

G4ThreeVector uz3 = G4ThreeVector(0,0,(55+3.515)*mm);
   
G4ThreeVector position3 = uz3;

G4Transform3D transform3 = G4Transform3D(Ra,position3);

G4VPhysicalVolume* ppart3=new G4PVPlacement(transform3,lpart3,"ppart3",lworld,false,0);



//part4

G4ThreeVector uz4 = G4ThreeVector(0,0,55.5*mm);
   
G4ThreeVector position4 = uz4;

G4Transform3D transform4 = G4Transform3D(Ra,position4);

G4VPhysicalVolume* ppart4=new G4PVPlacement(transform4,lpart4,"ppart4",lworld,false,0);



//part5

G4ThreeVector uz5 = G4ThreeVector(0,0,50*mm);   //11.530 mm
   
G4ThreeVector position5 = uz5;

G4Transform3D transform5 = G4Transform3D(Ra,position5);

G4VPhysicalVolume* ppart5 =new G4PVPlacement(transform5,lpart5,"fAbsorberPV",lworld,false,0);
fAbsorberPV=ppart5;


//part6

G4ThreeVector uz6 = G4ThreeVector(0,0,50*mm);
   
G4ThreeVector position6 = uz6;

G4Transform3D transform6 = G4Transform3D(Ra,position6);

G4VPhysicalVolume* ppart6=new G4PVPlacement(transform6,lpart6,"fGapPV",lworld,false,0);
fGapPV=ppart6;




//Visualization attributes

G4VisAttributes* Red = new G4VisAttributes( G4Colour(255/255. ,0/255. ,0/255.));
   
G4VisAttributes* Yellow= new G4VisAttributes( G4Colour(255/255. ,255/255. ,0/255.));
   
G4VisAttributes* LightBlue = new G4VisAttributes( G4Colour(0/255.   ,204/255. ,204/255.));

G4VisAttributes* LightGreen = new G4VisAttributes( G4Colour(153/255. ,255/255. ,153/255.));



  //
  // print parameters
  //
//  G4cout
  //  << G4endl 
    //<< "------------------------------------------------------------" << G4endl
    //<< "---> The calorimeter is " << nofLayers << " layers of: [ "
   // << absoThickness/mm << "mm of " << absorberMaterial->GetName() 
    //<< " + "
    //<< gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
    //<< "------------------------------------------------------------" << G4endl;
  
                                        
   //Visualization attributes
  
  lworld->SetVisAttributes (G4VisAttributes::GetInvisible());

  //auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  Red->SetVisibility(true);
  LightBlue->SetVisibility(true);
  lpart6->SetVisAttributes(Red);
  lpart5->SetVisAttributes(Red);
  lpart1->SetVisAttributes(LightBlue);
  lpart2->SetVisAttributes(LightBlue);
  lpart3->SetVisAttributes(LightBlue);
  lpart4->SetVisAttributes(LightBlue);
  //
  // Always return the physical World
  //
  return pworld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABDetectorConstruction::ConstructSDandField()
{ 
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
