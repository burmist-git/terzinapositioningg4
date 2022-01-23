//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"

#include "globals.hh"

#include "math.h"

//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct(){
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  //bool overlapsChecking = false;
  bool overlapsChecking = true;
  bool buildSmallBoxIncenter = false;
  //bool buildSingleSiPMArray = true;
  //bool buildSingleSiPMArray = false;

  //
  //Satellite envelope
  //
  G4double satellite_envelope_length_x = 552.0*mm;
  G4double satellite_envelope_length_y = 552.0*mm;
  G4double satellite_envelope_length_z = 500.0*mm;

  //
  //Telescope envelope (cylinder)
  //
  G4double telescope_envelope_cylinder_l = 350.0*mm;
  G4double telescope_envelope_cylinder_R = 197.0*mm;
  G4double telescope_envelope_cylinder_x0 = 0.0*mm;
  G4double telescope_envelope_cylinder_y0 = 0.0*mm;
  G4double telescope_envelope_cylinder_z0 = 0.0*mm;
  G4double telescope_envelope_cylinder_theta = (180.0 - 67.5)*deg;
  G4double telescope_envelope_cylinder_phi   = 25.0*deg;
  
  //     
  // World
  //
  G4double world_sizeXY = 100*cm;
  G4double world_sizeZ  = 200*cm;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld = new G4Box("World", 0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,world_mat,"World");                                   
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,                     //no rotation
						   G4ThreeVector(),       //at (0,0,0)
						   logicWorld,            //its logical volume
						   "World",               //its name
						   0,                     //its mother  volume
						   false,                 //no boolean operation
						   0,                     //copy number
						   overlapsChecking);     //overlaps checking

  //
  //Satellite envelope
  //
  G4Box* satellite_envelope_solid = new G4Box("satellite_envelope_solid",
					      satellite_envelope_length_x/2.0,
					      satellite_envelope_length_y/2.0,
					      satellite_envelope_length_z/2.0);
  //
  G4LogicalVolume* satellite_envelope_logic = new G4LogicalVolume(satellite_envelope_solid,world_mat,"satellite_envelope_logic");
  //
  G4RotationMatrix* satellite_envelope_rotation = new G4RotationMatrix();
  G4ThreeVector satellite_envelope_translation( 0.0, 0.0, 0.0);
  //G4VPhysicalVolume* satellite_envelope_phys
  new G4PVPlacement(satellite_envelope_rotation,    //no rotation
		    satellite_envelope_translation, //at (0,0,0)
		    satellite_envelope_logic,       //its logical volume
		    "satellite_envelope_phys",      //its name
		    logicWorld,                     //its mother  volume
		    false,                 //no boolean operation
		    0,                     //copy number
		    overlapsChecking);     //overlaps checking
  
  //
  //Telescope envelope (cylinder)
  //
  G4Tubs *telescope_envelope_cylinder_solid = new G4Tubs("telescope_envelope_cylinder_solid", //name
							 0.0,                                 //pRMin
							 telescope_envelope_cylinder_R,       //pRMax
							 telescope_envelope_cylinder_l/2.0,   //pDz
							 0,                                   //pSPhi
							 (twopi-0.1*deg));                      //pDPhi
  //
  G4LogicalVolume* telescope_envelope_cylinder_logic = new G4LogicalVolume(telescope_envelope_cylinder_solid,
									   world_mat,
									   "telescope_envelope_cylinder_logic");
  //
  G4RotationMatrix telescope_envelope_cylinder_rotation;
  G4ThreeVector telescope_envelope_cylinder_translation;
  G4Transform3D telescope_envelope_cylinder_transformation;
  telescope_envelope_cylinder_translation.setX(telescope_envelope_cylinder_x0);
  telescope_envelope_cylinder_translation.setY(telescope_envelope_cylinder_y0);
  telescope_envelope_cylinder_translation.setZ(telescope_envelope_cylinder_z0);
  telescope_envelope_cylinder_rotation.rotateY(telescope_envelope_cylinder_theta);
  telescope_envelope_cylinder_rotation.rotateZ(telescope_envelope_cylinder_phi);
  telescope_envelope_cylinder_transformation = G4Transform3D(telescope_envelope_cylinder_rotation, telescope_envelope_cylinder_translation);
  //G4VPhysicalVolume* telescope_envelope_cylinder_phys
  new G4PVPlacement(telescope_envelope_cylinder_transformation,  //Transformation
		    telescope_envelope_cylinder_logic,       //its logical volume
		    "telescope_envelope_cylinder_phys",      //its name
		    satellite_envelope_logic,                //its mother  volume
		    false,                                   //no boolean operation
		    0,                                       //copy number
		    overlapsChecking);                       //overlaps checking
  
  //
  //G4RotationMatrix* telescope_envelope_cylinder_rotation = new G4RotationMatrix();
  //telescope_envelope_cylinder_rotation->rotateY(telescope_envelope_cylinder_theta);
  //telescope_envelope_cylinder_rotation->rotateZ(telescope_envelope_cylinder_phi);
  //G4ThreeVector telescope_envelope_cylinder_translation( 0.0, 0.0, 0.0);
  //G4VPhysicalVolume* telescope_envelope_cylinder_phys = new G4PVPlacement(telescope_envelope_cylinder_rotation,    //no rotation
  //									  telescope_envelope_cylinder_translation, //at (0,0,0)
  //									  telescope_envelope_cylinder_logic,       //its logical volume
  //									  "telescope_envelope_cylinder_phys",      //its name
  //									  satellite_envelope_logic,                //its mother  volume
  //									  false,                                   //no boolean operation
  //									  0,                                       //copy number
  //									  overlapsChecking);                       //overlaps checking
 
  //
  // Small box for orientation 
  //
  G4VSolid *boxsmall_solid = new G4Box("boxsmall_solid", 1.0*cm, 3.0*cm, 10.0*cm);
  G4LogicalVolume *boxsmall_logical = new G4LogicalVolume(boxsmall_solid,world_mat,"boxsmall_solid");
  new G4PVPlacement(0,                      //no rotation
		    G4ThreeVector(90*cm/2.0,90*cm/2.0,89*cm),       //at (0,0,0)
		    boxsmall_logical,       //its logical volume
		    "World",                //its name
		    logicWorld,             //its mother  volume
		    false,                  //no boolean operation
		    0,                      //copy number
		    overlapsChecking);      //overlaps checking

  //
  // Small box incenter 
  //
  G4VSolid *boxsmall_centre_solid = new G4Box("boxsmall_centre_solid", 1.0*mm, 1.0*mm, 1.0*mm);
  G4LogicalVolume *boxsmall_centre_logical = new G4LogicalVolume(boxsmall_centre_solid,world_mat,"boxsmall_centre_solid");
  if(buildSmallBoxIncenter)
    new G4PVPlacement(0,                       //no rotation
		      G4ThreeVector(),         //at (0,0,0)
		      boxsmall_centre_logical, //its logical volume
		      "World",                 //its name
		      logicWorld,              //its mother  volume
		      false,                   //no boolean operation
		      0,                       //copy number
		      overlapsChecking);       //overlaps checking
  
  return physWorld;
}

void DetectorConstruction::load_fit_parameters_correction_lens_profile_function(std::string dataFileIn){
  std::ifstream fileIn(dataFileIn.c_str());
  G4double par;
  if (fileIn.is_open()){
    fileIn >> lens_profile_xmin;
    fileIn >> lens_profile_xmax;
    for(int i = 0 ;i<npar_pol8_symmetric;i++){
      fileIn >> par;
      par_pol8_symmetric[i] = par;
    }
    fileIn.close();
  }
  else G4cout<<"Unable to open file"<<G4endl;
}

G4double DetectorConstruction::lens_profile_function(G4double x){
  if(x < lens_profile_xmin || x > lens_profile_xmax )
    return -999.0;
  return par_pol8_symmetric[0] +
    par_pol8_symmetric[1]*pow(x,2) +
    par_pol8_symmetric[2]*pow(x,4) +
    par_pol8_symmetric[3]*pow(x,6) +
    par_pol8_symmetric[4]*pow(x,8);
}

bool DetectorConstruction::get_cons_ring_par(G4double x0, G4double y0, G4double x2, G4double y2, G4double dr, G4double &rMin, G4double &rMax, G4double &conh){
  if( (x2-x0)<0.0){
    G4cout<<"ERROR --> (x2-x0) < 0.0 "<<G4endl
	  <<"                  x2 = "<<x2<<G4endl
      	  <<"                  x0 = "<<x0<<G4endl;
    assert(0);
  }
  double k = (y2 - y0)/(x2-x0);
  double b = y0 - k*x0;
  rMin = x0 - dr;
  rMax = x2 + dr;
  double y0_dr = k*(x0 - dr) + b;
  double y2_dr = k*(x2 + dr) + b;
  conh = abs(y2_dr - y0_dr);
  if(k>0)
    return true;
  return false;
}
