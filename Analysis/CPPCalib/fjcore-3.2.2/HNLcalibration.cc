//----------------------------------------------------------------------
/// \file
/// \page Example01 01 - basic usage example
///
/// fastjet basic example program:
///   simplest illustration of the usage of the basic classes:
///   fjcore::PseudoJet, fjcore::JetDefinition and 
///   fjcore::ClusterSequence
///
/// run it with    : ./01-basic < single-event.dat
///
/// Source code: 01-basic.cc
//----------------------------------------------------------------------

//STARTHEADER
// $Id: 01-basic.cc 2684 2011-11-14 07:41:44Z soyez $
//
// Copyright (c) 2005-2011, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
#include "fjcore.hh"


#include <iostream> // needed for io
#include <cstdio>   // needed for io

using namespace std;

/// an example program showing how to use fastjet
int main(){
  
  // read in input particles
  //----------------------------------------------------------
  vector<fjcore::PseudoJet> input_particles;
  
  double px, py , pz, E;
 
        

  input_particles.push_back(fjcore::PseudoJet(141.760364,773.876369,149.110431,931.072972)); 
  input_particles.push_back(fjcore::PseudoJet(-141.760364,-773.876369,571.764090,972.571236)); 
  input_particles.push_back(fjcore::PseudoJet(104.453936,439.022975,-155.044174,477.169250)); 
  input_particles.push_back(fjcore::PseudoJet(37.306428,334.853394,304.154605,453.903722)); 
  


  // label the columns
  printf("%5s %15s %15s %15s\n","jet #", "pt", "eta", "phi");
 
  // print out the details for each jet
  for (unsigned int i = 0; i < input_particles.size(); i++) {
    printf("%5u %15.8f %15.8f %15.8f\n",
	   i, input_particles[i].perp(), input_particles[i].eta(), input_particles[i].phi() );
  }

  cout << "delta_R(3,0) " << input_particles[3].delta_R(input_particles[0]) << endl;
  cout << "delta_Phi(3,0) " << input_particles[3].delta_phi_to(input_particles[0]) << endl;
  cout << "mass(3+0) " << (input_particles[3]+input_particles[0]).m() << endl;

  return 0;
}
