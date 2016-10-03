#include "Riostream.h"
#include "TMath.h"
#include "TH1.h"
#include "THStack.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLine.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <limits>

void dup_remover()
{
  TFile *f = TFile::Open("trees/merge.root");
  TFile *newfile = new TFile("trees/2nu_2MeV_tree.root","recreate");

  TTree *oldtree = (TTree*)f->Get("snemodata");

  oldtree->SetBranchStatus("*",1);
  oldtree->SetBranchStatus("*version*",0);
  oldtree->SetBranchStatus("2e1g*",0);
  oldtree->SetBranchStatus("2e2g*",0);
  oldtree->SetBranchStatus("2e3g*",0);

  //Create a new file + a clone of old tree in new file

  TTree *newtree = oldtree->CloneTree();
  oldtree->SetBranchStatus("2e1g*",1);
  oldtree->SetBranchStatus("2e2g*",1);
  // oldtree->SetBranchStatus("2e3g*",1);

  //At this point, got all the non 2eNg topologies

  Double_t topo_2e1g_electron_minimal_energy;
  Double_t topo_2e1g_electron_maximal_energy;
  Double_t topo_2e1g_gamma_energy;
  Double_t topo_2e1g_electrons_energy_difference;
  Double_t topo_2e1g_electrons_energy_sum;
  Double_t topo_2e1g_electrons_gammas_energy_sum;
  Double_t topo_2e1g_electrons_internal_probability;
  Double_t topo_2e1g_electrons_external_probability;
  Double_t topo_2e1g_electron_min_gamma_internal_probability;
  Double_t topo_2e1g_electron_min_gamma_external_probability;
  Double_t topo_2e1g_electron_max_gamma_internal_probability;
  Double_t topo_2e1g_electron_max_gamma_external_probability;
  Double_t topo_2e1g_electrons_vertices_probability;
  Double_t topo_2e1g_electrons_angle;
  Double_t topo_2e1g_electrons_cos_angle;
  Double_t topo_2e1g_electron_Emin_track_length;
  Double_t topo_2e1g_electron_Emax_track_length;

  Double_t topo_2e2g_electron_minimal_energy;
  Double_t topo_2e2g_electron_maximal_energy;
  Double_t topo_2e2g_gamma_min_energy;
  Double_t topo_2e2g_gamma_max_energy;
  Double_t topo_2e2g_electrons_energy_difference;
  Double_t topo_2e2g_electrons_energy_sum;
  Double_t topo_2e2g_electrons_gammas_energy_sum;
  Double_t topo_2e2g_electrons_internal_probability;
  Double_t topo_2e2g_electrons_external_probability;
  Double_t topo_2e2g_electron_min_gamma_min_internal_probability;
  Double_t topo_2e2g_electron_min_gamma_max_internal_probability;
  Double_t topo_2e2g_electron_min_gamma_min_external_probability;
  Double_t topo_2e2g_electron_min_gamma_max_external_probability;
  Double_t topo_2e2g_electron_max_gamma_min_internal_probability;
  Double_t topo_2e2g_electron_max_gamma_max_internal_probability;
  Double_t topo_2e2g_electron_max_gamma_min_external_probability;
  Double_t topo_2e2g_electron_max_gamma_max_external_probability;
  Double_t topo_2e2g_electrons_vertices_probability;
  Double_t topo_2e2g_electrons_angle;
  Double_t topo_2e2g_electrons_cos_angle;
  Double_t topo_2e2g_electron_Emin_track_length;
  Double_t topo_2e2g_electron_Emax_track_length;

  // Double_t topo_2e3g_electron_minimal_energy;
  // Double_t topo_2e3g_electron_maximal_energy;
  // Double_t topo_2e3g_gamma_min_energy;
  // Double_t topo_2e3g_gamma_mid_energy;
  // Double_t topo_2e3g_gamma_max_energy;
  // Double_t topo_2e3g_electrons_energy_difference;
  // Double_t topo_2e3g_electrons_energy_sum;
  // Double_t topo_2e3g_electrons_gammas_energy_sum;
  // Double_t topo_2e3g_electrons_internal_probability;
  // Double_t topo_2e3g_electrons_external_probability;
  // Double_t topo_2e3g_electron_min_gamma_min_internal_probability;
  // Double_t topo_2e3g_electron_min_gamma_min_external_probability;
  // Double_t topo_2e3g_electron_min_gamma_mid_internal_probability;
  // Double_t topo_2e3g_electron_min_gamma_mid_external_probability;
  // Double_t topo_2e3g_electron_min_gamma_max_internal_probability;
  // Double_t topo_2e3g_electron_min_gamma_max_external_probability;
  // Double_t topo_2e3g_electron_max_gamma_min_internal_probability;
  // Double_t topo_2e3g_electron_max_gamma_min_external_probability;
  // Double_t topo_2e3g_electron_max_gamma_mid_internal_probability;
  // Double_t topo_2e3g_electron_max_gamma_mid_external_probability;
  // Double_t topo_2e3g_electron_max_gamma_max_internal_probability;
  // Double_t topo_2e3g_electron_max_gamma_max_external_probability;
  // Double_t topo_2e3g_electrons_vertices_probability;
  // Double_t topo_2e3g_electrons_angle;
  // Double_t topo_2e3g_electrons_cos_angle;
  // Double_t topo_2e3g_electron_Emin_track_length;
  // Double_t topo_2e3g_electron_Emax_track_length;

  oldtree->SetBranchAddress("2e1g_electron_minimal_energy",&topo_2e1g_electron_minimal_energy);
  oldtree->SetBranchAddress("2e1g_electron_maximal_energy",&topo_2e1g_electron_maximal_energy);
  oldtree->SetBranchAddress("2e1g_gamma_energy",&topo_2e1g_gamma_energy);
  oldtree->SetBranchAddress("2e1g_electrons_energy_difference",&topo_2e1g_electrons_energy_difference);
  oldtree->SetBranchAddress("2e1g_electrons_energy_sum",&topo_2e1g_electrons_energy_sum);
  oldtree->SetBranchAddress("2e1g_electrons_gammas_energy_sum",&topo_2e1g_electrons_gammas_energy_sum);
  oldtree->SetBranchAddress("2e1g_electrons_internal_probability",&topo_2e1g_electrons_internal_probability);
  oldtree->SetBranchAddress("2e1g_electrons_external_probability",&topo_2e1g_electrons_external_probability);
  oldtree->SetBranchAddress("2e1g_electron_min_gamma_internal_probability",&topo_2e1g_electron_min_gamma_internal_probability);
  oldtree->SetBranchAddress("2e1g_electron_min_gamma_external_probability",&topo_2e1g_electron_min_gamma_external_probability);
  oldtree->SetBranchAddress("2e1g_electron_max_gamma_internal_probability",&topo_2e1g_electron_max_gamma_internal_probability);
  oldtree->SetBranchAddress("2e1g_electron_max_gamma_external_probability",&topo_2e1g_electron_max_gamma_external_probability);
  oldtree->SetBranchAddress("2e1g_electrons_vertices_probability",&topo_2e1g_electrons_vertices_probability);
  oldtree->SetBranchAddress("2e1g_electrons_angle",&topo_2e1g_electrons_angle);
  oldtree->SetBranchAddress("2e1g_electrons_cos_angle",&topo_2e1g_electrons_cos_angle);
  oldtree->SetBranchAddress("2e1g_electron_Emin_track_length",&topo_2e1g_electron_Emin_track_length);
  oldtree->SetBranchAddress("2e1g_electron_Emax_track_length",&topo_2e1g_electron_Emax_track_length);

  oldtree->SetBranchAddress("2e2g_electron_minimal_energy",&topo_2e2g_electron_minimal_energy);
  oldtree->SetBranchAddress("2e2g_electron_maximal_energy",&topo_2e2g_electron_maximal_energy);
  oldtree->SetBranchAddress("2e2g_gamma_min_energy",&topo_2e2g_gamma_min_energy);
  oldtree->SetBranchAddress("2e2g_gamma_max_energy",&topo_2e2g_gamma_max_energy);
  oldtree->SetBranchAddress("2e2g_electrons_energy_difference",&topo_2e2g_electrons_energy_difference);
  oldtree->SetBranchAddress("2e2g_electrons_energy_sum",&topo_2e2g_electrons_energy_sum);
  oldtree->SetBranchAddress("2e2g_electrons_gammas_energy_sum",&topo_2e2g_electrons_gammas_energy_sum);
  oldtree->SetBranchAddress("2e2g_electrons_internal_probability",&topo_2e2g_electrons_internal_probability);
  oldtree->SetBranchAddress("2e2g_electrons_external_probability",&topo_2e2g_electrons_external_probability);
  oldtree->SetBranchAddress("2e2g_electron_min_gamma_min_internal_probability",&topo_2e2g_electron_min_gamma_min_internal_probability);
  oldtree->SetBranchAddress("2e2g_electron_min_gamma_max_internal_probability",&topo_2e2g_electron_min_gamma_max_internal_probability);
  oldtree->SetBranchAddress("2e2g_electron_min_gamma_min_external_probability",&topo_2e2g_electron_min_gamma_min_external_probability);
  oldtree->SetBranchAddress("2e2g_electron_min_gamma_max_external_probability",&topo_2e2g_electron_min_gamma_max_external_probability);
  oldtree->SetBranchAddress("2e2g_electron_max_gamma_min_internal_probability",&topo_2e2g_electron_max_gamma_min_internal_probability);
  oldtree->SetBranchAddress("2e2g_electron_max_gamma_max_internal_probability",&topo_2e2g_electron_max_gamma_max_internal_probability);
  oldtree->SetBranchAddress("2e2g_electron_max_gamma_min_external_probability",&topo_2e2g_electron_max_gamma_min_external_probability);
  oldtree->SetBranchAddress("2e2g_electron_max_gamma_max_external_probability",&topo_2e2g_electron_max_gamma_max_external_probability);
  oldtree->SetBranchAddress("2e2g_electrons_vertices_probability",&topo_2e2g_electrons_vertices_probability);
  oldtree->SetBranchAddress("2e2g_electrons_angle",&topo_2e2g_electrons_angle);
  oldtree->SetBranchAddress("2e2g_electrons_cos_angle",&topo_2e2g_electrons_cos_angle);
  oldtree->SetBranchAddress("2e2g_electron_Emin_track_length",&topo_2e2g_electron_Emin_track_length);
  oldtree->SetBranchAddress("2e2g_electron_Emax_track_length",&topo_2e2g_electron_Emax_track_length);

  // oldtree->SetBranchAddress("2e3g_electron_minimal_energy",&topo_2e3g_electron_minimal_energy);
  // oldtree->SetBranchAddress("2e3g_electron_maximal_energy",&topo_2e3g_electron_maximal_energy);
  // oldtree->SetBranchAddress("2e3g_gamma_min_energy",&topo_2e3g_gamma_min_energy);
  // oldtree->SetBranchAddress("2e3g_gamma_mid_energy",&topo_2e3g_gamma_mid_energy);
  // oldtree->SetBranchAddress("2e3g_gamma_max_energy",&topo_2e3g_gamma_max_energy);
  // oldtree->SetBranchAddress("2e3g_electrons_energy_difference",&topo_2e3g_electrons_energy_difference);
  // oldtree->SetBranchAddress("2e3g_electrons_energy_sum",&topo_2e3g_electrons_energy_sum);
  // oldtree->SetBranchAddress("2e3g_electrons_gammas_energy_sum",&topo_2e3g_electrons_gammas_energy_sum);
  // oldtree->SetBranchAddress("2e3g_electrons_internal_probability",&topo_2e3g_electrons_internal_probability);
  // oldtree->SetBranchAddress("2e3g_electrons_external_probability",&topo_2e3g_electrons_external_probability);
  // oldtree->SetBranchAddress("2e3g_electron_min_gamma_min_internal_probability",&topo_2e3g_electron_min_gamma_min_internal_probability);
  // oldtree->SetBranchAddress("2e3g_electron_min_gamma_min_external_probability",&topo_2e3g_electron_min_gamma_min_external_probability);
  // oldtree->SetBranchAddress("2e3g_electron_min_gamma_mid_internal_probability",&topo_2e3g_electron_min_gamma_mid_internal_probability);
  // oldtree->SetBranchAddress("2e3g_electron_min_gamma_mid_external_probability",&topo_2e3g_electron_min_gamma_mid_external_probability);
  // oldtree->SetBranchAddress("2e3g_electron_min_gamma_max_internal_probability",&topo_2e3g_electron_min_gamma_max_internal_probability);
  // oldtree->SetBranchAddress("2e3g_electron_min_gamma_max_external_probability",&topo_2e3g_electron_min_gamma_max_external_probability);
  // oldtree->SetBranchAddress("2e3g_electron_max_gamma_min_internal_probability",&topo_2e3g_electron_max_gamma_min_internal_probability);
  // oldtree->SetBranchAddress("2e3g_electron_max_gamma_min_external_probability",&topo_2e3g_electron_max_gamma_min_external_probability);
  // oldtree->SetBranchAddress("2e3g_electron_max_gamma_mid_internal_probability",&topo_2e3g_electron_max_gamma_mid_internal_probability);
  // oldtree->SetBranchAddress("2e3g_electron_max_gamma_mid_external_probability",&topo_2e3g_electron_max_gamma_mid_external_probability);
  // oldtree->SetBranchAddress("2e3g_electron_max_gamma_max_internal_probability",&topo_2e3g_electron_max_gamma_max_internal_probability);
  // oldtree->SetBranchAddress("2e3g_electron_max_gamma_max_external_probability",&topo_2e3g_electron_max_gamma_max_external_probability);
  // oldtree->SetBranchAddress("2e3g_electrons_vertices_probability",&topo_2e3g_electrons_vertices_probability);
  // oldtree->SetBranchAddress("2e3g_electrons_angle",&topo_2e3g_electrons_angle);
  // oldtree->SetBranchAddress("2e3g_electrons_cos_angle",&topo_2e3g_electrons_cos_angle);
  // oldtree->SetBranchAddress("2e3g_electron_Emin_track_length",&topo_2e3g_electron_Emin_track_length);
  // oldtree->SetBranchAddress("2e3g_electron_Emax_track_length",&topo_2e3g_electron_Emax_track_length);

  Double_t new_topo_2e1g_electron_minimal_energy;
  Double_t new_topo_2e1g_electron_maximal_energy;
  Double_t new_topo_2e1g_gamma_energy;
  Double_t new_topo_2e1g_electrons_energy_difference;
  Double_t new_topo_2e1g_electrons_energy_sum;
  Double_t new_topo_2e1g_electrons_gammas_energy_sum;
  Double_t new_topo_2e1g_electrons_internal_probability;
  Double_t new_topo_2e1g_electrons_external_probability;
  Double_t new_topo_2e1g_electron_min_gamma_internal_probability;
  Double_t new_topo_2e1g_electron_min_gamma_external_probability;
  Double_t new_topo_2e1g_electron_max_gamma_internal_probability;
  Double_t new_topo_2e1g_electron_max_gamma_external_probability;
  Double_t new_topo_2e1g_electrons_vertices_probability;
  Double_t new_topo_2e1g_electrons_angle;
  Double_t new_topo_2e1g_electrons_cos_angle;
  Double_t new_topo_2e1g_electron_Emin_track_length;
  Double_t new_topo_2e1g_electron_Emax_track_length;

  Double_t new_topo_2e2g_electron_minimal_energy;
  Double_t new_topo_2e2g_electron_maximal_energy;
  Double_t new_topo_2e2g_gamma_min_energy;
  Double_t new_topo_2e2g_gamma_max_energy;
  Double_t new_topo_2e2g_electrons_energy_difference;
  Double_t new_topo_2e2g_electrons_energy_sum;
  Double_t new_topo_2e2g_electrons_gammas_energy_sum;
  Double_t new_topo_2e2g_electrons_internal_probability;
  Double_t new_topo_2e2g_electrons_external_probability;
  Double_t new_topo_2e2g_electron_min_gamma_min_internal_probability;
  Double_t new_topo_2e2g_electron_min_gamma_max_internal_probability;
  Double_t new_topo_2e2g_electron_min_gamma_min_external_probability;
  Double_t new_topo_2e2g_electron_min_gamma_max_external_probability;
  Double_t new_topo_2e2g_electron_max_gamma_min_internal_probability;
  Double_t new_topo_2e2g_electron_max_gamma_max_internal_probability;
  Double_t new_topo_2e2g_electron_max_gamma_min_external_probability;
  Double_t new_topo_2e2g_electron_max_gamma_max_external_probability;
  Double_t new_topo_2e2g_electrons_vertices_probability;
  Double_t new_topo_2e2g_electrons_angle;
  Double_t new_topo_2e2g_electrons_cos_angle;
  Double_t new_topo_2e2g_electron_Emin_track_length;
  Double_t new_topo_2e2g_electron_Emax_track_length;

  // Double_t new_topo_2e3g_electron_minimal_energy;
  // Double_t new_topo_2e3g_electron_maximal_energy;
  // Double_t new_topo_2e3g_gamma_min_energy;
  // Double_t new_topo_2e3g_gamma_mid_energy;
  // Double_t new_topo_2e3g_gamma_max_energy;
  // Double_t new_topo_2e3g_electrons_energy_difference;
  // Double_t new_topo_2e3g_electrons_energy_sum;
  // Double_t new_topo_2e3g_electrons_gammas_energy_sum;
  // Double_t new_topo_2e3g_electrons_internal_probability;
  // Double_t new_topo_2e3g_electrons_external_probability;
  // Double_t new_topo_2e3g_electron_min_gamma_min_internal_probability;
  // Double_t new_topo_2e3g_electron_min_gamma_min_external_probability;
  // Double_t new_topo_2e3g_electron_min_gamma_mid_internal_probability;
  // Double_t new_topo_2e3g_electron_min_gamma_mid_external_probability;
  // Double_t new_topo_2e3g_electron_min_gamma_max_internal_probability;
  // Double_t new_topo_2e3g_electron_min_gamma_max_external_probability;
  // Double_t new_topo_2e3g_electron_max_gamma_min_internal_probability;
  // Double_t new_topo_2e3g_electron_max_gamma_min_external_probability;
  // Double_t new_topo_2e3g_electron_max_gamma_mid_internal_probability;
  // Double_t new_topo_2e3g_electron_max_gamma_mid_external_probability;
  // Double_t new_topo_2e3g_electron_max_gamma_max_internal_probability;
  // Double_t new_topo_2e3g_electron_max_gamma_max_external_probability;
  // Double_t new_topo_2e3g_electrons_vertices_probability;
  // Double_t new_topo_2e3g_electrons_angle;
  // Double_t new_topo_2e3g_electrons_cos_angle;
  // Double_t new_topo_2e3g_electron_Emin_track_length;
  // Double_t new_topo_2e3g_electron_Emax_track_length;

  newtree->Branch("2e1g_electron_minimal_energy",&new_topo_2e1g_electron_minimal_energy,"new_topo_2e1g_electron_minimal_energy/D");
  newtree->Branch("2e1g_electron_maximal_energy",&new_topo_2e1g_electron_maximal_energy,"new_topo_2e1g_electron_maximal_energy/D");
  newtree->Branch("2e1g_gamma_energy",&new_topo_2e1g_gamma_energy,"new_topo_2e1g_gamma_energy/D");
  newtree->Branch("2e1g_electrons_energy_difference",&new_topo_2e1g_electrons_energy_difference,"new_topo_2e1g_electrons_energy_difference/D");
  newtree->Branch("2e1g_electrons_energy_sum",&new_topo_2e1g_electrons_energy_sum,"new_topo_2e1g_electrons_energy_sum/D");
  newtree->Branch("2e1g_electrons_gammas_energy_sum",&new_topo_2e1g_electrons_gammas_energy_sum,"new_topo_2e1g_electrons_gammas_energy_sum/D");
  newtree->Branch("2e1g_electrons_internal_probability",&new_topo_2e1g_electrons_internal_probability,"new_topo_2e1g_electrons_internal_probability/D");
  newtree->Branch("2e1g_electrons_external_probability",&new_topo_2e1g_electrons_external_probability,"new_topo_2e1g_electrons_external_probability/D");
  newtree->Branch("2e1g_electron_min_gamma_internal_probability",&new_topo_2e1g_electron_min_gamma_internal_probability,"new_topo_2e1g_electron_min_gamma_internal_probability/D");
  newtree->Branch("2e1g_electron_min_gamma_external_probability",&new_topo_2e1g_electron_min_gamma_external_probability,"new_topo_2e1g_electron_min_gamma_external_probability/D");
  newtree->Branch("2e1g_electron_max_gamma_internal_probability",&new_topo_2e1g_electron_max_gamma_internal_probability,"new_topo_2e1g_electron_max_gamma_internal_probability/D");
  newtree->Branch("2e1g_electron_max_gamma_external_probability",&new_topo_2e1g_electron_max_gamma_external_probability,"new_topo_2e1g_electron_max_gamma_external_probability/D");
  newtree->Branch("2e1g_electrons_vertices_probability",&new_topo_2e1g_electrons_vertices_probability,"new_topo_2e1g_electrons_vertices_probability/D");
  newtree->Branch("2e1g_electrons_angle",&new_topo_2e1g_electrons_angle,"new_topo_2e1g_electrons_angle/D");
  newtree->Branch("2e1g_electrons_cos_angle",&new_topo_2e1g_electrons_cos_angle,"new_topo_2e1g_electrons_cos_angle/D");
  newtree->Branch("2e1g_electron_Emin_track_length",&new_topo_2e1g_electron_Emin_track_length,"new_topo_2e1g_electron_Emin_track_length/D");
  newtree->Branch("2e1g_electron_Emax_track_length",&new_topo_2e1g_electron_Emax_track_length,"new_topo_2e1g_electron_Emax_track_length/D");

  newtree->Branch("2e2g_electron_minimal_energy",&new_topo_2e2g_electron_minimal_energy,"new_topo_2e2g_electron_minimal_energy/D");
  newtree->Branch("2e2g_electron_maximal_energy",&new_topo_2e2g_electron_maximal_energy,"new_topo_2e2g_electron_maximal_energy/D");
  newtree->Branch("2e2g_gamma_min_energy",&new_topo_2e2g_gamma_min_energy,"new_topo_2e2g_gamma_min_energy/D");
  newtree->Branch("2e2g_gamma_max_energy",&new_topo_2e2g_gamma_max_energy,"new_topo_2e2g_gamma_max_energy/D");
  newtree->Branch("2e2g_electrons_energy_difference",&new_topo_2e2g_electrons_energy_difference,"new_topo_2e2g_electrons_energy_difference/D");
  newtree->Branch("2e2g_electrons_energy_sum",&new_topo_2e2g_electrons_energy_sum,"new_topo_2e2g_electrons_energy_sum/D");
  newtree->Branch("2e2g_electrons_gammas_energy_sum",&new_topo_2e2g_electrons_gammas_energy_sum,"new_topo_2e2g_electrons_gammas_energy_sum/D");
  newtree->Branch("2e2g_electrons_internal_probability",&new_topo_2e2g_electrons_internal_probability,"new_topo_2e2g_electrons_internal_probability/D");
  newtree->Branch("2e2g_electrons_external_probability",&new_topo_2e2g_electrons_external_probability,"new_topo_2e2g_electrons_external_probability/D");
  newtree->Branch("2e2g_electron_min_gamma_min_internal_probability",&new_topo_2e2g_electron_min_gamma_min_internal_probability,"new_topo_2e2g_electron_min_gamma_min_internal_probability/D");
  newtree->Branch("2e2g_electron_min_gamma_max_internal_probability",&new_topo_2e2g_electron_min_gamma_max_internal_probability,"new_topo_2e2g_electron_min_gamma_max_internal_probability/D");
  newtree->Branch("2e2g_electron_min_gamma_min_external_probability",&new_topo_2e2g_electron_min_gamma_min_external_probability,"new_topo_2e2g_electron_min_gamma_min_external_probability/D");
  newtree->Branch("2e2g_electron_min_gamma_max_external_probability",&new_topo_2e2g_electron_min_gamma_max_external_probability,"new_topo_2e2g_electron_min_gamma_max_external_probability/D");
  newtree->Branch("2e2g_electron_max_gamma_min_internal_probability",&new_topo_2e2g_electron_max_gamma_min_internal_probability,"new_topo_2e2g_electron_max_gamma_min_internal_probability/D");
  newtree->Branch("2e2g_electron_max_gamma_max_internal_probability",&new_topo_2e2g_electron_max_gamma_max_internal_probability,"new_topo_2e2g_electron_max_gamma_max_internal_probability/D");
  newtree->Branch("2e2g_electron_max_gamma_min_external_probability",&new_topo_2e2g_electron_max_gamma_min_external_probability,"new_topo_2e2g_electron_max_gamma_min_external_probability/D");
  newtree->Branch("2e2g_electron_max_gamma_max_external_probability",&new_topo_2e2g_electron_max_gamma_max_external_probability,"new_topo_2e2g_electron_max_gamma_max_external_probability/D");
  newtree->Branch("2e2g_electrons_vertices_probability",&new_topo_2e2g_electrons_vertices_probability,"new_topo_2e2g_electrons_vertices_probability/D");
  newtree->Branch("2e2g_electrons_angle",&new_topo_2e2g_electrons_angle,"new_topo_2e2g_electrons_angle/D");
  newtree->Branch("2e2g_electrons_cos_angle",&new_topo_2e2g_electrons_cos_angle,"new_topo_2e2g_electrons_cos_angle/D");
  newtree->Branch("2e2g_electron_Emin_track_length",&new_topo_2e2g_electron_Emin_track_length,"new_topo_2e2g_electron_Emin_track_length/D");
  newtree->Branch("2e2g_electron_Emax_track_length",&new_topo_2e2g_electron_Emax_track_length,"new_topo_2e2g_electron_Emax_track_length/D");

  // newtree->Branch("2e3g_electron_minimal_energy",&new_topo_2e3g_electron_minimal_energy,"new_topo_2e3g_electron_minimal_energy/D");
  // newtree->Branch("2e3g_electron_maximal_energy",&new_topo_2e3g_electron_maximal_energy,"new_topo_2e3g_electron_maximal_energy/D");
  // newtree->Branch("2e3g_gamma_min_energy",&new_topo_2e3g_gamma_min_energy,"new_topo_2e3g_gamma_min_energy/D");
  // newtree->Branch("2e3g_gamma_mid_energy",&new_topo_2e3g_gamma_mid_energy,"new_topo_2e3g_gamma_mid_energy/D");
  // newtree->Branch("2e3g_gamma_max_energy",&new_topo_2e3g_gamma_max_energy,"new_topo_2e3g_gamma_max_energy/D");
  // newtree->Branch("2e3g_electrons_energy_difference",&new_topo_2e3g_electrons_energy_difference,"new_topo_2e3g_electrons_energy_difference/D");
  // newtree->Branch("2e3g_electrons_energy_sum",&new_topo_2e3g_electrons_energy_sum,"new_topo_2e3g_electrons_energy_sum/D");
  // newtree->Branch("2e3g_electrons_gammas_energy_sum",&new_topo_2e3g_electrons_gammas_energy_sum,"new_topo_2e3g_electrons_gammas_energy_sum/D");
  // newtree->Branch("2e3g_electrons_internal_probability",&new_topo_2e3g_electrons_internal_probability,"new_topo_2e3g_electrons_internal_probability/D");
  // newtree->Branch("2e3g_electrons_external_probability",&new_topo_2e3g_electrons_external_probability,"new_topo_2e3g_electrons_external_probability/D");
  // newtree->Branch("2e3g_electron_min_gamma_min_internal_probability",&new_topo_2e3g_electron_min_gamma_min_internal_probability,"new_topo_2e3g_electron_min_gamma_min_internal_probability/D");
  // newtree->Branch("2e3g_electron_min_gamma_min_external_probability",&new_topo_2e3g_electron_min_gamma_min_external_probability,"new_topo_2e3g_electron_min_gamma_min_external_probability/D");
  // newtree->Branch("2e3g_electron_min_gamma_mid_internal_probability",&new_topo_2e3g_electron_min_gamma_mid_internal_probability,"new_topo_2e3g_electron_min_gamma_mid_internal_probability/D");
  // newtree->Branch("2e3g_electron_min_gamma_mid_external_probability",&new_topo_2e3g_electron_min_gamma_mid_external_probability,"new_topo_2e3g_electron_min_gamma_mid_external_probability/D");
  // newtree->Branch("2e3g_electron_min_gamma_max_internal_probability",&new_topo_2e3g_electron_min_gamma_max_internal_probability,"new_topo_2e3g_electron_min_gamma_max_internal_probability/D");
  // newtree->Branch("2e3g_electron_min_gamma_max_external_probability",&new_topo_2e3g_electron_min_gamma_max_external_probability,"new_topo_2e3g_electron_min_gamma_max_external_probability/D");
  // newtree->Branch("2e3g_electron_max_gamma_min_internal_probability",&new_topo_2e3g_electron_max_gamma_min_internal_probability,"new_topo_2e3g_electron_max_gamma_min_internal_probability/D");
  // newtree->Branch("2e3g_electron_max_gamma_min_external_probability",&new_topo_2e3g_electron_max_gamma_min_external_probability,"new_topo_2e3g_electron_max_gamma_min_external_probability/D");
  // newtree->Branch("2e3g_electron_max_gamma_mid_internal_probability",&new_topo_2e3g_electron_max_gamma_mid_internal_probability,"new_topo_2e3g_electron_max_gamma_mid_internal_probability/D");
  // newtree->Branch("2e3g_electron_max_gamma_mid_external_probability",&new_topo_2e3g_electron_max_gamma_mid_external_probability,"new_topo_2e3g_electron_max_gamma_mid_external_probability/D");
  // newtree->Branch("2e3g_electron_max_gamma_max_internal_probability",&new_topo_2e3g_electron_max_gamma_max_internal_probability,"new_topo_2e3g_electron_max_gamma_max_internal_probability/D");
  // newtree->Branch("2e3g_electron_max_gamma_max_external_probability",&new_topo_2e3g_electron_max_gamma_max_external_probability,"new_topo_2e3g_electron_max_gamma_max_external_probability/D");
  // newtree->Branch("2e3g_electrons_vertices_probability",&new_topo_2e3g_electrons_vertices_probability,"new_topo_2e3g_electrons_vertices_probability/D");
  // newtree->Branch("2e3g_electrons_angle",&new_topo_2e3g_electrons_angle,"new_topo_2e3g_electrons_angle/D");
  // newtree->Branch("2e3g_electrons_cos_angle",&new_topo_2e3g_electrons_cos_angle,"new_topo_2e3g_electrons_cos_angle/D");
  // newtree->Branch("2e3g_electron_Emin_track_length",&new_topo_2e3g_electron_Emin_track_length,"new_topo_2e3g_electron_Emin_track_length/D");
  // newtree->Branch("2e3g_electron_Emax_track_length",&new_topo_2e3g_electron_Emax_track_length,"new_topo_2e3g_electron_Emax_track_length/D");

  double previous_value_2e1g = 0;
  double previous_value_2e2g = 0;
  // double previous_value_2e3g = 0;

  // double previous_value_2e1g = 1e999;
  // double previous_value_2e2g = 1e999;
  // double previous_value_2e3g = 1e999;

  int count = 0;
  int nentries = oldtree->GetEntriesFast();

  for(int i = 0; i< nentries; ++i) {
    oldtree->GetEntry(i);
    if(i%10000 == 0)
      std::cout << i << std::endl;
    // std::cout << i << "  " << topo_2e1g_electron_minimal_energy << std::endl;

    // if(topo_2e2g_electron_minimal_energy < 100 && topo_2e2g_electron_minimal_energy != 0)
    //   std::cout << i << "debug  " << topo_2e2g_electron_minimal_energy << std::endl;

    //hack, by lack of a proper solution
    if(topo_2e1g_electron_minimal_energy != previous_value_2e1g &&
       topo_2e2g_electron_minimal_energy == previous_value_2e2g &&
       true) {
       // topo_2e3g_electron_minimal_energy == previous_value_2e3g) {
      // std::cout << i << "  " << topo_2e1g_electron_minimal_energy << std::endl;

      new_topo_2e1g_electron_minimal_energy = topo_2e1g_electron_minimal_energy;
      new_topo_2e1g_electron_maximal_energy = topo_2e1g_electron_maximal_energy;
      new_topo_2e1g_gamma_energy = topo_2e1g_gamma_energy;
      new_topo_2e1g_electrons_energy_difference = topo_2e1g_electrons_energy_difference;
      new_topo_2e1g_electrons_energy_sum = topo_2e1g_electrons_energy_sum;
      new_topo_2e1g_electrons_gammas_energy_sum = topo_2e1g_electrons_gammas_energy_sum;
      new_topo_2e1g_electrons_internal_probability = topo_2e1g_electrons_internal_probability;
      new_topo_2e1g_electrons_external_probability = topo_2e1g_electrons_external_probability;
      new_topo_2e1g_electron_min_gamma_internal_probability = topo_2e1g_electron_min_gamma_internal_probability;
      new_topo_2e1g_electron_min_gamma_external_probability = topo_2e1g_electron_min_gamma_external_probability;
      new_topo_2e1g_electron_max_gamma_internal_probability = topo_2e1g_electron_max_gamma_internal_probability;
      new_topo_2e1g_electron_max_gamma_external_probability = topo_2e1g_electron_max_gamma_external_probability;
      new_topo_2e1g_electrons_vertices_probability = topo_2e1g_electrons_vertices_probability;
      new_topo_2e1g_electrons_angle = topo_2e1g_electrons_angle;
      new_topo_2e1g_electrons_cos_angle = topo_2e1g_electrons_cos_angle;
      new_topo_2e1g_electron_Emin_track_length = topo_2e1g_electron_Emin_track_length;
      new_topo_2e1g_electron_Emax_track_length = topo_2e1g_electron_Emax_track_length;

      new_topo_2e2g_electron_minimal_energy = 1e999;
      new_topo_2e2g_electron_maximal_energy = 1e999;
      new_topo_2e2g_gamma_min_energy = 1e999;
      new_topo_2e2g_gamma_max_energy = 1e999;
      new_topo_2e2g_electrons_energy_difference = 1e999;
      new_topo_2e2g_electrons_energy_sum = 1e999;
      new_topo_2e2g_electrons_gammas_energy_sum = 1e999;
      new_topo_2e2g_electrons_internal_probability = 1e999;
      new_topo_2e2g_electrons_external_probability = 1e999;
      new_topo_2e2g_electron_min_gamma_min_internal_probability = 1e999;
      new_topo_2e2g_electron_min_gamma_max_internal_probability = 1e999;
      new_topo_2e2g_electron_min_gamma_min_external_probability = 1e999;
      new_topo_2e2g_electron_min_gamma_max_external_probability = 1e999;
      new_topo_2e2g_electron_max_gamma_min_internal_probability = 1e999;
      new_topo_2e2g_electron_max_gamma_max_internal_probability = 1e999;
      new_topo_2e2g_electron_max_gamma_min_external_probability = 1e999;
      new_topo_2e2g_electron_max_gamma_max_external_probability = 1e999;
      new_topo_2e2g_electrons_vertices_probability = 1e999;
      new_topo_2e2g_electrons_angle = 1e999;
      new_topo_2e2g_electrons_cos_angle = 1e999;
      new_topo_2e2g_electron_Emin_track_length = 1e999;
      new_topo_2e2g_electron_Emax_track_length = 1e999;

      // new_topo_2e3g_electron_minimal_energy = 1e999;
      // new_topo_2e3g_electron_maximal_energy = 1e999;
      // new_topo_2e3g_gamma_min_energy = 1e999;
      // new_topo_2e3g_gamma_mid_energy = 1e999;
      // new_topo_2e3g_gamma_max_energy = 1e999;
      // new_topo_2e3g_electrons_energy_difference = 1e999;
      // new_topo_2e3g_electrons_energy_sum = 1e999;
      // new_topo_2e3g_electrons_gammas_energy_sum = 1e999;
      // new_topo_2e3g_electrons_internal_probability = 1e999;
      // new_topo_2e3g_electrons_external_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_min_internal_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_min_external_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_mid_internal_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_mid_external_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_max_internal_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_max_external_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_min_internal_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_min_external_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_mid_internal_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_mid_external_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_max_internal_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_max_external_probability = 1e999;
      // new_topo_2e3g_electrons_vertices_probability = 1e999;
      // new_topo_2e3g_electrons_angle = 1e999;
      // new_topo_2e3g_electrons_cos_angle = 1e999;
      // new_topo_2e3g_electron_Emin_track_length = 1e999;
      // new_topo_2e3g_electron_Emax_track_length = 1e999;

    }
    else if(topo_2e1g_electron_minimal_energy == previous_value_2e1g &&
            topo_2e2g_electron_minimal_energy != previous_value_2e2g &&
            true) {
      // topo_2e3g_electron_minimal_energy == previous_value_2e3g) {

      new_topo_2e1g_electron_minimal_energy = 1e999;
      new_topo_2e1g_electron_maximal_energy = 1e999;
      new_topo_2e1g_gamma_energy = 1e999;
      new_topo_2e1g_electrons_energy_difference = 1e999;
      new_topo_2e1g_electrons_energy_sum = 1e999;
      new_topo_2e1g_electrons_gammas_energy_sum = 1e999;
      new_topo_2e1g_electrons_internal_probability = 1e999;
      new_topo_2e1g_electrons_external_probability = 1e999;
      new_topo_2e1g_electron_min_gamma_internal_probability = 1e999;
      new_topo_2e1g_electron_min_gamma_external_probability = 1e999;
      new_topo_2e1g_electron_max_gamma_internal_probability = 1e999;
      new_topo_2e1g_electron_max_gamma_external_probability = 1e999;
      new_topo_2e1g_electrons_vertices_probability = 1e999;
      new_topo_2e1g_electrons_angle = 1e999;
      new_topo_2e1g_electrons_cos_angle = 1e999;
      new_topo_2e1g_electron_Emin_track_length = 1e999;
      new_topo_2e1g_electron_Emax_track_length = 1e999;

      new_topo_2e2g_electron_minimal_energy = topo_2e2g_electron_minimal_energy;
      new_topo_2e2g_electron_maximal_energy = topo_2e2g_electron_maximal_energy;
      new_topo_2e2g_gamma_min_energy = topo_2e2g_gamma_min_energy;
      new_topo_2e2g_gamma_max_energy = topo_2e2g_gamma_max_energy;
      new_topo_2e2g_electrons_energy_difference = topo_2e2g_electrons_energy_difference;
      new_topo_2e2g_electrons_energy_sum = topo_2e2g_electrons_energy_sum;
      new_topo_2e2g_electrons_gammas_energy_sum = topo_2e2g_electrons_gammas_energy_sum;
      new_topo_2e2g_electrons_internal_probability = topo_2e2g_electrons_internal_probability;
      new_topo_2e2g_electrons_external_probability = topo_2e2g_electrons_external_probability;
      new_topo_2e2g_electron_min_gamma_min_internal_probability = topo_2e2g_electron_min_gamma_min_internal_probability;
      new_topo_2e2g_electron_min_gamma_max_internal_probability = topo_2e2g_electron_min_gamma_max_internal_probability;
      new_topo_2e2g_electron_min_gamma_min_external_probability = topo_2e2g_electron_min_gamma_min_external_probability;
      new_topo_2e2g_electron_min_gamma_max_external_probability = topo_2e2g_electron_min_gamma_max_external_probability;
      new_topo_2e2g_electron_max_gamma_min_internal_probability = topo_2e2g_electron_max_gamma_min_internal_probability;
      new_topo_2e2g_electron_max_gamma_max_internal_probability = topo_2e2g_electron_max_gamma_max_internal_probability;
      new_topo_2e2g_electron_max_gamma_min_external_probability = topo_2e2g_electron_max_gamma_min_external_probability;
      new_topo_2e2g_electron_max_gamma_max_external_probability = topo_2e2g_electron_max_gamma_max_external_probability;
      new_topo_2e2g_electrons_vertices_probability = topo_2e2g_electrons_vertices_probability;
      new_topo_2e2g_electrons_angle = topo_2e2g_electrons_angle;
      new_topo_2e2g_electrons_cos_angle = topo_2e2g_electrons_cos_angle;
      new_topo_2e2g_electron_Emin_track_length = topo_2e2g_electron_Emin_track_length;
      new_topo_2e2g_electron_Emax_track_length = topo_2e2g_electron_Emax_track_length;

      // new_topo_2e3g_electron_minimal_energy = 1e999;
      // new_topo_2e3g_electron_maximal_energy = 1e999;
      // new_topo_2e3g_gamma_min_energy = 1e999;
      // new_topo_2e3g_gamma_mid_energy = 1e999;
      // new_topo_2e3g_gamma_max_energy = 1e999;
      // new_topo_2e3g_electrons_energy_difference = 1e999;
      // new_topo_2e3g_electrons_energy_sum = 1e999;
      // new_topo_2e3g_electrons_gammas_energy_sum = 1e999;
      // new_topo_2e3g_electrons_internal_probability = 1e999;
      // new_topo_2e3g_electrons_external_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_min_internal_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_min_external_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_mid_internal_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_mid_external_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_max_internal_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_max_external_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_min_internal_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_min_external_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_mid_internal_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_mid_external_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_max_internal_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_max_external_probability = 1e999;
      // new_topo_2e3g_electrons_vertices_probability = 1e999;
      // new_topo_2e3g_electrons_angle = 1e999;
      // new_topo_2e3g_electrons_cos_angle = 1e999;
      // new_topo_2e3g_electron_Emin_track_length = 1e999;
      // new_topo_2e3g_electron_Emax_track_length = 1e999;

    }
    else if(topo_2e1g_electron_minimal_energy == previous_value_2e1g &&
            topo_2e2g_electron_minimal_energy == previous_value_2e2g &&
            true) {
            // topo_2e3g_electron_minimal_energy != previous_value_2e3g) {

      new_topo_2e1g_electron_minimal_energy = 1e999;
      new_topo_2e1g_electron_maximal_energy = 1e999;
      new_topo_2e1g_gamma_energy = 1e999;
      new_topo_2e1g_electrons_energy_difference = 1e999;
      new_topo_2e1g_electrons_energy_sum = 1e999;
      new_topo_2e1g_electrons_gammas_energy_sum = 1e999;
      new_topo_2e1g_electrons_internal_probability = 1e999;
      new_topo_2e1g_electrons_external_probability = 1e999;
      new_topo_2e1g_electron_min_gamma_internal_probability = 1e999;
      new_topo_2e1g_electron_min_gamma_external_probability = 1e999;
      new_topo_2e1g_electron_max_gamma_internal_probability = 1e999;
      new_topo_2e1g_electron_max_gamma_external_probability = 1e999;
      new_topo_2e1g_electrons_vertices_probability = 1e999;
      new_topo_2e1g_electrons_angle = 1e999;
      new_topo_2e1g_electrons_cos_angle = 1e999;
      new_topo_2e1g_electron_Emin_track_length = 1e999;
      new_topo_2e1g_electron_Emax_track_length = 1e999;

      new_topo_2e2g_electron_minimal_energy = 1e999;
      new_topo_2e2g_electron_maximal_energy = 1e999;
      new_topo_2e2g_gamma_min_energy = 1e999;
      new_topo_2e2g_gamma_max_energy = 1e999;
      new_topo_2e2g_electrons_energy_difference = 1e999;
      new_topo_2e2g_electrons_energy_sum = 1e999;
      new_topo_2e2g_electrons_gammas_energy_sum = 1e999;
      new_topo_2e2g_electrons_internal_probability = 1e999;
      new_topo_2e2g_electrons_external_probability = 1e999;
      new_topo_2e2g_electron_min_gamma_min_internal_probability = 1e999;
      new_topo_2e2g_electron_min_gamma_max_internal_probability = 1e999;
      new_topo_2e2g_electron_min_gamma_min_external_probability = 1e999;
      new_topo_2e2g_electron_min_gamma_max_external_probability = 1e999;
      new_topo_2e2g_electron_max_gamma_min_internal_probability = 1e999;
      new_topo_2e2g_electron_max_gamma_max_internal_probability = 1e999;
      new_topo_2e2g_electron_max_gamma_min_external_probability = 1e999;
      new_topo_2e2g_electron_max_gamma_max_external_probability = 1e999;
      new_topo_2e2g_electrons_vertices_probability = 1e999;
      new_topo_2e2g_electrons_angle = 1e999;
      new_topo_2e2g_electrons_cos_angle = 1e999;
      new_topo_2e2g_electron_Emin_track_length = 1e999;
      new_topo_2e2g_electron_Emax_track_length = 1e999;

      // new_topo_2e3g_electron_minimal_energy = topo_2e3g_electron_minimal_energy;
      // new_topo_2e3g_electron_maximal_energy = topo_2e3g_electron_maximal_energy;
      // new_topo_2e3g_gamma_min_energy = topo_2e3g_gamma_min_energy;
      // new_topo_2e3g_gamma_mid_energy = topo_2e3g_gamma_mid_energy;
      // new_topo_2e3g_gamma_max_energy = topo_2e3g_gamma_max_energy;
      // new_topo_2e3g_electrons_energy_difference = topo_2e3g_electrons_energy_difference;
      // new_topo_2e3g_electrons_energy_sum = topo_2e3g_electrons_energy_sum;
      // new_topo_2e3g_electrons_gammas_energy_sum = topo_2e3g_electrons_gammas_energy_sum;
      // new_topo_2e3g_electrons_internal_probability = topo_2e3g_electrons_internal_probability;
      // new_topo_2e3g_electrons_external_probability = topo_2e3g_electrons_external_probability;
      // new_topo_2e3g_electron_min_gamma_min_internal_probability = topo_2e3g_electron_min_gamma_min_internal_probability;
      // new_topo_2e3g_electron_min_gamma_min_external_probability = topo_2e3g_electron_min_gamma_min_external_probability;
      // new_topo_2e3g_electron_min_gamma_mid_internal_probability = topo_2e3g_electron_min_gamma_mid_internal_probability;
      // new_topo_2e3g_electron_min_gamma_mid_external_probability = topo_2e3g_electron_min_gamma_mid_external_probability;
      // new_topo_2e3g_electron_min_gamma_max_internal_probability = topo_2e3g_electron_min_gamma_max_internal_probability;
      // new_topo_2e3g_electron_min_gamma_max_external_probability = topo_2e3g_electron_min_gamma_max_external_probability;
      // new_topo_2e3g_electron_max_gamma_min_internal_probability = topo_2e3g_electron_max_gamma_min_internal_probability;
      // new_topo_2e3g_electron_max_gamma_min_external_probability = topo_2e3g_electron_max_gamma_min_external_probability;
      // new_topo_2e3g_electron_max_gamma_mid_internal_probability = topo_2e3g_electron_max_gamma_mid_internal_probability;
      // new_topo_2e3g_electron_max_gamma_mid_external_probability = topo_2e3g_electron_max_gamma_mid_external_probability;
      // new_topo_2e3g_electron_max_gamma_max_internal_probability = topo_2e3g_electron_max_gamma_max_internal_probability;
      // new_topo_2e3g_electron_max_gamma_max_external_probability = topo_2e3g_electron_max_gamma_max_external_probability;
      // new_topo_2e3g_electrons_vertices_probability = topo_2e3g_electrons_vertices_probability;
      // new_topo_2e3g_electrons_angle = topo_2e3g_electrons_angle;
      // new_topo_2e3g_electrons_cos_angle = topo_2e3g_electrons_cos_angle;
      // new_topo_2e3g_electron_Emin_track_length = topo_2e3g_electron_Emin_track_length;
      // new_topo_2e3g_electron_Emax_track_length = topo_2e3g_electron_Emax_track_length;

    }
    else if(topo_2e1g_electron_minimal_energy == previous_value_2e1g &&
            topo_2e2g_electron_minimal_energy == previous_value_2e2g &&
            true) {
      // topo_2e3g_electron_minimal_energy == previous_value_2e3g) {

      new_topo_2e1g_electron_minimal_energy = 1e999;
      new_topo_2e1g_electron_maximal_energy = 1e999;
      new_topo_2e1g_gamma_energy = 1e999;
      new_topo_2e1g_electrons_energy_difference = 1e999;
      new_topo_2e1g_electrons_energy_sum = 1e999;
      new_topo_2e1g_electrons_gammas_energy_sum = 1e999;
      new_topo_2e1g_electrons_internal_probability = 1e999;
      new_topo_2e1g_electrons_external_probability = 1e999;
      new_topo_2e1g_electron_min_gamma_internal_probability = 1e999;
      new_topo_2e1g_electron_min_gamma_external_probability = 1e999;
      new_topo_2e1g_electron_max_gamma_internal_probability = 1e999;
      new_topo_2e1g_electron_max_gamma_external_probability = 1e999;
      new_topo_2e1g_electrons_vertices_probability = 1e999;
      new_topo_2e1g_electrons_angle = 1e999;
      new_topo_2e1g_electrons_cos_angle = 1e999;
      new_topo_2e1g_electron_Emin_track_length = 1e999;
      new_topo_2e1g_electron_Emax_track_length = 1e999;

      new_topo_2e2g_electron_minimal_energy = 1e999;
      new_topo_2e2g_electron_maximal_energy = 1e999;
      new_topo_2e2g_gamma_min_energy = 1e999;
      new_topo_2e2g_gamma_max_energy = 1e999;
      new_topo_2e2g_electrons_energy_difference = 1e999;
      new_topo_2e2g_electrons_energy_sum = 1e999;
      new_topo_2e2g_electrons_gammas_energy_sum = 1e999;
      new_topo_2e2g_electrons_internal_probability = 1e999;
      new_topo_2e2g_electrons_external_probability = 1e999;
      new_topo_2e2g_electron_min_gamma_min_internal_probability = 1e999;
      new_topo_2e2g_electron_min_gamma_max_internal_probability = 1e999;
      new_topo_2e2g_electron_min_gamma_min_external_probability = 1e999;
      new_topo_2e2g_electron_min_gamma_max_external_probability = 1e999;
      new_topo_2e2g_electron_max_gamma_min_internal_probability = 1e999;
      new_topo_2e2g_electron_max_gamma_max_internal_probability = 1e999;
      new_topo_2e2g_electron_max_gamma_min_external_probability = 1e999;
      new_topo_2e2g_electron_max_gamma_max_external_probability = 1e999;
      new_topo_2e2g_electrons_vertices_probability = 1e999;
      new_topo_2e2g_electrons_angle = 1e999;
      new_topo_2e2g_electrons_cos_angle = 1e999;
      new_topo_2e2g_electron_Emin_track_length = 1e999;
      new_topo_2e2g_electron_Emax_track_length = 1e999;

      // new_topo_2e3g_electron_minimal_energy = 1e999;
      // new_topo_2e3g_electron_maximal_energy = 1e999;
      // new_topo_2e3g_gamma_min_energy = 1e999;
      // new_topo_2e3g_gamma_mid_energy = 1e999;
      // new_topo_2e3g_gamma_max_energy = 1e999;
      // new_topo_2e3g_electrons_energy_difference = 1e999;
      // new_topo_2e3g_electrons_energy_sum = 1e999;
      // new_topo_2e3g_electrons_gammas_energy_sum = 1e999;
      // new_topo_2e3g_electrons_internal_probability = 1e999;
      // new_topo_2e3g_electrons_external_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_min_internal_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_min_external_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_mid_internal_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_mid_external_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_max_internal_probability = 1e999;
      // new_topo_2e3g_electron_min_gamma_max_external_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_min_internal_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_min_external_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_mid_internal_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_mid_external_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_max_internal_probability = 1e999;
      // new_topo_2e3g_electron_max_gamma_max_external_probability = 1e999;
      // new_topo_2e3g_electrons_vertices_probability = 1e999;
      // new_topo_2e3g_electrons_angle = 1e999;
      // new_topo_2e3g_electrons_cos_angle = 1e999;
      // new_topo_2e3g_electron_Emin_track_length = 1e999;
      // new_topo_2e3g_electron_Emax_track_length = 1e999;

    }
    // else {
    //   if(i=!0) {
    //     std::cout << " Should not happen... "  << i << std::endl;
    //     std::cout << "            2e1g  " << topo_2e1g_electron_minimal_energy << " and " << previous_value_2e1g << std::endl;
    //     std::cout << "            2e2g  " << topo_2e2g_electron_minimal_energy << " and " << previous_value_2e2g << std::endl;
    //   // std::cout << "            2e3g  " << topo_2e3g_electron_minimal_energy << " and " << previous_value_2e3g << std::endl;
    //     std::cout << count << std::endl;
    //     return;
    //   }
    // }

    TBranch * b_2e1g_electron_minimal_energy = newtree->GetBranch("2e1g_electron_minimal_energy");
    TBranch * b_2e1g_electron_maximal_energy = newtree->GetBranch("2e1g_electron_maximal_energy");
    TBranch * b_2e1g_gamma_energy = newtree->GetBranch("2e1g_gamma_energy");
    TBranch * b_2e1g_electrons_energy_difference = newtree->GetBranch("2e1g_electrons_energy_difference");
    TBranch * b_2e1g_electrons_energy_sum = newtree->GetBranch("2e1g_electrons_energy_sum");
    TBranch * b_2e1g_electrons_gammas_energy_sum = newtree->GetBranch("2e1g_electrons_gammas_energy_sum");
    TBranch * b_2e1g_electrons_internal_probability = newtree->GetBranch("2e1g_electrons_internal_probability");
    TBranch * b_2e1g_electrons_external_probability = newtree->GetBranch("2e1g_electrons_external_probability");
    TBranch * b_2e1g_electron_min_gamma_internal_probability = newtree->GetBranch("2e1g_electron_min_gamma_internal_probability");
    TBranch * b_2e1g_electron_min_gamma_external_probability = newtree->GetBranch("2e1g_electron_min_gamma_external_probability");
    TBranch * b_2e1g_electron_max_gamma_internal_probability = newtree->GetBranch("2e1g_electron_max_gamma_internal_probability");
    TBranch * b_2e1g_electron_max_gamma_external_probability = newtree->GetBranch("2e1g_electron_max_gamma_external_probability");
    TBranch * b_2e1g_electrons_vertices_probability = newtree->GetBranch("2e1g_electrons_vertices_probability");
    TBranch * b_2e1g_electrons_angle = newtree->GetBranch("2e1g_electrons_angle");
    TBranch * b_2e1g_electrons_cos_angle = newtree->GetBranch("2e1g_electrons_cos_angle");
    TBranch * b_2e1g_electron_Emin_track_length = newtree->GetBranch("2e1g_electron_Emin_track_length");
    TBranch * b_2e1g_electron_Emax_track_length = newtree->GetBranch("2e1g_electron_Emax_track_length");

    TBranch * b_2e2g_electron_minimal_energy = newtree->GetBranch("2e2g_electron_minimal_energy");
    TBranch * b_2e2g_electron_maximal_energy = newtree->GetBranch("2e2g_electron_maximal_energy");
    TBranch * b_2e2g_gamma_min_energy = newtree->GetBranch("2e2g_gamma_min_energy");
    TBranch * b_2e2g_gamma_max_energy = newtree->GetBranch("2e2g_gamma_max_energy");
    TBranch * b_2e2g_electrons_energy_difference = newtree->GetBranch("2e2g_electrons_energy_difference");
    TBranch * b_2e2g_electrons_energy_sum = newtree->GetBranch("2e2g_electrons_energy_sum");
    TBranch * b_2e2g_electrons_gammas_energy_sum = newtree->GetBranch("2e2g_electrons_gammas_energy_sum");
    TBranch * b_2e2g_electrons_internal_probability = newtree->GetBranch("2e2g_electrons_internal_probability");
    TBranch * b_2e2g_electrons_external_probability = newtree->GetBranch("2e2g_electrons_external_probability");
    TBranch * b_2e2g_electron_min_gamma_min_internal_probability = newtree->GetBranch("2e2g_electron_min_gamma_min_internal_probability");
    TBranch * b_2e2g_electron_min_gamma_max_internal_probability = newtree->GetBranch("2e2g_electron_min_gamma_max_internal_probability");
    TBranch * b_2e2g_electron_min_gamma_min_external_probability = newtree->GetBranch("2e2g_electron_min_gamma_min_external_probability");
    TBranch * b_2e2g_electron_min_gamma_max_external_probability = newtree->GetBranch("2e2g_electron_min_gamma_max_external_probability");
    TBranch * b_2e2g_electron_max_gamma_min_internal_probability = newtree->GetBranch("2e2g_electron_max_gamma_min_internal_probability");
    TBranch * b_2e2g_electron_max_gamma_max_internal_probability = newtree->GetBranch("2e2g_electron_max_gamma_max_internal_probability");
    TBranch * b_2e2g_electron_max_gamma_min_external_probability = newtree->GetBranch("2e2g_electron_max_gamma_min_external_probability");
    TBranch * b_2e2g_electron_max_gamma_max_external_probability = newtree->GetBranch("2e2g_electron_max_gamma_max_external_probability");
    TBranch * b_2e2g_electrons_vertices_probability = newtree->GetBranch("2e2g_electrons_vertices_probability");;
    TBranch * b_2e2g_electrons_angle = newtree->GetBranch("2e2g_electrons_angle");
    TBranch * b_2e2g_electrons_cos_angle = newtree->GetBranch("2e2g_electrons_cos_angle");
    TBranch * b_2e2g_electron_Emin_track_length = newtree->GetBranch("2e2g_electron_Emin_track_length");
    TBranch * b_2e2g_electron_Emax_track_length = newtree->GetBranch("2e2g_electron_Emax_track_length");

    // TBranch * b_2e3g_electron_minimal_energy = newtree->GetBranch("2e3g_electron_minimal_energy");
    // TBranch * b_2e3g_electron_maximal_energy = newtree->GetBranch("2e3g_electron_maximal_energy");
    // TBranch * b_2e3g_gamma_min_energy = newtree->GetBranch("2e3g_gamma_min_energy");
    // TBranch * b_2e3g_gamma_mid_energy = newtree->GetBranch("2e3g_gamma_mid_energy");
    // TBranch * b_2e3g_gamma_max_energy = newtree->GetBranch("2e3g_gamma_max_energy");
    // TBranch * b_2e3g_electrons_energy_difference = newtree->GetBranch("2e3g_electrons_energy_difference");
    // TBranch * b_2e3g_electrons_energy_sum = newtree->GetBranch("2e3g_electrons_energy_sum");
    // TBranch * b_2e3g_electrons_gammas_energy_sum = newtree->GetBranch("2e3g_electrons_gammas_energy_sum");
    // TBranch * b_2e3g_electrons_internal_probability = newtree->GetBranch("2e3g_electrons_internal_probability");
    // TBranch * b_2e3g_electrons_external_probability = newtree->GetBranch("2e3g_electrons_external_probability");
    // TBranch * b_2e3g_electron_min_gamma_min_internal_probability = newtree->GetBranch("2e3g_electron_min_gamma_min_internal_probability");
    // TBranch * b_2e3g_electron_min_gamma_min_external_probability = newtree->GetBranch("2e3g_electron_min_gamma_min_external_probability");
    // TBranch * b_2e3g_electron_min_gamma_mid_internal_probability = newtree->GetBranch("2e3g_electron_min_gamma_mid_internal_probability");
    // TBranch * b_2e3g_electron_min_gamma_mid_external_probability = newtree->GetBranch("2e3g_electron_min_gamma_mid_external_probability");
    // TBranch * b_2e3g_electron_min_gamma_max_internal_probability = newtree->GetBranch("2e3g_electron_min_gamma_max_internal_probability");
    // TBranch * b_2e3g_electron_min_gamma_max_external_probability = newtree->GetBranch("2e3g_electron_min_gamma_max_external_probability");
    // TBranch * b_2e3g_electron_max_gamma_min_internal_probability = newtree->GetBranch("2e3g_electron_max_gamma_min_internal_probability");
    // TBranch * b_2e3g_electron_max_gamma_min_external_probability = newtree->GetBranch("2e3g_electron_max_gamma_min_external_probability");
    // TBranch * b_2e3g_electron_max_gamma_mid_internal_probability = newtree->GetBranch("2e3g_electron_max_gamma_mid_internal_probability");
    // TBranch * b_2e3g_electron_max_gamma_mid_external_probability = newtree->GetBranch("2e3g_electron_max_gamma_mid_external_probability");
    // TBranch * b_2e3g_electron_max_gamma_max_internal_probability = newtree->GetBranch("2e3g_electron_max_gamma_max_internal_probability");
    // TBranch * b_2e3g_electron_max_gamma_max_external_probability = newtree->GetBranch("2e3g_electron_max_gamma_max_external_probability");
    // TBranch * b_2e3g_electrons_vertices_probability = newtree->GetBranch("2e3g_electrons_vertices_probability");
    // TBranch * b_2e3g_electrons_angle = newtree->GetBranch("2e3g_electrons_angle");
    // TBranch * b_2e3g_electrons_cos_angle = newtree->GetBranch("2e3g_electrons_cos_angle");
    // TBranch * b_2e3g_electron_Emin_track_length = newtree->GetBranch("2e3g_electron_Emin_track_length");
    // TBranch * b_2e3g_electron_Emax_track_length = newtree->GetBranch("2e3g_electron_Emax_track_length");

    b_2e1g_electron_minimal_energy->Fill();
    b_2e1g_electron_maximal_energy->Fill();
    b_2e1g_gamma_energy->Fill();
    b_2e1g_electrons_energy_difference->Fill();
    b_2e1g_electrons_energy_sum->Fill();
    b_2e1g_electrons_gammas_energy_sum->Fill();
    b_2e1g_electrons_internal_probability->Fill();
    b_2e1g_electrons_external_probability->Fill();
    b_2e1g_electron_min_gamma_internal_probability->Fill();
    b_2e1g_electron_min_gamma_external_probability->Fill();
    b_2e1g_electron_max_gamma_internal_probability->Fill();
    b_2e1g_electron_max_gamma_external_probability->Fill();
    b_2e1g_electrons_vertices_probability->Fill();
    b_2e1g_electrons_angle->Fill();
    b_2e1g_electrons_cos_angle->Fill();
    b_2e1g_electron_Emin_track_length->Fill();
    b_2e1g_electron_Emax_track_length->Fill();

    b_2e2g_electron_minimal_energy->Fill();
    b_2e2g_electron_maximal_energy->Fill();
    b_2e2g_gamma_min_energy->Fill();
    b_2e2g_gamma_max_energy->Fill();
    b_2e2g_electrons_energy_difference->Fill();
    b_2e2g_electrons_energy_sum->Fill();
    b_2e2g_electrons_gammas_energy_sum->Fill();
    b_2e2g_electrons_internal_probability->Fill();
    b_2e2g_electrons_external_probability->Fill();
    b_2e2g_electron_min_gamma_min_internal_probability->Fill();
    b_2e2g_electron_min_gamma_max_internal_probability->Fill();
    b_2e2g_electron_min_gamma_min_external_probability->Fill();
    b_2e2g_electron_min_gamma_max_external_probability->Fill();
    b_2e2g_electron_max_gamma_min_internal_probability->Fill();
    b_2e2g_electron_max_gamma_max_internal_probability->Fill();
    b_2e2g_electron_max_gamma_min_external_probability->Fill();
    b_2e2g_electron_max_gamma_max_external_probability->Fill();
    b_2e2g_electrons_vertices_probability->Fill();
    b_2e2g_electrons_angle->Fill();
    b_2e2g_electrons_cos_angle->Fill();
    b_2e2g_electron_Emin_track_length->Fill();
    b_2e2g_electron_Emax_track_length->Fill();

    // b_2e3g_electron_minimal_energy->Fill();
    // b_2e3g_electron_maximal_energy->Fill();
    // b_2e3g_gamma_min_energy->Fill();
    // b_2e3g_gamma_mid_energy->Fill();
    // b_2e3g_gamma_max_energy->Fill();
    // b_2e3g_electrons_energy_difference->Fill();
    // b_2e3g_electrons_energy_sum->Fill();
    // b_2e3g_electrons_gammas_energy_sum->Fill();
    // b_2e3g_electrons_internal_probability->Fill();
    // b_2e3g_electrons_external_probability->Fill();
    // b_2e3g_electron_min_gamma_min_internal_probability->Fill();
    // b_2e3g_electron_min_gamma_min_external_probability->Fill();
    // b_2e3g_electron_min_gamma_mid_internal_probability->Fill();
    // b_2e3g_electron_min_gamma_mid_external_probability->Fill();
    // b_2e3g_electron_min_gamma_max_internal_probability->Fill();
    // b_2e3g_electron_min_gamma_max_external_probability->Fill();
    // b_2e3g_electron_max_gamma_min_internal_probability->Fill();
    // b_2e3g_electron_max_gamma_min_external_probability->Fill();
    // b_2e3g_electron_max_gamma_mid_internal_probability->Fill();
    // b_2e3g_electron_max_gamma_mid_external_probability->Fill();
    // b_2e3g_electron_max_gamma_max_internal_probability->Fill();
    // b_2e3g_electron_max_gamma_max_external_probability->Fill();
    // b_2e3g_electrons_vertices_probability->Fill();
    // b_2e3g_electrons_angle->Fill();
    // b_2e3g_electrons_cos_angle->Fill();
    // b_2e3g_electron_Emin_track_length->Fill();
    // b_2e3g_electron_Emax_track_length->Fill();

    // newtree->Fill(); // not the tree otherwise twice the events
    previous_value_2e1g = topo_2e1g_electron_minimal_energy;
    previous_value_2e2g = topo_2e2g_electron_minimal_energy;
    // previous_value_2e3g = topo_2e3g_electron_minimal_energy;
    count++;
  }

  // TFile *newfile = new TFile("test_new_tree.root","recreate");
  // TTree *newtree = oldtree->CloneTree();

  // newtree->Print();
  newfile->Write();
  delete f;
  delete newfile;

  return;
}
