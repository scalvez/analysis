#include <iostream>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
// #include <gSystem.h>
#include <map>
#include "analysis_config.h"
#include "channel_selection.h"
// #include "sensitivity_measurements.h"

extern std::map < TString , double > quantity_efficiency;
// extern std::map < TString , TH1F* > quantity_pdf;

void channel_selection(TString isotope, std::vector<TString> quantities_pdf, bool normalize)
{
  // const std::string sep = "_";

  // const std::string topology_1e   = "1e";
  // const std::string topology_1e1a = "1e1a";
  // const std::string topology_2e   = "2e";
  // const std::string topology_1e1p = "1e1p";
  // const std::string topology_2p   = "2p";
  // const std::string topology_1e1g = "1e1g";
  // const std::string topology_1e2g = "1e2g";
  // const std::string topology_1e3g = "1e3g";
  // const std::string topology_2e1g = "2e1g";
  // const std::string topology_2e2g = "2e2g";
  // const std::string topology_2e3g = "2e3g";

  // TString input_file = "../data/trees/" + isotope + "/merge.root";
  // TString output_file = "../data/trees/" + isotope + "_pdf.root";

  TString input_file = "$SW_WORK_DIR/test_data_selection/" + isotope + "/root/merge.root";
  TString output_file = "/tmp/calvez/snemo.d/test_data_selection/test_pdf.root";

  TFile *f = TFile::Open(input_file);
  TTree *tree = (TTree*)f->Get("snemodata");
  TFile *f_output= new TFile(output_file, "RECREATE");

  double isotope_mc_size = get_isotope_mc_size(input_file);

  // TCut cut = get_channel_cut(channel);
  for(unsigned int j=0; j<quantities_pdf.size(); ++j) {

    int nbins;
    double xmin, xmax;
    //default values
    nbins = 100;
    xmin = 0;
    xmax = 5;

    TString qty = quantities_pdf.at(j);

    get_histogram_options(qty, nbins, xmin, xmax);

    TH1F* h = new TH1F(qty,qty,nbins,xmin,xmax);

    // //Remove un-initialized events (MC sim before May 26th 16)
    // //this, or concatenate char and convert to string
    // if(qty.Contains("2e1g")) {
    //   TCut cut = "2e1g_electrons_gammas_energy_sum != 0";
    //   tree->Project(qty,qty,cut);
    // }
    // else if(qty.Contains("2e2g")) {
    //   TCut cut = "2e2g_electrons_gammas_energy_sum != 0";
    //   tree->Project(qty,qty,cut);
    // }
    // else if(qty.Contains("2e3g")) {
    //   TCut cut = "2e3g_electrons_gammas_energy_sum != 0";
    //   tree->Project(qty,qty,cut);
    // }
    // else
    //   tree->Project(qty,qty);

    // //No cuts, projected into h via its name
    // tree->Project(qty,qty);

    // TCut cut = "2e_electrons_vertex_location == 0:2e_electrons_vertex_location == 0";
    // TCut cut = "2e_vertices_distance_y <= 0";

    // tree->Project(qty,qty,"2e_electrons_internal_probability >= 0.04 && 2e_electrons_vertices_distance_y <= 60 && 2e_electrons_vertices_distance_z <= 70 && 2e_electrons_vertex_location == 0 && 2e_electrons_energy_sum >= 2");

    //get rid of 3.3 MeV event
    if(isotope.Contains("2nu"))
      tree->Project(qty,qty,"2e_electrons_internal_probability >= 0.04 && 2e_electrons_vertices_distance_y <= 60 && 2e_electrons_vertices_distance_z <= 70 && 2e_electrons_vertex_location == 0 && 2e_electrons_energy_sum <= 3");
    else
    tree->Project(qty,qty,"2e_electrons_internal_probability >= 0.04 && 2e_electrons_vertices_distance_y <= 60 && 2e_electrons_vertices_distance_z <= 70 && 2e_electrons_vertex_location == 0");

    // tree->Project(qty,qty,"2e_electrons_internal_probability >= 0.0");

    // h->ClearUnderflowAndOverflow(); Unavailable with this version of ROOT
    h->SetBinContent(0,0);
    h->SetBinContent(nbins+1,0);

    // h->SetDirectory(0);
    // Don't care about the weights for now

    // tree->Project("h","1e1g_electron_gamma_energy_sum","","",1000);

    TString isotope_quantity = isotope + "_" + qty;
    // std::cout << "inserting " << isotope_quantity  <<  "  " <<  h->Integral(1,h->GetXaxis()->GetNbins())/isotope_mc_size << std::endl;
    if(isotope.Contains("2nu"))
      quantity_efficiency.insert(std::pair<TString,double>(isotope_quantity,h->Integral(1,h->GetXaxis()->GetNbins())/isotope_mc_size/25.));
    else
    quantity_efficiency.insert(std::pair<TString,double>(isotope_quantity,h->Integral(1,h->GetXaxis()->GetNbins())/isotope_mc_size));

    h->Sumw2();

    // if(normalize)
    //   h->Scale(1./h->Integral(1,h->GetXaxis()->GetNbins()));
    h->SetDrawOption("APL");
    h->Write();

    // std::cout << "Inserting pair " << std::endl;
    // // quantity_pdf.insert(std::pair<TString,TH1F*>(isotope_quantity,h));
    // // quantity_pdf[isotope_quantity] = h;
    // std::cout << "Inserted pair " << std::endl;
    // std::cout << "Checking histo reading " << std::endl;
    // std::cout << "               " << isotope_ <<  "    " << quantity_pdf.at(isotope_quantity)->GetBinContent(20) << std::endl;

  }

  f->Close();
  f_output->Close();

  return;
}
