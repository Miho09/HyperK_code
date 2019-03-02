

#include <stdio.h>
#include <stdlib.h>
#include <TH3.h>
#include <TH1.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
using std::vector;
using std::sort;



void reset_test_PMT(char *filename=NULL) {
  /* A simple script to plot aspects of phototube hits
   * This code is rather cavalier; I should be checking return values, etc.
   * First revision 6-24-10 David Webber
   *
   * I like to run this macro as
   * $ root -l -x 'read_PMT.C("../wcsim.root")'
   */



  gROOT->Reset();
  char* wcsimdirenv;
  wcsimdirenv = getenv ("WCSIMDIR");
  if(wcsimdirenv !=  NULL){
    gSystem->Load("${WCSIMDIR}/libWCSimRoot.so");
  }else{
    gSystem->Load("../libWCSimRoot.so");
  }
  gStyle->SetOptStat(1);

  TFile *f;
  if (filename==NULL){
    f = new TFile("../wcsim.root");
  }else{
    f = new TFile(filename);
  }
  if (!f->IsOpen()){
    cout << "Error, could not open input file: " << filename << endl;
    return -1;
  }

  TTree  *wcsimT = f->Get("wcsimT");

  WCSimRootEvent *wcsimroothyperevent = new WCSimRootEvent();
  wcsimT->SetBranchAddress("wcsimrootevent",&wcsimroothyperevent);

  TTree  *wcsimGeoT = (TTree*) f->Get("wcsimGeoT");

  WCSimRootGeom *wcsimrootgeom = 0;
  wcsimGeoT->SetBranchAddress("wcsimrootgeom",&wcsimrootgeom);
  // cout << "wcsimrootgeom value: " << wcsimrootgeom << endl;
  // cout << "getentry: " << wcsimGeoT->GetEntries() << endl;
  wcsimGeoT->GetEntry(0);

  // Force deletion to prevent memory leak when issuing multiple
  // calls to GetEvent()
  wcsimT->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

  //--------------------------
  // As you can see, there are lots of ways to get the number of hits.
  cout << "Nb of entries " << wcsimT->GetEntries() << endl;

  //-----------------------


  TH2F *histo_2 = new TH2F("Uncertainty", "Chi Squared",100,0,300000, 100, 0, 50);
  histo_2->SetXTitle("Chi Squared");
  histo_2->SetYTitle("Average");


  const long unsigned int nbEntries = wcsimT->GetEntries();

  vector <double> arr(0);
  vector<vector<double> >  hist(110,arr);
  vector<vector<vector<double> > > mother_hist(nbEntries, hist);

  vector <double> arr_2(0);
  vector<vector<double> >  hist_2(110,arr);
  vector<vector<vector<double> > > PMT_hist(nbEntries, hist_2);


  for(long unsigned int iEntry = 0; iEntry < nbEntries; iEntry++){
    // Point to event iEntry inside WCSimTree
    wcsimT->GetEvent(iEntry);
    // Nb of Trigger inside the event
    const unsigned int nbTriggers = wcsimroothyperevent->GetNumberOfEvents();
    const unsigned int nbSubTriggers = wcsimroothyperevent->GetNumberOfSubEvents();

    for(long unsigned int iTrig = 0; iTrig < nbTriggers; iTrig++){
      WCSimRootTrigger *wcsimrootevent = wcsimroothyperevent->GetTrigger(iTrig);

      // RAW HITS
      int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
        for (int i = 0; i < ncherenkovdigihits; i++){
          WCSimRootCherenkovDigiHit *hit = (WCSimRootCherenkovDigiHit*)
          (wcsimrootevent->GetCherenkovDigiHits()->At(i));

          WCSimRootCherenkovDigiHit *cDigiHit = wcsimrootevent->GetCherenkovDigiHits()->At(i);

          double PhotoElectrons = cDigiHit->GetQ();

	  int tubeID = hit->GetTubeId() -1;
	  double tube[3];
	  int cylLoc = geo->GetPMT(tubeID).GetCylLoc();
	  tube[0] = geo->GetPMT(tubeID).GetPosition(0);
	  tube[1] = geo->GetPMT(tubeID).GetPosition(1);
	  tube[2] = geo->GetPMT(tubeID).GetPosition(2);

          double timing = hit->GetT();
          // WCSimRootPMT pmt = wcsimrootgeom->GetPMT(tubeId);
          // cout << "me " << endl;
          // cout << "pmt ID: " << tubeId << endl;

          double pmtZ=tube[2];

          double real_z = pmtZ + 2750.;
          double index = real_z / 5500.;
          double new_index = index * 110.;
          double f_index = floor(new_index);

          double K_value = 0.0004;
          double d_w = 54.2523298687;
          double L_att = 52015.;
          double expon = exp(d_w / L_att);
          double f_theta = pow(cos(TMath::Pi() * 47/180), 2);

          double q_coor = (K_value * PhotoElectrons * d_w * expon) / f_theta;
          mother_hist[iEntry][f_index].push_back(PhotoElectrons);
          PMT_hist[iEntry][f_index].push_back(1.0);




          } // END FOR RAW HITS


    } // END FOR iTRIG

  } // END FOR iENTRY

vector <double> indi(110);
vector<vector<double> >  all_Q_coor(nbEntries,indi);

for (int per_event = 0; per_event < mother_hist.size(); ++per_event){

  for (int bins=0; bins < mother_hist[per_event].size(); ++bins){
    double N_pmt = PMT_hist[per_event][bins].size();
    // cout << "N_pmt" << N_pmt << endl;
    double sum_q = 0;

    for (int loop_PMT=0; loop_PMT < mother_hist[per_event][bins].size(); ++loop_PMT){

      // cout << "sum_q entries: " << sum_q << endl;
      sum_q += mother_hist[per_event][bins][loop_PMT];
    }
    if (sum_q == 0){
      double test_Q = 0;
    }
    else {
      test_Q = sum_q / N_pmt;
      // cout << "test_Q: " << test_Q << endl;
    }
    // cout << "sum_q: " << sum_q << endl;
    // cout << "test_Q: " << test_Q << endl;
    all_Q_coor[per_event][bins] = test_Q;
    if (N_pmt == 0){
      all_Q_coor[per_event][bins] = 0;
    }

  }

}

// for (int perevt = 0; perevt < all_Q_coor.size(); ++perevt;){
//   for (int bins = 0; bins < all_Q_coor[perevt].size(); ++bins){
//     cout << "entries: " << all_Q_coor[perevt][bins] << endl;
//   }
// }

// for (int perevt = 0; perevt < PMT_hist.size(); ++perevt){
//   // cout << "first vector size: " << PMT_hist.size() << endl;
//   for (int bins = 0; bins < PMT_hist[perevt].size(); ++bins){
//     // cout << "second vector size: " << PMT_hist[perevt].size() << endl;
//
//     for (int loop = 0; loop < PMT_hist[perevt][bins].size(); ++loop){
//       // cout << "third vector size: " << PMT_hist[perevt][bins].size() << endl;
//       cout << "size: " << PMT_hist[perevt][bins].size() << endl;
//       // cout << "loop: " << loop << "elements: " << PMT_hist[perevt][bins][loop] << endl;
//     }
//   }
// }


vector <double> indi_2(110);
vector<vector<double> >  uncer_all_Q_coor(nbEntries,indi_2);


  for (int per_evt=0; per_evt < mother_hist.size(); ++per_evt){
    for (int bins=0; bins < mother_hist[per_evt].size(); ++bins){
      double N_pmt = pow(PMT_hist[per_evt][bins].size(), 2);

      double sum_q_raw = 0;

      for (int loop_PMT=0; loop_PMT < mother_hist[per_evt][bins].size(); ++loop_PMT){

        sum_q_raw += mother_hist[per_evt][bins][loop_PMT];
      }

        uncer_all_Q_coor[per_evt][bins]= sum_q_raw / N_pmt;
        if (N_pmt == 0){
          uncer_all_Q_coor[per_evt][bins] = 0;

        }
    }

}

// for (int per_evt=0; per_evt < uncer_all_Q_coor.size(); ++per_evt){
//   for (int bins=0; bins < uncer_all_Q_coor[per_evt].size(); ++bins){
//     cout << "entries: " << uncer_all_Q_coor[per_evt][bins] << endl;
//   }
// }


double numer = 0;
double denom = 0;
vector <double> average_hist(uncer_all_Q_coor.size());
  for (int evt=0; evt < uncer_all_Q_coor.size(); ++evt){

    for (int bin_l = 3; bin_l < 106; ++bin_l){

      if (uncer_all_Q_coor[evt][bin_l]==0){
        numer = 0;
        denom = 0;
      }

      else {
        // cout << all_Q_coor[evt][bin_l] << endl;
        numer += all_Q_coor[evt][bin_l] / uncer_all_Q_coor[evt][bin_l];
        // cout << "numer: " << numer << endl;

        denom += 1 / uncer_all_Q_coor[evt][bin_l];
        // cout << "denom: " << denom << endl;
       }

      }

      double div = numer / denom;
      // cout << "div: " << div << endl;
      average_hist[evt] = div;

    }

  for (int evt=0; evt < all_Q_coor.size(); ++evt){
    double rh = 0;
    for (int bin_l = 3; bin_l < 106; ++bin_l){
      int N_bins = all_Q_coor[evt].size();

      if (uncer_all_Q_coor[evt][bin_l]==0){
        rh += 0;
      }
      else {
        rh += pow((all_Q_coor[evt][bin_l] - average_hist[evt]) / uncer_all_Q_coor[evt][bin_l], 2);
      }

    }

    double chi = rh / (N_bins - 6);
    // cout << "chi: " << chi << endl;
    histo_2->Fill(chi, average_hist[evt]);

  }

  TH1 *temp;
    float win_scale=0.75;
    int n_wide=1;
    int n_high=1;
    TCanvas *c1 = new TCanvas("c1","c1",800*n_wide*win_scale,800*n_high*win_scale);
    c1->Divide(n_wide,n_high);

    c1->cd(1);
    histo_2->Draw("p");
    histo_2->SetMarkerStyle(30);
    histo_2->SetMarkerSize(1.0);
    histo_2->SetMarkerColor(6);

  }
