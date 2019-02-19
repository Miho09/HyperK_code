

#include <stdio.h>
#include <stdlib.h>
#include <TH3.h>
#include <TH1.h>
#include <math.h>
#include <vector>
#include <iostream>

#include <iterator>
// #include <tuple>
using std::vector;
// using std::tuple;
// using std::get;



void carte_polar(char *filename=NULL) {
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

  TH1D *histo = new TH1D("Q_coor","Q_coor", 100,0,110);
  histo->SetXTitle("Q_coor per 50cm bin");

  TH1D *histo_2 = new TH1D("Uncertainty", "Q_coor Uncertainty squared",100,0,1);
  histo_2->SetXTitle("Q_coor uncertainty");


  const long unsigned int nbEntries = wcsimT->GetEntries();

  vector <double> arr(0);
  vector<vector<double> >  hist(110,arr);
  vector<vector<vector<double> > > mother_hist(nbEntries, hist);

  vector <double> uncer(0);
  vector<vector<double> >  uncert_hist(110,uncer);
  vector<vector<vector<double> > > mother_uncert_hist(nbEntries, uncert_hist);

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

          WCSimRootTrack * track = (WCSimRootTrack*) wcsimrootevent->GetTracks()->At(0);

          vtxX = wcsimrootevent->GetVtx(0);
          vtxY = wcsimrootevent->GetVtx(1);
          vtxZ = wcsimrootevent->GetVtx(2);
          dirX = track->GetDir(0);
          dirY = track->GetDir(1);
          dirZ = track->GetDir(2);
          energy = track->GetE();
          theta = abs( atan( dirX/dirZ ) );
          phi = abs( atan( dirY/dirX) );


          WCSimRootCherenkovDigiHit *hit = (WCSimRootCherenkovDigiHit*)
          (wcsimrootevent->GetCherenkovDigiHits()->At(i));

          WCSimRootCherenkovDigiHit *cDigiHit = wcsimrootevent->GetCherenkovDigiHits()->At(i);

          double PhotoElectrons = cDigiHit->GetQ();
          // cout << "1" << endl;
          int tubeId = hit->GetTubeId();
          double timing = hit->GetT();
          WCSimRootPMT pmt = wcsimrootgeom->GetPMT(tubeId);

          double pmtZ = pmt.GetPosition(2);
          double real_z = pmtZ + 2750.;
          double index = real_z / 5500.;
          double new_index = index * 110.;
          double f_index = floor(new_index);
          // cout << "f_index: " << f_index << endl;

          double K_value = 0.0004;
          // double d_w = (real_z) / (cos(TMath::Pi() * 43/180));
          double d_w = 54.2523298687;
          double L_att = 52015.;
          double expon = exp(d_w / L_att);
          double f_theta = pow(cos(TMath::Pi() * 47/180), 2);

          double q_coor = (K_value * PhotoElectrons * d_w * expon) / f_theta;

          mother_hist[iEntry][f_index].push_back(PhotoElectrons);
          mother_uncert_hist[iEntry][f_index].push_back(PhotoElectrons);


          } // END FOR RAW HITS


    } // END FOR iTRIG

  } // END FOR iENTRY


// for (int bins = 0; bins < hist.size(); ++bins){
//   for (int entry = 0; entry < hist[bins].size(); ++entry){
//     cout << hist[bins][entry] << endl;
//   }
//   // cout << hist[i]
// }

// for (int bins = 0; bins < mother_uncert_hist.size(); ++bins){
//   for (int entry = 0; entry < mother_uncert_hist[bins].size(); ++entry){
//     cout << "length of mother_hist: " << mother_hist[bins].size() << endl;
//     cout << "length of mother_uncert_hist: " << mother_uncert_hist[bins].size() << endl;
//
//     // cout << "length of mother_hist: " << mother_hist[bins][entry].size() << endl;
//     // cout << "length of mother_uncert_hist: " << mother_uncert_hist[bins][entry].size() << endl;
//     //
//     // for (int find = 0; find < mother_uncert_hist[bins][entry].size(); ++find){
//     //   cout << "length of mother_uncer_hist: " << mother_uncert_hist[bins][entry][find].size() << endl;
//     //   cout << "mother_hist: " << mother_hist[bins][entry][find] << endl;
//     // }
//   }
// }


vector <double> indi(110);
vector<vector<double> >  all_Q_coor(nbEntries,indi);

for (int per_event = 0; per_event < mother_hist.size(); ++per_event){
  for (int bins=0; bins < mother_hist[per_event].size(); ++bins){
    double N_pmt = mother_hist[per_event][bins].size();
    cout << "N_pmt: " << N_pmt << endl;
    double sum_q = 0;

    for (int loop_PMT=0; loop_PMT < mother_hist[per_event][bins].size(); ++loop_PMT){
      sum_q += mother_hist[per_event][bins][loop_PMT];
    }
    if (sum_q == 0){
      double test_Q = 0;
    }
    else {
      test_Q = sum_q / N_pmt;
    }
    // cout << "test_Q: " << test_Q << endl;
    all_Q_coor[per_event][bins]= sum_q / N_pmt;

    if (N_pmt == 0){
      all_Q_coor[per_event][bins] = 0;
    }

  }

}

vector <double> indi_2(110);
vector<vector<double> >  uncer_all_Q_coor(nbEntries,indi_2);


  for (int per_evt=0; per_evt < mother_uncert_hist.size(); ++per_evt){
    for (int bins=0; bins < mother_uncert_hist[per_evt].size(); ++bins){
      double N_pmt = pow(mother_uncert_hist[per_evt][bins].size(),2);

      double sum_q_coor_2 = 0;
      double sum_q_raw = 0;
      double sum_q_tot = 0;

      for (int loop_PMT=0; loop_PMT < mother_uncert_hist[per_evt][bins].size(); ++loop_PMT){

        sum_q_coor_2 += pow(mother_hist[per_evt][bins][loop_PMT],2);

        sum_q_raw += mother_uncert_hist[per_evt][bins][loop_PMT];
        // cout << "1" << endl;

        // cout << "sum_q_raw: " << endl;
        sum_q_tot += sum_q_coor_2 / sum_q_raw;
        // cout << "7" << endl;

        // cout << "sum_q_tot: " << sum_q_tot << endl;
      }
        uncer_all_Q_coor[per_evt][bins]= sum_q_tot / N_pmt;
        if (N_pmt == 0){
          uncer_all_Q_coor[per_evt][bins] = 0;
          // cout << "8" << endl;

        }
    }

  }
  // cout << "9" << endl;
  // for (int per_evt = 0; per_evt < all_Q_coor.size(); ++per_evt){
  //   for (int entry = 0; entry < all_Q_coor[per_evt].size(); ++entry){
  //     cout << "length of all_Q_coor: " << all_Q_coor[per_evt][entry] << endl;
  //
  //     cout << "length of uncer_all_Q_coor: " << uncer_all_Q_coor[per_evt][entry] << endl;
  //       // histo -> Fill(all_Q_coor[per_evt][entry]);
  //       // histo_2 -> Fill(uncer_all_Q_coor[per_evt][entry]);
  //     }
  //   }


  for (int per_evt = 0; per_evt < all_Q_coor.size(); ++per_evt){
    for (int entry = 0; entry < all_Q_coor[per_evt].size(); ++entry){
        histo -> Fill(all_Q_coor[per_evt][entry]);
        histo_2 -> Fill(uncer_all_Q_coor[per_evt][entry]);
      }
    }



  // for (int per_e = 0; per_e < all_Q_coor.size(); ++per_e){
  //   for (int bins = 0; bins < all_Q_coor[per_e].size(); ++bins){
  //      histo -> Fill(all_Q_coor[per_e][bins])
  //   }
  // }
  // for (int per_e = 0; per_e < mother_uncert_hist.size(); ++per_e){
  //   for (int bins = 0; bins < mother_uncert_hist[per_e].size(); ++bins){
  //      // histo_2 -> Fill(uncer_all_Q_coor[per_e][bins])
  //      for (int loop_PMT=0; loop_PMT < mother_uncert_hist[per_e][bins].size(); ++bins){
  //        cout << "contents: " << mother_uncert_hist[per_e][bins] << endl;
  //      }
  //   }
  // }
  TH1 *temp;
    float win_scale=0.75;
    int n_wide=2;
    int n_high=1;
    TCanvas *c1 = new TCanvas("c1","c1",800*n_wide*win_scale,800*n_high*win_scale);
    c1->Divide(n_wide,n_high);

    // c1->cd(1);
    // QvsT->Draw("colz");

    c1->cd(1);
    histo->Draw("colz");

    c1->cd(2);
    histo_2->Draw("colz");




  }
