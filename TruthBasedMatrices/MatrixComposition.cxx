#include <cassert>
#include <cmath>
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcpart.h"
#include "DataFormat/mceventweight.h"
#include "MatrixComposition.h"

namespace larlite {

  MatrixComposition::EventSample::EventSample(std::string _name,
                                              size_t nbins,
                                              double elo, double ehi,
                                              size_t nweights)
      : name(_name), enu(nullptr), cov(nullptr) {
    enu = new TH1D(("enu_" + name).c_str(),
                   ";E_{#nu} [GeV];Entries per bin",
                   nbins, elo, ehi);

    Resize(nweights);
  }

  MatrixComposition::EventSample::~EventSample() {
    delete cov;
    delete enu;
  }

  TGraphErrors* MatrixComposition::EventSample::EnuCollapsed() {
    size_t nbins = enu->GetNbinsX();

    // Compute the mean and standard deviation across universes using
    // Welford's method, cf. Art of Computer Programming (D. Knuth)
    double x[nbins];
    double y[nbins];
    double y0[nbins];
    double s[nbins];
    double s0[nbins];

    for (size_t i=0; i<nbins; i++) {
      x[i] = enu->GetBinCenter(i);
      y[i] = y0[i] = enu_syst[0]->GetBinContent(i+1);
      s[i] = s0[i] = 0.0;

      for (size_t j=1; j<enu_syst.size(); j++) {
        double v = enu_syst[j]->GetBinContent(i);
        y[i] = y0[i] + (v - y0[i]) / j;
        s[i] = s0[i] + (v - y0[i]) * (v - y[i]);

        y0[i] = y[i];
        s0[i] = s[i];
      }

      s[i] = sqrt(s[i] / enu_syst.size());
    }

    return new TGraphErrors(nbins, x, y, nullptr, s);
  }

  void MatrixComposition::EventSample::Resize(size_t nweights) {
    enu_syst.clear();
    for (size_t i=0; i<nweights; i++) {
      std::string hname = Form("enu_%s_%zu", name.c_str(), i);
      TH1D* h = (TH1D*) enu->Clone(hname.c_str());
      enu_syst.push_back(h);
    }
  }

  TH2D* MatrixComposition::EventSample::CovarianceMatrix(
      TH1D* nom, std::vector<TH1D*> syst) {
    int nbins = nom->GetNbinsX();

    TH2D* _cov = new TH2D("cov", "", nbins, 0, nbins, nbins, 0, nbins);

    for (int i=1; i<nbins+1; i++) {
      for (int j=1; j<nbins+1; j++) {
        double vij = 0;
        for (size_t k = 0; k<syst.size(); k++) {
          double vi = nom->GetBinContent(i) - syst[k]->GetBinContent(i);
          double vj = nom->GetBinContent(j) - syst[k]->GetBinContent(j);
          vij += vi * vj / syst.size();
        }
        _cov->SetBinContent(i, j, vij);
      }
    }

    return _cov;
  }

  TH2D* MatrixComposition::EventSample::CovarianceMatrix() {
    delete cov;
    cov = CovarianceMatrix(enu, enu_syst);
    cov->SetName(("cov_" + name).c_str());
    return cov;
  }

  TH2D* MatrixComposition::EventSample::CorrelationMatrix(TH2D* _cov) {
    TH2D* cor = (TH2D*) _cov->Clone("cor");

    for (int i=1; i<_cov->GetNbinsX()+1; i++) {
      for (int j=1; j<_cov->GetNbinsY()+1; j++) {
        double vij = _cov->GetBinContent(i, j);
        double si = sqrt(_cov->GetBinContent(i, i));
        double sj = sqrt(_cov->GetBinContent(j, j));
        cor->SetBinContent(i, j, vij / (si * sj));
      }
    }

    return cor;
  }

  TH2D* MatrixComposition::EventSample::CorrelationMatrix() {
    // Compute the covariance matrix first, if we haven't already
    if (!cov) {
      CovarianceMatrix();
    }

    TH2D* cor = CorrelationMatrix(cov);
    cor->SetName(("cor_" + name).c_str());
    return cor;
  }


  bool MatrixComposition::initialize() {
    samples.push_back(new EventSample("numu"));
    samples.push_back(new EventSample("nue"));

    // FIXME: Move to config file.
    use_weights = {
       "genie_qema_Genie"
      //,"genie_ncelAxial_Genie"
      //,"genie_qevec_Genie"
      //,"genie_ccresAxial_Genie"
      //,"genie_ccresVector_Genie"
    };

    return true;
  }
  
  bool MatrixComposition::analyze(storage_manager* storage) {
    /// Currently this pulls just the MCTruth information 
    auto mcEvent_v = storage->get_data<event_mctruth>("generator");

    std::vector<double> weights;

    // Iterate through all the events we have
    for (auto event : *mcEvent_v) {

      /// Grab the MC event weights
      auto mcWeight_v = storage->get_data<event_mceventweight>("eventweight")->at(0);
      auto mcWeight = mcWeight_v.GetWeights();

      // Iterate through all the weighting functions to compute a set of
      // total weights for this event. mcWeight is a mapping from reweighting
      // function name to a vector of weights for each "universe."
      for (auto const& it : mcWeight) {
        if (weights.empty()) {
          weights.resize(it.second.size(), 1.0);
        }
        else {
          assert(weights.size() == it.second.size());
        }

        // Compute universe-wise product of all requsted weights
        if (use_weights.find("*") != use_weights.end() ||
            use_weights.find(it.first) != use_weights.end()) {
          for (size_t i=0; i<weights.size(); i++) {
            weights[i] *= it.second[i];
          }
        }
      }

      // OK, let's start looking at the neutrino interaction
      auto const& nu = event.GetNeutrino().Nu();
      double nuEnergy = nu.Momentum(0).E();

      // Determine which event sample this event corresponds to. Currently
      // this is based on neutrino PDG code, but this can be extended.
      EventSample* sample = nullptr;

      if (nu.PdgCode() == 12) {
        sample = samples[1];
      }
      else if (nu.PdgCode() == 14) {
        sample = samples[0];
      }

      if (!sample) {
        return true;
      }

      // Fill histograms for this event sample
      if (sample->enu_syst.empty()) {
        sample->Resize(weights.size());
      }
      else {
        assert(sample->enu_syst.size() == weights.size());
      }

      // Neutrino Energy
      sample->enu->Fill(nuEnergy);

      // Fill weights
      for (int i=0; i<weights.size(); i++) {
        sample->enu_syst[i]->Fill(nuEnergy, weights[i]);
      }
    }

    return true;
  }

  bool MatrixComposition::finalize() {
    size_t total_bins = 0;

    // Write out sample=-wise distributions
    for (size_t i=0; i<samples.size(); i++) {
      samples[i]->enu->Write();
      total_bins += samples[i]->enu->GetNbinsX();

      TH2D* cov = samples[i]->CovarianceMatrix();
      cov->Write();

      TH2D* cor = samples[i]->CorrelationMatrix();
      cor->Write();

      TGraphErrors* g = samples[i]->EnuCollapsed();
      g->SetName(("err_" + samples[i]->name).c_str());
      g->Write();
    }

    // Global (sample-to-sample) distributions
    // Build glued-together energy spectra for the nominal and each systematics
    // universe, and feed into the correlation matrix calculator.
    total_bins -= samples.size();
    TH1D hg("hg", ";E_{#nu};Entries per bin", total_bins, 0, total_bins);
    std::vector<TH1D*> hgsys;
    for (size_t i=0; i<samples[0]->enu_syst.size(); i++) {
      hgsys.push_back(new TH1D(Form("hg%zu", i), "", total_bins, 0, total_bins));
    }
    size_t ibin = 0;
    for (size_t i=0; i<samples.size(); i++) {
      for(size_t j=1; j<samples[i]->enu->GetNbinsX()+1; j++) {
        hg.SetBinContent(ibin, samples[i]->enu->GetBinContent(j));
        for (size_t k=0; k<hgsys.size(); k++) {
          hgsys[k]->SetBinContent(ibin, samples[i]->enu_syst[k]->GetBinContent(j));
        }
        ibin++;
      }
    }

    hg.Write();

    TH2D* gcov = EventSample::CovarianceMatrix(&hg, hgsys);
    gcov->Write();

    TH2D* gcor = EventSample::CorrelationMatrix(gcov);
    gcor->Write();

    return true;
  }

}  // namespace larlite

