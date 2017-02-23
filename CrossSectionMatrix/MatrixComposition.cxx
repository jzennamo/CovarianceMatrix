#ifndef LARLITE_MATRIXCOMPOSITION_CXX
#define LARLITE_MATRIXCOMPOSITION_CXX

#include "MatrixComposition.h"

namespace larlite {

  bool MatrixComposition::initialize() {
    // These are defined in the header

    // Weights to be applied 
    weight_nue  = new TH1D("nue_Weights",  "", 100, -3, 3);
    weight_numu = new TH1D("numu_Weights", "", 100, -3, 3);

    Enu_nue =  new TH1D("Enu_nue" , "", 20, 0.05, 1.5);
    Enu_numu = new TH1D("Enu_numu", "", 20, 0.05, 1.5);

    weights_v_energy_numu = new TH2D("weights_v_energy_numu", "", 20, 0.05, 1.5, 200, -100, 100);
    
    // Covariance Matrix
    cov_nue =  new TH2D("nue_Cov",  "", Enu_nue->GetNbinsX(), 1, Enu_nue->GetNbinsX(), Enu_nue->GetNbinsX(), 1, Enu_nue->GetNbinsX());
    cov_numu = new TH2D("numu_Cov", "", Enu_numu->GetNbinsX(), 1, Enu_numu->GetNbinsX(), Enu_numu->GetNbinsX(), 1, Enu_numu->GetNbinsX());
    cov_numunue = new TH2D("numunue_Cov", "", Enu_numu->GetNbinsX(), 1, Enu_numu->GetNbinsX(), Enu_nue->GetNbinsX(), 1, Enu_nue->GetNbinsX());

    // Correlation Matrix
    corr_nue =  new TH2D("nue_Corr",  "", Enu_nue->GetNbinsX(), 1, Enu_nue->GetNbinsX(), Enu_nue->GetNbinsX(), 1, Enu_nue->GetNbinsX());
    corr_numu = new TH2D("numu_Corr", "", Enu_numu->GetNbinsX(), 1, Enu_numu->GetNbinsX(), Enu_numu->GetNbinsX(), 1, Enu_numu->GetNbinsX());
    corr_numunue = new TH2D("numunue_Corr", "", Enu_numu->GetNbinsX(), 1, Enu_numu->GetNbinsX(), Enu_nue->GetNbinsX(), 1, Enu_nue->GetNbinsX());
    

    // Now we want to resize our vector of weighted histograms
    //   This is a hack-y hack to get it to work quickly 
    syst_Enu_nue.resize(5000);
    syst_Enu_numu.resize(5000);
    
    for(int i = 0; i < 5000; i++){
      
      syst_Enu_nue[i] =  new TH1D(Form("Enu_nue_Uni%d",i+1) , "", 20, 0.05, 1.5);
      syst_Enu_numu[i] = new TH1D(Form("Enu_numu_Uni%d",i+1) , "", 20, 0.05, 1.5);
      
      } 


    return true;
  }
  
  bool MatrixComposition::analyze(storage_manager* storage) {
    
    /// This is how we grab the data products we want to work with

    /// Currently this pulls just the MCTruth information 
    auto mcEvent_v = storage->get_data<event_mctruth>("generator");

    std::vector< double > weights;

    // Iterate through all the events in our 
    for(auto event : *mcEvent_v){
      auto mcWeight_v = storage->get_data<event_mceventweight>("eventweight")->at(0);
      auto mcWeight = mcWeight_v.GetWeights();
      
      for(auto const& it : mcWeight){     
	if(it.first == "piplussplines_PrimaryHadronSplines"){
	  weights = it.second; 
	}
      }
      //// Weights been got. 

      
      auto const& nu = event.GetNeutrino().Nu();
      
      // Electrons 
      if(nu.PdgCode() == 12){
	
	Enu_nue->Fill(nu.Momentum(0).E());

	for(int i = 0; i < weights.size(); i++){
	  weight_nue->Fill(weights[i]);
	  syst_Enu_nue[i]->Fill(nu.Momentum(0).E(),weights[i]);	  
	} 				
      }
      
      // Muons 
      if(nu.PdgCode() == 14){

	Enu_numu->Fill(nu.Momentum(0).E());

	for(int i = 0; i < weights.size(); i++){
	  weight_numu->Fill(weights[i]);
	  syst_Enu_numu[i]->Fill(nu.Momentum(0).E(),weights[i]);	  
	  weights_v_energy_numu->Fill(nu.Momentum(0).E(), weights[i]);	  
	} 				
      }
    }

    


    
    return true;
  }

  bool MatrixComposition::finalize() {

    /// COVARIANCE MATRIX Calculation
 	  /*
	    i && j     = energy bins
	    n          = number of weights
	    N^cv_i     = number of events in bin i for the central value  
	    N^syst_i,m = number of events in bin i for the systematic variation in universe "m"
	    E_ij       = the covariance (square of the uncertainty) for bins i,j

	    E_ij = (1/n) Sum( ( N^cv_i - N^syst_i,m) * ( N^cv_j - N^syst_j,m), m)

	  */

    ///NUMU-NUMU
    for(int i = 1; i < Enu_numu->GetNbinsX()+1; i++){
      for(int j = 1; j < Enu_numu->GetNbinsX()+1; j++){
	
	double element = 0; 
	
	  for(int w = 0; w < 5000; w++){

	  double element_i = (Enu_numu->GetBinContent(i)-syst_Enu_numu[w]->GetBinContent(i));
	  double element_j = (Enu_numu->GetBinContent(j)-syst_Enu_numu[w]->GetBinContent(j));
	  
	  element +=  element_i*element_j;

	  }
	  element /= 5000;
	  
	  cov_numu->SetBinContent(i,j,element);
      }        
    }

    for(int i = 1; i < Enu_numu->GetNbinsX()+1; i++){
      for(int j = 1; j < Enu_numu->GetNbinsX()+1; j++){
	
	double Corr = cov_numu->GetBinContent(i,j)/sqrt(cov_numu->GetBinContent(i,i)*cov_numu->GetBinContent(j,j));
	
	corr_numu->SetBinContent(i,j, Corr);
      }
    }

    ///NUE-NUE
    for(int i = 1; i < Enu_nue->GetNbinsX()+1; i++){
      for(int j = 1; j < Enu_nue->GetNbinsX()+1; j++){
	
	double element = 0; 
	
	  for(int w = 0; w < 5000; w++){

	  double element_i = (Enu_nue->GetBinContent(i)-syst_Enu_nue[w]->GetBinContent(i));
	  double element_j = (Enu_nue->GetBinContent(j)-syst_Enu_nue[w]->GetBinContent(j));
	  
	  element +=  element_i*element_j;

	  }
	  element /= 5000;
	  
	  cov_nue->SetBinContent(i,j,element);
      }        
    }

    for(int i = 1; i < Enu_nue->GetNbinsX()+1; i++){
      for(int j = 1; j < Enu_nue->GetNbinsX()+1; j++){
	
	double Corr = cov_nue->GetBinContent(i,j)/sqrt(cov_nue->GetBinContent(i,i)*cov_nue->GetBinContent(j,j));
	
	corr_nue->SetBinContent(i,j, Corr);
      }
    }


    /// NUMU-NUE
    for(int i = 1; i < Enu_numu->GetNbinsX()+1; i++){
      for(int j = 1; j < Enu_nue->GetNbinsX()+1; j++){
	
	double element = 0; 
	
	  for(int w = 0; w < 5000; w++){

	  double element_i = (Enu_numu->GetBinContent(i)-syst_Enu_numu[w]->GetBinContent(i));
	  double element_j = (Enu_nue->GetBinContent(j)-syst_Enu_nue[w]->GetBinContent(j));
	  
	  element +=  element_i*element_j;

	  }
	  element /= 5000;
	  
	  cov_numunue->SetBinContent(i,j,element);
      }        
    }

    for(int i = 1; i < Enu_numu->GetNbinsX()+1; i++){
      for(int j = 1; j < Enu_nue->GetNbinsX()+1; j++){
	
	double Corr = cov_numunue->GetBinContent(i,j)/sqrt(cov_numunue->GetBinContent(i,i)*cov_numunue->GetBinContent(j,j));
	
	corr_numunue->SetBinContent(i,j, Corr);
      }
    }
        
    weight_nue->Write();
    weight_numu->Write();
    Enu_nue->Write();
    Enu_numu->Write();
    corr_numu->Write();
    corr_nue->Write();
    corr_numunue->Write();
    weights_v_energy_numu->Write();


    return true;
  }

}
#endif
