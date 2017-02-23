/**
 * \file MatrixComposition.h
 *
 * \ingroup CrossSectionMatrix
 * 
 * \brief Class def header for a class MatrixComposition
 *
 * @author jzennamo
 */

/** \addtogroup CrossSectionMatrix

    @{*/

#ifndef LARLITE_MATRIXCOMPOSITION_H
#define LARLITE_MATRIXCOMPOSITION_H

#include "Analysis/ana_base.h"
#include "DataFormat/mcpart.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mceventweight.h"

#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"


namespace larlite {
  /**
     \class MatrixComposition
     User custom analysis class made by SHELL_USER_NAME
   */
  class MatrixComposition : public ana_base{
  
  public:

    /// Default constructor
    MatrixComposition(){ _name="MatrixComposition"; _fout=0;}

    /// Default destructor
    virtual ~MatrixComposition(){}

    /** IMPLEMENT in MatrixComposition.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MatrixComposition.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MatrixComposition.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();


    TH1D* weight_nue; 
    TH1D* weight_numu;

    TH1D* Enu_nue; 
    std::vector<TH1D*> syst_Enu_nue;

    TH1D* Enu_numu; 
    std::vector<TH1D*> syst_Enu_numu;

    TH2D* cov_nue;
    TH2D* cov_numu;
    TH2D* cov_numunue;

    TH2D* corr_nue;
    TH2D* corr_numu;
    TH2D* corr_numunue;

    TH2D* weights_v_energy_numu;


  protected:
    
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
