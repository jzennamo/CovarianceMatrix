#ifndef LARLITE_MATRIXCOMPOSITION_H
#define LARLITE_MATRIXCOMPOSITION_H

/**
 * \file MatrixComposition.h
 *
 * \ingroup CrossSectionMatrix
 * 
 * \brief Class def header for a class MatrixComposition
 *
 * \author mastbaum@uchicago.edu, jzennamo
 */

/** \addtogroup CrossSectionMatrix

    @{*/

#include "Analysis/ana_base.h"
#include <set>
#include <string>
#include <vector>

class TGraphErrors;
class TH1D;
class TH2D;

namespace larlite {
  /**
   * \class MatrixComposition
   */
  class MatrixComposition : public ana_base {
  
  public:

    class EventSample {
      public:
        /**
         * Constructor.
         *
         * Note: You can optionally specify the number of systematics
         * universes up front with the nweights parameter. If this isn't
         * known until runtime, call Resize() later.
         *
         * \param _name String name for the event sample
         * \param nbins Number of energy bins
         * \param elo Lower limit of energy spectrum
         * \param ehi Upper limit of energy spectrum
         * \param nweights Number of systematics universes (see note)
         */
        EventSample(std::string _name="sample", size_t nbins=14,
                    double elo=0.2, double ehi=3.0, size_t nweights=0);

        /** Destructor. */
        ~EventSample();

        /**
         * Get the energy spectrum as a graph with error bars representing
         * the systematic uncertainty.
         */
        TGraphErrors* EnuCollapsed();

        /** Set the number of universes. */
        void Resize(size_t nweights);

        /**
         * Covariance Matrix
         *
	 *   i && j     = energy bins
	 *   n          = number of weights
	 *   N^cv_i     = number of events in bin i for the central value
	 *   N^syst_i,m = number of events in bin i for the systematic
	 *                variation in universe "m"
	 *   E_ij       = the covariance (square of the uncertainty) for
	 *                bins i,j
	 *   E_ij = (1/n) Sum((N^cv_i - N^syst_i,m)*(N^cv_j - N^syst_j,m), m)
	 */
        static TH2D* CovarianceMatrix(TH1D* nom, std::vector<TH1D*> syst);

        /** Covariance matrix using internal histograms */
        TH2D* CovarianceMatrix();

        /** Correlation matrix: Corr[ij] = Cov[ij]/Sqrt(Cov[ii]*Cov[jj]) */
        static TH2D* CorrelationMatrix(TH2D* _cov);

        /** Correlation matrix using internal histograms */
        TH2D* CorrelationMatrix();

        std::string name;  //!< String name for this event sample
        TH1D* enu;  //!< "Nominal" energy spectrum
        std::vector<TH1D*> enu_syst;  //!< Spectra for each systematic universe

      protected:
        TH2D* cov;  //!< Cached covariance matrix
    };

    /** Default constructor */
    MatrixComposition() { _name = "MatrixComposition"; _fout = 0; }

    /** Default destructor */
    virtual ~MatrixComposition() {}

    /** Initialization method to be called before the analysis event loop. */
    virtual bool initialize();

    /** Analyze data event-by-event */
    virtual bool analyze(storage_manager* storage);

    /** Finalize method to be called after all events processed. */
    virtual bool finalize();

  protected:
    std::set<std::string> use_weights;  //!< Weight functions to use
    std::vector<EventSample*> samples;  //!< Event samples
  };

}  // namespace larlite

#endif  // LARLITE_MATRIXCOMPOSITION_H

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 

