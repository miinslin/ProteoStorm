#ifndef SPECTRAL_LIBRARY_H
#define SPECTRAL_LIBRARY_H

#include "SpecSet.h"
#include "PeptideSpectrumMatch.h"
//#include "mzxml.h"

#include <vector>
#include <tr1/tuple>
#include <string>

namespace specnets
{
    typedef std::tr1::shared_ptr<PeptideSpectrumMatch> psmPtr;

    typedef tr1::tuple<Spectrum *, float, int, std::string, float, float, float, psmPtr, float> score_results_tuple;

    class SpectralLibrary : public SpecSet
    {
        public:
            enum MatchScoreType {
                MatchScoreType_DotProduct = 0, MatchScoreType_DotProduct_DeltaD = 1, MatchScoreType_DotProduct_DeltaD_DotBias = 2, MatchScoreType_DotProduct_SLGF = 3
            };

            SpectralLibrary(unsigned int sz = 0){
                specs.resize(sz);
            }

            int createlibrary(  float envelope_score_filter,
                                float pvalue_filter,
                                float fdr_filter,
                                MS2ScoringModel &model,
                                vector<string> &ionsToExtract,
                                string allIons,
                                string aminoacidexclusions,
                                vector<int> charge_filter,
                                bool filter_best_rep);

            int createlibrary(  float envelope_score_filter,
                                float pvalue_filter,
                                MS2ScoringModel &model,
                                vector<string> &ionsToExtract,
                                string allIons,
                                string aminoacidexclusions,
                                vector<int> charge_filter,
                                bool filter_best_rep);

            int projection( string target_annotation,
                            MS2ScoringModel model,
                            vector<string> ions_to_extract,
                            string allIons,
                            Spectrum & outputspectrum);
            virtual int search(Spectrum &query_spec, PeptideSpectrumMatchSet &output_psms, float parentmass_tolerance, float ion_tolerance, int top_hits, int analog_search = 0, float score_threshold = 0.0, int library_search_quality = 1000, int minimum_shared_peaks=1000);
            virtual int search_target_decoy(SpectralLibrary &decoy,
                                            Spectrum query_spec,
                                            psmPtr output_psm,
                                            float parentmz_tolerance,
                                            vector<Spectrum *> target_library_ptr,
                                            vector<Spectrum *> decoy_library_ptr,
                                            int scoring_method);




            int search_target_decoy_SLGF(SpectralLibrary &decoy,
                                            Spectrum query_spec,
                                            psmPtr output_psm,
                                            float parentmz_tolerance,
                                            vector<Spectrum *> target_library_ptr,
                                            vector<Spectrum *> decoy_library_ptr,
                                            int scoring_method,
                                            MS2ScoringModel &model,
                                            vector<string> &ionsToExtract,
                                            string allIons);

            int search_target_decoy_SLGF(SpectralLibrary &decoy,
                                             Spectrum query_spec,
                                             psmPtr output_psm,
                                             float parentmz_tolerance,
                                             vector<Spectrum *> target_library_ptr,
                                             vector<Spectrum *> decoy_library_ptr,
                                             SpectralLibrary &target_isocombined,
                                             SpectralLibrary &decoy_isocombined,
                                             vector<Spectrum *> target_isocombined_library_ptr,
                                             vector<Spectrum *> decoy_isocombined_library_ptr,
                                             int scoring_method,
                                             MS2ScoringModel &model,
                                             vector<string> &ionsToExtract,
                                             string allIons);


            int search_target_decoy_SLGFNew(SpectralLibrary &decoy,
                                             Spectrum query_spec,
                                             vector<psmPtr> & output_psms,
                                             int top_psm_number,
                                             float parentmz_tolerance,
                                             vector<Spectrum *> target_library_ptr,
                                             vector<Spectrum *> decoy_library_ptr,
                                             SpectralLibrary &target_isocombined,
                                             SpectralLibrary &decoy_isocombined,
                                             vector<Spectrum *> target_isocombined_library_ptr,
                                             vector<Spectrum *> decoy_isocombined_library_ptr,
                                             int scoring_method,
                                             MS2ScoringModel &model,
                                             vector<string> &ionsToExtract,
                                             string allIons,
                                             int abundance);

            int search_target_decoy_specset_SLGF(SpectralLibrary &decoy,
                                                                         SpecSet searchable_spectra,
                                                                         float parentmz_tolerance,
                                                                         vector<Spectrum *> target_library_ptr,
                                                                         vector<Spectrum *> decoy_library_ptr,
                                                                         int scoring_method,
                                                                         MS2ScoringModel &model,
                                                                         vector<string> &ionsToExtract,
                                                                         string allIons,
                                                                         int do_score_threshold,
                                                                         float score_threshold,
                                                                         PeptideSpectrumMatchSet &output_psm_set);

            int search_target_decoy_specset_SLGF(SpectralLibrary &decoy,
                                                                         SpecSet searchable_spectra,
                                                                         float parentmz_tolerance,
                                                                         vector<Spectrum *> target_library_ptr,
                                                                         vector<Spectrum *> decoy_library_ptr,
                                                                         SpectralLibrary &target_isocombined,
                                                                         SpectralLibrary &decoy_isocombined,
                                                                         vector<Spectrum *> target_isocombined_library_ptr,
                                                                         vector<Spectrum *> decoy_isocombined_library_ptr,
                                                                         int scoring_method,
                                                                         MS2ScoringModel &model,
                                                                         vector<string> &ionsToExtract,
                                                                         string allIons,
                                                                         int do_score_threshold,
                                                                         float score_threshold,
                                                                         PeptideSpectrumMatchSet &output_psm_set,
                                                                         int output_psm_count,
                                                                         int abundance,
                                                                         int start_seach_idx,
                                                                         int end_search_idx);


            int add_update_spectrum_to_Library(Spectrum & addition_spectrum);

            unsigned int LoadSpecSet_additionalmgf(const char * filename);

            vector<vector<float> > find_global_average_spectrum(MS2ScoringModel model,
                                                                vector<string> ions_to_extract,
                                                                string allIons);

            SpectralLibrary create_decoy_spectral_library(MS2ScoringModel model, vector<string> ions_to_extract, string allIons);

        protected:
            void filter_no_psm();
            void filter_fdr(float fdr_filter);
            void filter_multiple_interpretations();
            void filter_envelope_filter(float envelope_score_filter);
            void filter_pvalue(float pvalue_filter);
            void filter_best_representative(MS2ScoringModel &model,
                                            vector<string> &ionsToExtract,
                                            string allIons);
            void filter_aminoacids(string aminoacidexclusions);
            void filter_forcharge(vector<int> charge);
            void filter_forfragmentation(vector<int> accepted_fragmentation);


            void get_consensus( vector<Spectrum *> spectrum_group,
                                MS2ScoringModel model,
                                vector<string> ionsToExtract,
                                string allIons,
                                Spectrum &consensus_spectrum);

            int get_possible_projections(string target_annotation,
                                         MS2ScoringModel model,
                                         vector<string> ions_to_extract,
                                         string allIons, SpecSet &specs,
                                         vector<float> &projections_set_cosine_depression);

            string create_decoy_peptide(string peptide, int charge);

            void generate_prefix_suffix_peptide(string peptide, vector<string> &prefix, vector<string> &suffix);

    };

    class SpectralLibraryGFSearch : public SpectralLibrary
    {
        public:
            void train_distribution(SpectralLibrary &library, MS2ScoringModel model, vector<string> ions_to_extract, string allIons);

            void generate_distributions(MS2ScoringModel model, vector<string> ions_to_extract, string allIons);
            void save_distributions(string output_path);
            void load_distributions(string input_path);






            //virtual int search(Spectrum query_spec, psmPtr output_psm, float parentmz_tolerance);
            //virtual int search_target_decoy(SpectralLibrary &decoy, Spectrum query_spec, psmPtr output_psm, float parentmz_tolerance, MS2ScoringModel &model, vector<string> &ionsToExtract, string allIons);

        private:
            void train_intensity_bands(SpectralLibrary &library, MS2ScoringModel model, vector<string> ions_to_extract, string allIons);

            void train_cosine_distribution(SpectralLibrary &library, MS2ScoringModel model, vector<string> ions_to_extract, string allIons);

            void train_library_dynamicprogramming(MS2ScoringModel model, vector<string> ions_to_extract, string allIons, SpectralLibrary &library);

            void train_library_montecarlo(MS2ScoringModel model, vector<string> ions_to_extract, string allIons, SpectralLibrary &library);

            void train_debug_data(SpectralLibrary &library);

            vector<vector<pair<float, float> > > histograms;
            vector< pair < float, float > > deletion_probability;
    };
}

#endif
