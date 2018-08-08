#include "PeptideSpectrumMatchSet.h"
#include "SpectralLibrary.h"
#include "projectionutils.h"
#include "spectrum_window_filter.h"
#include "spectrum.h"

//System Includes
#include <map>
#include <vector>
#include <algorithm>

namespace specnets
{

    int SpectralLibrary::createlibrary( float envelope_score_filter,
                                        float pvalue_filter,
                                        MS2ScoringModel &model,
                                        vector<string> &ionsToExtract,
                                        string allIons,
                                        string aminoacidexclusions,
                                        vector<int> charge_filter,
                                        bool filter_best_rep){
        return createlibrary(envelope_score_filter,
                             pvalue_filter,
                             -1,
                             model,
                             ionsToExtract,
                             allIons,
                             aminoacidexclusions,
                             charge_filter,
                             filter_best_rep);

    }

    int SpectralLibrary::createlibrary( float envelope_score_filter,
                                        float pvalue_filter,
                                        float fdr_filter,
                                        MS2ScoringModel &model,
                                        vector<string> &ionsToExtract,
                                        string allIons,
                                        string aminoacidexclusions,
                                        vector<int> charge_filter,
                                        bool filter_best_rep){

        DEBUG_MSG("ORIGINAL LIB SIZE: "<< specs.size());

        //Creating annotation ends
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            list<psmPtr>::iterator it;
            for ( it=specs[spectrum_idx].psmList.begin() ; it != specs[spectrum_idx].psmList.end(); it++ ){
                string annot = (*it)->m_annotation;
                string new_annot = create_annotation_ends(annot);
                (*it)->m_annotation = new_annot;
            }
        }

        cout<<"After Creating Annotation Ends: "<<specs.size()<<endl;

        filter_no_psm();
        cout<<"Filtering No PSM: "<<specs.size()<<endl;

        if(aminoacidexclusions.length() > 0){
            filter_aminoacids(aminoacidexclusions);
            filter_no_psm();
            cout<<"Filter Amino Acids: "<<specs.size()<<endl;
        }


        if(envelope_score_filter > -1){
            filter_envelope_filter(envelope_score_filter);
            filter_no_psm();
            cout<<"Envelope Filter: "<<specs.size()<<endl;
        }


        if(pvalue_filter > -1){
            filter_pvalue(pvalue_filter);
            filter_no_psm();
            cout<<"P Value Filter: "<<specs.size()<<endl;
        }



        if(fdr_filter > -1){
            filter_fdr(fdr_filter);
            filter_no_psm();
            cout<<"FDR Filter of "<<fdr_filter<<": "<<specs.size()<<endl;
        }


        if(false){
            filter_multiple_interpretations();
            filter_no_psm();
            cout<<"Filtering Multiple Interpretations: "<<specs.size()<<endl;
        }

        if(charge_filter.size() != 0){
            filter_forcharge(charge_filter);
            cout<<"Charge Filter: "<<specs.size()<<endl;
        }

        if(false){
            vector<int> accepted_fragmentation;
            accepted_fragmentation.push_back(Spectrum::FragType_CID);
            filter_forfragmentation(accepted_fragmentation);
            cout<<"Fragmentation Filter: "<<specs.size()<<endl;
        }

        if(filter_best_rep){
            filter_best_representative(model, ionsToExtract, allIons);
            cout<<"After Best Rep: "<<specs.size()<<endl;
        }

        bool filter_min_peaks = false;
        if(filter_min_peaks){
            //filter_simple_spectra(model, ionsToExtract, allIons);
            //filter_no_psm();
            cout<<"After Spectra Simplicity: "<<specs.size()<<endl;
        }


        return 0;
    }

    int SpectralLibrary::projection(string target_annotation,
                                    MS2ScoringModel model,
                                    vector<string> ions_to_extract,
                                    string allIons,
                                    Spectrum & outputspectrum){

        //Making sure the annotation looks nice
        target_annotation = create_annotation_ends(target_annotation);

        SpecSet projection_spectra;
        vector<float> projections_set_cosine_depression;
        get_possible_projections(target_annotation,
                                 model,
                                 ions_to_extract,
                                 allIons,
                                 projection_spectra,
                                 projections_set_cosine_depression);

        if(projection_spectra.size() == 0) return -1;

        vector<Spectrum *> spectrum_group;
        for(int spectrum_idx = 0; spectrum_idx < projection_spectra.size(); spectrum_idx++)
            spectrum_group.push_back(&projection_spectra[spectrum_idx]);

        get_consensus(spectrum_group, model, ions_to_extract, allIons, outputspectrum);

        outputspectrum.psmList.clear();
        psmPtr psm(new PeptideSpectrumMatch());
        outputspectrum.psmList.push_back(psm);
        psm->m_annotation = target_annotation;
        psm->m_spectrum = &outputspectrum;


        return 0;
    }


    int SpectralLibrary::search_target_decoy(SpectralLibrary &decoy,
                                             Spectrum query_spec,
                                             psmPtr output_psm,
                                             float parentmz_tolerance,
                                             vector<Spectrum *> target_library_ptr,
                                             vector<Spectrum *> decoy_library_ptr,
                                             int scoring_method){
        AAJumps aajumps(1);
        vector<float> masses;
        int charge;


        PeptideSpectrumMatchSet search_results_decoy;
        vector<score_results_tuple> scores_tuple_decoy;
        //psmPtr psm(new PeptideSpectrumMatch);

        int decoy_start_search_idx;
        int decoy_end_search_idx;
        spectrum_ptr_startend(decoy_library_ptr, query_spec.parentMZ, parentmz_tolerance, decoy_start_search_idx, decoy_end_search_idx);

        for(int library_idx = decoy_start_search_idx; library_idx <= decoy_end_search_idx; library_idx++){

            float library_mass = decoy_library_ptr[library_idx]->parentMZ;
            charge = decoy_library_ptr[library_idx]->parentCharge;

            if(abs(query_spec.parentMZ - library_mass) > parentmz_tolerance) continue;
            if(query_spec.parentCharge != 0 && charge != 0 && query_spec.parentCharge != charge) continue;

            float sim = full_spectrum_similarity(*decoy_library_ptr[library_idx], query_spec);

            float dot_bias = full_spectrum_dotbias(*decoy_library_ptr[library_idx], query_spec, sim);

            score_results_tuple similarity_tuple;
            decoy.specs[library_idx].scan = decoy_library_ptr[library_idx]->scan;
            tr1::get<0>(similarity_tuple) = decoy_library_ptr[library_idx];
            tr1::get<1>(similarity_tuple) = sim;
            tr1::get<2>(similarity_tuple) = decoy_library_ptr[library_idx]->psmList.front()->m_dbIndex;
            tr1::get<3>(similarity_tuple) = (string)decoy_library_ptr[library_idx]->psmList.front()->m_annotation;
            tr1::get<4>(similarity_tuple) = dot_bias;
            tr1::get<7>(similarity_tuple) = decoy_library_ptr[library_idx]->psmList.front();
            scores_tuple_decoy.push_back(similarity_tuple);
        }



        //Finding the best
        //sort(scores.begin(), scores.end(), search_results_comparator);
        sort(scores_tuple_decoy.begin(), scores_tuple_decoy.end(), search_results_comparator);
        for(int i = 0; i < scores_tuple_decoy.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum->scan = psm->m_spectrum->scan;
            psm->m_spectrum = tr1::get<0>(scores_tuple_decoy[i]);
            psm->m_score = tr1::get<1>(scores_tuple_decoy[i]);
            psm->m_dbIndex = tr1::get<2>(scores_tuple_decoy[i]);
            psm->m_annotation = tr1::get<3>(scores_tuple_decoy[i]);

            if(tr1::get<7>(scores_tuple_decoy[i])->m_organism.size() > 0){
                psm->m_organism.push_back(tr1::get<7>(scores_tuple_decoy[i])->m_organism[tr1::get<7>(scores_tuple_decoy[i])->m_organism.size()-1]);
            }

            if(tr1::get<7>(scores_tuple_decoy[i])->m_compound_name.size() > 0){
                psm->m_compound_name.push_back(tr1::get<7>(scores_tuple_decoy[i])->m_compound_name[tr1::get<7>(scores_tuple_decoy[i])->m_compound_name.size()-1]);
            }

            search_results_decoy.push_back(psm);
        }


        PeptideSpectrumMatchSet search_results;
        vector<score_results_tuple> scores_tuple;

        int target_start_search_idx;
        int target_end_search_idx;
        spectrum_ptr_startend(target_library_ptr, query_spec.parentMZ, parentmz_tolerance, target_start_search_idx, target_end_search_idx);

        for(int library_idx = target_start_search_idx; library_idx <= target_end_search_idx; library_idx++){

            float library_mass = target_library_ptr[library_idx]->parentMZ;
            charge = target_library_ptr[library_idx]->parentCharge;

            //DEBUG_MSG(query_spec.parentMZ<<"\t"<<library_mass<<"\t"<<charge<<"\t"<<query_spec.parentCharge<<"\t"<<target_library_ptr[library_idx]->psmList.front()->m_compound_name[0]);

            if(abs(query_spec.parentMZ - library_mass) > parentmz_tolerance) continue;
            if(query_spec.parentCharge != 0 && charge != 0 && query_spec.parentCharge != charge) continue;
            float sim = full_spectrum_similarity(*target_library_ptr[library_idx], query_spec);
            float dot_bias = full_spectrum_dotbias(*target_library_ptr[library_idx], query_spec, sim);

            score_results_tuple similarity_tuple;

            tr1::get<0>(similarity_tuple) = target_library_ptr[library_idx];
            tr1::get<1>(similarity_tuple) = sim;
            tr1::get<2>(similarity_tuple) = target_library_ptr[library_idx]->psmList.front()->m_dbIndex;
            tr1::get<3>(similarity_tuple) = (string)target_library_ptr[library_idx]->psmList.front()->m_annotation;
            tr1::get<4>(similarity_tuple) = dot_bias;
            tr1::get<7>(similarity_tuple) = target_library_ptr[library_idx]->psmList.front();

            scores_tuple.push_back(similarity_tuple);
        }


        //Finding the best
        sort(scores_tuple.begin(), scores_tuple.end(), search_results_comparator);
        for(int i = 0; i < scores_tuple.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum = tr1::get<0>(scores_tuple[i]);
            psm->m_score = tr1::get<1>(scores_tuple[i]);
            psm->m_dbIndex = tr1::get<2>(scores_tuple[i]);
            psm->m_annotation = tr1::get<3>(scores_tuple[i]);

            if(tr1::get<7>(scores_tuple[i])->m_organism.size() > 0){
                psm->m_organism.push_back(tr1::get<7>(scores_tuple[i])->m_organism[tr1::get<7>(scores_tuple[i])->m_organism.size()-1]);
            }

            if(tr1::get<7>(scores_tuple[i])->m_compound_name.size() > 0){
                psm->m_compound_name.push_back(tr1::get<7>(scores_tuple[i])->m_compound_name[tr1::get<7>(scores_tuple[i])->m_compound_name.size()-1]);
            }

            search_results.push_back(psm);
        }

        float target_top_scoring = 0.f;
        float target_second_scoring = 0.f;
        float decoy_top_scoring = 0.f;
        float decoy_second_scoring = 0.f;

        if(search_results.size() > 0){
            target_top_scoring = search_results[0]->m_score;
        }

        if(search_results.size() > 1){
            target_second_scoring = search_results[1]->m_score;
        }

        if(search_results_decoy.size() > 0){
            decoy_top_scoring = search_results_decoy[0]->m_score;
        }

        if(search_results_decoy.size() > 1){
            decoy_second_scoring = search_results_decoy[1]->m_score;
        }


        if(target_top_scoring > decoy_top_scoring && search_results.size() > 0){
            float dot_bias = tr1::get<4>(scores_tuple[0]);
            float dot_product = target_top_scoring;
            float deltaD = target_top_scoring - max(target_second_scoring, decoy_top_scoring);
            float match_score = target_top_scoring * 0.6 + deltaD * 0.4 - dot_bias;

            //cout<<"Dot Bias\t"<<dot_bias<<"\t"<<search_results[0]->m_annotation<<"\tdot\t"<<dot_product<<"\tdeltaD\t"<<deltaD<<endl;

            switch(scoring_method){
                case SpectralLibrary::MatchScoreType_DotProduct:
                    match_score = dot_product;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD:
                    match_score = dot_product * 0.6 + deltaD * 0.4;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD_DotBias:
                    match_score = dot_product * 0.6 + deltaD * 0.4 - dot_bias;
                    break;
                default:
                    break;
            }

            output_psm->m_spectrum      = search_results[0]->m_spectrum;
            output_psm->m_score         = search_results[0]->m_score;
            output_psm->m_dbIndex       = search_results[0]->m_dbIndex;
            output_psm->m_annotation    = search_results[0]->m_annotation;
            output_psm->m_spectrumFile  = search_results[0]->m_spectrum->fileName;
            output_psm->m_organism      = search_results[0]->m_organism;
            output_psm->m_compound_name = search_results[0]->m_compound_name;
            //output_psm->m_organism.insert(output_psm->m_organism.end(), search_results[0]->m_organism.begin(), search_results[0]->m_organism.end());
            //output_psm->m_compound_name.insert(output_psm->m_compound_name.end(), search_results[0]->m_compound_name.begin(), search_results[0]->m_compound_name.end());

            output_psm->m_score         = match_score;
            output_psm->m_isDecoy = false;

            return 0;
        }
        if(target_top_scoring < decoy_top_scoring && search_results_decoy.size() > 0){
            float dot_bias = tr1::get<4>(scores_tuple_decoy[0]);
            float dot_product = decoy_top_scoring;
            float deltaD = decoy_top_scoring - max(target_top_scoring, decoy_second_scoring);
            float match_score = decoy_top_scoring * 0.6 + deltaD * 0.4 - dot_bias;

            switch(scoring_method){
                case SpectralLibrary::MatchScoreType_DotProduct:
                    match_score = dot_product;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD:
                    match_score = dot_product * 0.6 + deltaD * 0.4;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD_DotBias:
                    match_score = dot_product * 0.6 + deltaD * 0.4 - dot_bias;
                    break;
                default:
                    break;
            }

            output_psm->m_spectrum      = search_results_decoy[0]->m_spectrum;
            output_psm->m_score         = search_results_decoy[0]->m_score;
            output_psm->m_dbIndex       = search_results_decoy[0]->m_dbIndex;
            output_psm->m_annotation    = search_results_decoy[0]->m_annotation;
            output_psm->m_spectrumFile  = search_results_decoy[0]->m_spectrum->fileName;
            output_psm->m_organism      = search_results_decoy[0]->m_organism;
            output_psm->m_compound_name = search_results_decoy[0]->m_compound_name;
            output_psm->m_score         = match_score;
            output_psm->m_isDecoy = true;

            return 0;
        }

        return -1;
    }


    /*! \brief Default Spectral Library Search
        Simply returns top_hits using full spectrum cosine scoring
     */

    int SpectralLibrary::search(Spectrum &query_spec, PeptideSpectrumMatchSet &output_psms, float parentmass_tolerance, float ion_tolerance, int top_hits, int analog_search, float score_threshold, int library_search_quality, int minimum_shared_peaks){
        //AAJumps aajumps(1);
        vector<float> masses;
        int charge;

        PeptideSpectrumMatchSet search_results;
        vector<score_results_tuple> scores_tuple;

        //psmPtr psm(new PeptideSpectrumMatch);

        for(int library_idx = 0; library_idx < specs.size(); library_idx++){
            //masses.clear();
            //aajumps.getPRMMasses(specs[library_idx].psmList.front()->m_annotation.c_str(), masses);


            charge = specs[library_idx].parentCharge;
            //float library_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
            float library_mass = specs[library_idx].parentMZ;
            float query_mass = query_spec.parentMZ;

            float mz_error = abs(query_mass - library_mass);
            float mass_error = abs(specs[library_idx].parentMass - query_spec.parentMass);


            if(mz_error > parentmass_tolerance and analog_search == 0) continue;

            if(mass_error > analog_search and analog_search > 0) continue;

            if(specs[library_idx].psmList.size() > 0){
                if(specs[library_idx].psmList.front()->m_library_quality > library_search_quality){
                    continue;
                }
            }


            float ms_error_ppm = mz_error/library_mass * 1000000;

            unsigned int shared_peaks = 0;


            //float sim = full_spectrum_similarity(specs[library_idx], query_spec, shared_peaks);
            float score1, score2;
            float sim = specs[library_idx].scoreMatch(query_spec, ion_tolerance, shared_peaks, score1, score2);

            if(sim != sim) sim = 0.f;

            //cout<<"DEBUG_BACKGROUND:\t"<<query_mass<<"\t"<<library_mass<<"\t"<<mz_error<<"\t"<<ms_error_ppm<<"\t"<<query_spec.scan<<"\t"<<specs[library_idx].scan<<"\t"<<abs((int)((int)specs[library_idx].scan - (int)query_spec.scan))<<"\t"<<sim<<"\t"<<shared_peaks<<endl;


            if(sim < score_threshold) continue;

            if(shared_peaks < minimum_shared_peaks) continue;

            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum = &query_spec;
            psm->m_score = sim;
            psm->m_dbIndex = specs[library_idx].psmList.front()->m_dbIndex;
            psm->m_annotation = (string)specs[library_idx].psmList.front()->m_annotation;
            psm->m_spectrumFile = query_spec.fileName;
            psm->m_organism = specs[library_idx].psmList.front()->m_organism;
            psm->m_compound_name = specs[library_idx].psmList.front()->m_compound_name;
            psm->m_scanNum = query_spec.scan;
            psm->m_library_name = get_only_filename(specs[library_idx].psmList.front()->m_spectrumFile);
            psm->m_spectrumID = specs[library_idx].psmList.front()->m_spectrumID;
            psm->m_mz_error_ppm = ms_error_ppm;
            psm->m_shared_peaks = shared_peaks;
            psm->m_mz = query_spec.parentMZ;


            //Metadata temp for Metabolomic searches
            psm->m_fdr = specs[library_idx].ITOL;                                       //RT of the library
            psm->m_pValue = query_spec.retention_time;                                  //RT of the query
            psm->m_strict_envelope_score = specs[library_idx].getTotalIonCurrent();     //TIC of library
            psm->m_unstrict_envelope_score = query_spec.getTotalIonCurrent();           //TIC of query
            psm->m_startMass = specs[library_idx].parentMZ;                             //MZ of match
            psm->m_parentmass_difference = mass_error;                                  //Mass Difference
            stringstream ss;
            ss << specs[library_idx].scan;
            psm->m_protein = specs[library_idx].psmList.front()->m_notes + ":" + ss.str();               //filename



            if(specs[library_idx].psmList.front()->m_charge == 0){
                psm->m_charge = query_spec.parentCharge;
            }
            else{
                psm->m_charge = specs[library_idx].psmList.front()->m_charge;
            }

            search_results.push_back(psm);


            /*
            score_results_tuple similarity_tuple;

            tr1::get<0>(similarity_tuple) = &specs[library_idx];
            tr1::get<1>(similarity_tuple) = sim;
            tr1::get<2>(similarity_tuple) = library_idx;
            tr1::get<3>(similarity_tuple) = (string)specs[library_idx].psmList.front()->m_annotation;
            scores_tuple.push_back(similarity_tuple);
            */
        }

        //Finding the best
        //sort(scores.begin(), scores.end(), search_results_comparator);
        sort(search_results.m_psmSet.begin(), search_results.m_psmSet.end(), search_results_comparator_psmPtr);



        for(int i = 0; i < min((int)search_results.size(), top_hits); i++){
            output_psms.push_back(search_results[i]);
        }

        return 0;

        /*sort(scores_tuple.begin(), scores_tuple.end(), search_results_comparator);

        for(int i = 0; i < scores_tuple.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum = tr1::get<0>(scores_tuple[i]);
            psm->m_score = tr1::get<1>(scores_tuple[i]);
            psm->m_scanNum = tr1::get<2>(scores_tuple[i]) + 1;
            psm->m_annotation = tr1::get<3>(scores_tuple[i]);
            search_results.push_back(psm);
        }*/

        if(search_results.size() > 0){
            //float match_score = 0.6*search_results[0]->m_score;
            //if(search_results.size() == 1) match_score += 0.4 * search_results[0]->m_score;
            //else match_score += 0.4 * search_results[1]->m_score;

            //output_psm->m_spectrum      = search_results[0]->m_spectrum;
            //output_psm->m_score         = search_results[0]->m_score;
            //output_psm->m_scanNum       = search_results[0]->m_scanNum;
            //output_psm->m_annotation    = search_results[0]->m_annotation;
            //output_psm->m_score         = match_score;

            return 0;
        }

        /*for(int score_idx = 0; score_idx < min((int)scores_tuple.size(), 1); score_idx++){
            cout<<tr1::get<0>(scores_tuple[score_idx])->psmList.front()->m_annotation<<"\t"<<tr1::get<1>(scores_tuple[score_idx])<<"\t0BasedIdx\t"<<tr1::get<2>(scores_tuple[score_idx])<<endl;
        }*/

        return -1;
    }

    /*! \brief Loads an mgf file and adds it to current specset

     */
    unsigned int SpectralLibrary::LoadSpecSet_additionalmgf(const char * filename){
        SpecSet tempspecs;
        tempspecs.LoadSpecSet_mgf(filename);

        DEBUG_MSG("Loaded Spec MGF");

        for(int specidx = 0; specidx < tempspecs.size(); specidx++){
            specs.push_back(tempspecs[specidx]);
            list<psmPtr>::iterator it;
            for ( it=specs[specs.size()-1].psmList.begin() ; it != specs[specs.size()-1].psmList.end(); it++ ){
                (*it)->m_spectrum = &(specs[specs.size()-1]);
            }
        }

        return 0;
    }

    void SpectralLibrary::filter_no_psm(){
        SpectralLibrary temp_lib;
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            if(specs[spectrum_idx].psmList.size() >= 1){
                temp_lib.specs.push_back(specs[spectrum_idx]);
            }
        }
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();

        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }
    }

    /*! \brief Filters out PSMs that are not the best scoring PValue

     */
    void SpectralLibrary::filter_multiple_interpretations(){
        SpectralLibrary temp_lib;
        //Filtering out spectra with multiple interpretations, or 0 interpretations
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            if(specs[spectrum_idx].psmList.size() > 1){
                specs[spectrum_idx].psmList.clear();
            }
        }
    }

    /*! \brief Filters out spectra that exceed the fdr threshold

     */
    void SpectralLibrary::filter_fdr(float fdr_filter){
        SpectralLibrary temp_lib;
        //Filtering based on pvalues

        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            list<psmPtr>::iterator it;
            //DEBUG_MSG("MYFDR\t"<<specs[spectrum_idx].psmList.front()->m_fdr);
            while(1){
                bool rerun = false;
                for ( it=specs[spectrum_idx].psmList.begin() ; it != specs[spectrum_idx].psmList.end(); it++ ){
                    if( (*it)->m_fdr > fdr_filter){
                        specs[spectrum_idx].psmList.erase(it);
                        rerun = true;
                        break;
                    }
                }
                if(!rerun){
                    break;
                }
            }
        }
    }

    /*! \brief Filters out spectra that exceed the envelope score

     */
    void SpectralLibrary::filter_envelope_filter(float envelope_score_filter){
        SpectralLibrary temp_lib;
        //Filtering based on envelope scores
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            //cout<<"Envelope: "<<specs[spectrum_idx].psmList.front()->m_strict_envelope_score<<endl;

            list<psmPtr>::iterator it;
            while(1){
                bool rerun = false;
                for ( it=specs[spectrum_idx].psmList.begin() ; it != specs[spectrum_idx].psmList.end(); it++ ){
                    if( (*it)->m_strict_envelope_score > envelope_score_filter){
                        specs[spectrum_idx].psmList.erase(it);
                        rerun = true;
                        break;
                    }
                }
                if(!rerun){
                    break;
                }
            }

            /*

            if(specs[spectrum_idx].psmList.front()->m_strict_envelope_score < envelope_score_filter){
                temp_lib.specs.push_back(specs[spectrum_idx]);
            }*/
        }
        /*
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();

        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }*/
    }

    /*! \brief Filters out spectra that exceed the pvalue threshold

     */
    void SpectralLibrary::filter_pvalue(float pvalue_filter){
        SpectralLibrary temp_lib;
        //Filtering based on pvalues




        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            list<psmPtr>::iterator it;
            while(1){
                bool rerun = false;
                for ( it=specs[spectrum_idx].psmList.begin() ; it != specs[spectrum_idx].psmList.end(); it++ ){
                    if( (*it)->m_pValue > pvalue_filter){
                        specs[spectrum_idx].psmList.erase(it);
                        rerun = true;
                        break;
                    }
                }
                if(!rerun){
                    break;
                }
            }
            /*if(specs[spectrum_idx].psmList.front()->m_pValue < pvalue_filter){
                temp_lib.specs.push_back(specs[spectrum_idx]);
            }*/
        }
        /*
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();

        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }*/
    }


    void SpectralLibrary::filter_forcharge(vector<int> charge){
        SpectralLibrary temp_lib;
        //Filtering based on pvalues
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){

            if(specs[spectrum_idx].parentCharge == 0){
                specs[spectrum_idx].parentCharge = specs[spectrum_idx].psmList.front()->m_charge;
            }

            int spec_charge = specs[spectrum_idx].psmList.front()->m_charge;
            if(spec_charge == -1) spec_charge = specs[spectrum_idx].parentCharge;
            bool valid_charge = false;
            for(int charge_idx = 0; charge_idx < charge.size(); charge_idx++){
                if(spec_charge == charge[charge_idx]){
                    valid_charge = true;
                    break;
                }
            }
            if(valid_charge){
                temp_lib.specs.push_back(specs[spectrum_idx]);
            }
        }
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();

        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }
    }

    void SpectralLibrary::filter_forfragmentation(vector<int> accepted_fragmentation){
        SpectralLibrary temp_lib;
        //Filtering based on fragmentation
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            bool valid_fragmentation = false;
            int spec_fragmentation = specs[spectrum_idx].msFragType;
            for(int fragmentation_idx = 0; fragmentation_idx < accepted_fragmentation.size(); fragmentation_idx++){
                if(spec_fragmentation == accepted_fragmentation[fragmentation_idx]){
                    valid_fragmentation = true;
                    break;
                }
            }
            if(valid_fragmentation){
                temp_lib.specs.push_back(specs[spectrum_idx]);
            }
        }
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();

        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }
    }


    void SpectralLibrary::filter_aminoacids(string aminoacidexclusions){
        SpectralLibrary temp_lib;
        //Filtering based on pvalues
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            //string annotation = specs[spectrum_idx].psmList.front()->m_annotation;

            list<psmPtr>::iterator it;
            while(1){
                bool rerun = false;
                for ( it=specs[spectrum_idx].psmList.begin() ; it != specs[spectrum_idx].psmList.end(); it++ ){
                    if( (*it)->m_annotation.find_first_of(aminoacidexclusions) != string::npos){
                        specs[spectrum_idx].psmList.erase(it);
                        rerun = true;
                        break;
                    }
                }
                if(!rerun){
                    break;
                }
            }
            //if(annotation.find_first_of(aminoacidexclusions) == string::npos){
            //    temp_lib.specs.push_back(specs[spectrum_idx]);
            //}
        }
        /*
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();

        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }*/
    }


    /*! \brief Finds spectrum that best represents each peptide

     */
    void SpectralLibrary::filter_best_representative(   MS2ScoringModel &model,
                                                        vector<string> &ionsToExtract,
                                                        string allIons){
        //Grouping up peptides with the same annotation
        map<string, vector<Spectrum *> > same_annotation_clusters;
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            string annotation = specs[spectrum_idx].psmList.front()->m_annotation;
            same_annotation_clusters[annotation].push_back(&specs[spectrum_idx]);
        }

        cout<<"Unique Annotations:"<<same_annotation_clusters.size()<<endl;

        //Choosing representative peptide
        SpectralLibrary temp_lib;
        map<string, vector<Spectrum *> >::iterator it;
        for ( it=same_annotation_clusters.begin() ; it != same_annotation_clusters.end(); it++ ){
            Spectrum consensus_spectrum;
            get_consensus((*it).second, model, ionsToExtract, allIons, consensus_spectrum);
            temp_lib.specs.push_back(consensus_spectrum);
        }

        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();

        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }
    }

    void SpectralLibrary::get_consensus( vector<Spectrum *> spectrum_group,
                                MS2ScoringModel model,
                                vector<string> ionsToExtract,
                                string allIons,
                                Spectrum &consensus_spectrum){
        //Setting psms to be correct
        for(int spec_idx = 0; spec_idx < spectrum_group.size(); spec_idx++){
            spectrum_group[spec_idx]->psmList.front()->m_spectrum = spectrum_group[spec_idx];
        }


        float best_score = 0.f;
        int best_idx = 0;
        for(int cluster_idx1 = 0; cluster_idx1 < spectrum_group.size(); cluster_idx1++){
            float spec_score = 0.f;
            for(int cluster_idx2 = 0; cluster_idx2 < spectrum_group.size(); cluster_idx2++){
                if(cluster_idx1 == cluster_idx2) continue;
                spec_score += spectrum_similarity(  spectrum_group[cluster_idx1]->psmList.front(),
                                                    spectrum_group[cluster_idx2]->psmList.front(),
                                                    spectrum_group[cluster_idx2]->psmList.front()->m_annotation.length() - 4,
                                                    model,
                                                    ionsToExtract,
                                                    allIons);
            }
            if (spec_score > best_score){
                best_idx = cluster_idx1;
                best_score = spec_score;
            }
        }
        consensus_spectrum = *(spectrum_group[best_idx]);

    }

    int SpectralLibrary::get_possible_projections(string target_annotation, MS2ScoringModel model, vector<string> ions_to_extract, string allIons, SpecSet &projection_specs, vector<float> &projections_set_cosine_depression){
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            string library_annotation = specs[spectrum_idx].psmList.front()->m_annotation;
            if(getStringDifference(target_annotation, library_annotation) == 1){
                Spectrum new_synth_spec;

                psmPtr psm = specs[spectrum_idx].psmList.front();
                string source_annotation = psm->m_annotation;

                vector<pair<float, float> > ion_mass_intensity_pair_vector;

                string annotation = cleanAnnotationEnds(source_annotation);
                psm->annotate(annotation,allIons,model,0,0,.45);

                extractIons(psm,annotation.length()-4,model,ions_to_extract,ion_mass_intensity_pair_vector, 0, 0);

                int diff_location = getDifferenceIndex(annotation, target_annotation) - 2;

                //cout<<"Dif Location: "<<diff_location<<endl;
                //cout<<"Source Annotation: "<<annotation<<endl;

                //Calculating the expected cosine depression
                map<string, float> aa_substitution_map = getAminoAcidSubstitutionLookup();
                //cout<<annotation[diff_location+2]<<"\t"<<target_annotation[diff_location+2]<<endl;
                float cosine_depression = getSubstitutionCosineDepression(annotation[diff_location+2], target_annotation[diff_location+2], aa_substitution_map);
                //cout<<"Cosine Depression: "<<cosine_depression<<endl;

                //If we have a cosine depression of less than a value, we skip this projection
                if(cosine_depression < 0.5){
                    //cout<<"Transformation Not Supported"<<endl;
                    //continue;
                }

                AAJumps aajumps(1);
                vector<float> masses;
                //cout<<"AAJUMPS annot: "<<annotation<<"\t"<<target_annotation<<endl;
                aajumps.getPRMMasses(annotation.c_str(), masses);

                int charge = 2;
                map<char, float> aamasses = getAminoAcidLookup();
                float library_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                //float library_mass = (getMass(annotation, aamasses)+ AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                aajumps.getPRMMasses(target_annotation.c_str(), masses);
                float target_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                //float target_mass = (getMass(target_annotation, aamasses)+ AAJumps::massH2O + AAJumps::massHion*charge)/charge;

                float mass_difference = library_mass - target_mass;

                //cout<<"Library Mass: "<<library_mass<<"\t"<<"Target mass: "<<target_mass<<"\t"<<mass_difference<<endl;
                for(int i = 0; i < ion_mass_intensity_pair_vector.size(); i++){
                    //cout<<ions_to_extract[i/(annotation.length()-4)]<<"\t"<<i%(annotation.length()-4) + 1<<"\t";
                    //cout<<ion_mass_intensity_pair_vector[i].first<<"\t"<<ion_mass_intensity_pair_vector[i].second<<endl;
                }


                for(int i = 0; i < ion_mass_intensity_pair_vector.size(); i++){
                    if(ions_to_extract[i/(annotation.length()-4)].find('b') != -1){
                        if( (i%(annotation.length()-4) + 1) < diff_location + 1) //For B ions don't need to edit
                            continue;

                        if(ions_to_extract[i/(annotation.length()-4)].find("++") != -1){
                            //cout<<ions_to_extract[i/(annotation.length()-4)]<<"\t"<<i%(annotation.length()-4) + 1<<"\t";
                            //cout<<ion_mass_intensity_pair_vector[i].first<<"\t"<<ion_mass_intensity_pair_vector[i].second<<endl;
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                        }
                        else{
                            //cout<<ions_to_extract[i/(annotation.length()-4)]<<"\t"<<i%(annotation.length()-4) + 1<<"\t";
                            //cout<<ion_mass_intensity_pair_vector[i].first<<"\t"<<ion_mass_intensity_pair_vector[i].second<<endl;
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference*2;
                        }
                    }

                    if(ions_to_extract[i/(annotation.length()-4)].find('y') != -1){
                        if( (i%(annotation.length()-4) + 1) < ((annotation.length()-4) - diff_location) ){ //For Y ions don't need to edit
                            //cout<<"Contineu: "<<(i%(annotation.length()-4) + 1)<<endl;
                            continue;
                        }

                        if(ions_to_extract[i/(annotation.length()-4)].find("++") != -1){
                            //cout<<ions_to_extract[i/(annotation.length()-4)]<<"\t"<<i%(annotation.length()-4) + 1<<"\t";
                            //cout<<ion_mass_intensity_pair_vector[i].first<<"\t"<<ion_mass_intensity_pair_vector[i].second<<endl;
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                        }
                        else{
                            //cout<<ions_to_extract[i/(annotation.length()-4)]<<"\t"<<i%(annotation.length()-4) + 1<<"\t";
                            //cout<<ion_mass_intensity_pair_vector[i].first<<"\t"<<ion_mass_intensity_pair_vector[i].second<<endl;
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference*2;
                        }
                    }

                    if(ions_to_extract[i/(annotation.length()-4)].find('P') != -1){
                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                    }
                }

                new_synth_spec.resize(0);
                sort( ion_mass_intensity_pair_vector.begin(), ion_mass_intensity_pair_vector.end(), mass_intensity_pair_mass_comp);
                int new_peaklist_size = 0;
                for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                    if(ion_mass_intensity_pair_vector[p].first > 0.f){
                        new_peaklist_size++;
                        //cout<<ion_mass_intensity_pair_vector[p].first<<" "<<ion_mass_intensity_pair_vector[p].second<<endl;
                    }
                }

                new_synth_spec.resize(new_peaklist_size);
                int new_peaklist_idx = 0;
                for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                    if(ion_mass_intensity_pair_vector[p].second > 0.f){
                        new_synth_spec[new_peaklist_idx][0] = ion_mass_intensity_pair_vector[p].first;
                        new_synth_spec[new_peaklist_idx][1] = ion_mass_intensity_pair_vector[p].second;
                        new_peaklist_idx++;
                    }
                }

                //adding projection to this
                projections_set_cosine_depression.push_back(cosine_depression);
                psmPtr target_psm(new PeptideSpectrumMatch);
                target_psm->m_annotation = target_annotation;
                new_synth_spec.psmList.push_back(target_psm);
                new_synth_spec.parentMZ = target_mass;
                new_synth_spec.parentCharge = specs[spectrum_idx].parentCharge;
                projection_specs.push_back(new_synth_spec);
                target_psm->m_spectrum = &projection_specs[projection_specs.size()-1];
            }
        }
        return 0;
    }


    vector<vector<float> > SpectralLibrary::find_global_average_spectrum(MS2ScoringModel model, vector<string> ions_to_extract, string allIons){
        psmPtr psm(new PeptideSpectrumMatch);
        vector<vector<float> > averages;

        for(int sizes = 0; sizes < 100; sizes++){
            vector<float> average;
            int psm_count = 0;
            for(int i = 0; i < specs.size(); i++){
                if(specs[i].psmList.front()->m_annotation.length()-4 == sizes){
                    psm_count++;
                    //cout<<"I: "<<i<<endl;
                    //cout<<m_psmSet[i]->m_annotation<<endl;
                    Spectrum tempspec = specs[i];
                    //Making Max to 1000
                    preprocess_spectrum_intensities_max_intensity(&tempspec, 1000.f);
                    //Applying SQRT
                    preprocess_spectrum_intensities(&tempspec, 0, 1);

                    psm->m_annotation = specs[i].psmList.front()->m_annotation;
                    psm->m_spectrum = &tempspec;

                    psm->annotate(psm->m_annotation,allIons,model,0,0,.45);

                    vector<float> ion_mass_intensity_pair_vector;

                    extractIons(psm,psm->m_annotation.length()-4,model,ions_to_extract,ion_mass_intensity_pair_vector, 0, 0);
                    norm_vector(ion_mass_intensity_pair_vector);

                    if(average.size() == 0){
                        for(int j = 0; j < ion_mass_intensity_pair_vector.size(); j++){
                            average.push_back(ion_mass_intensity_pair_vector[j]);
                        }
                        continue;
                    }
                    else{
                        for(int j = 0; j < ion_mass_intensity_pair_vector.size(); j++){
                            average[j] += (ion_mass_intensity_pair_vector[j]);
                        }
                    }
                }
            }
            cout<<sizes<<"\t";
            for(int j = 0; j < average.size(); j++){
                average[j] = average[j]/psm_count;
                cout<<average[j]<<"\t";
            }
            cout<<endl;
            averages.push_back(average);
        }

        /*
        average_spectrum->peakList.clear();
        sort( average.begin(), average.end(), mass_intensity_pair_mass_comp);
        int new_peaklist_size = 0;
        for(int p = 0; p < average.size(); p++){
            if(average[p].first > 0.f){
                new_peaklist_size++;
                //cout<<ion_mass_intensity_pair_vector[p].first<<" "<<ion_mass_intensity_pair_vector[p].second<<endl;
            }
        }

        average_spectrum->peakList.resize(new_peaklist_size);
        int new_peaklist_idx = 0;
        for(int p = 0; p < average.size(); p++){
            if(average[p].second > 0.f){
                average_spectrum->peakList[new_peaklist_idx][0] = average[p].first;
                average_spectrum->peakList[new_peaklist_idx][1] = average[p].second;
                new_peaklist_idx++;
            }
        }*/

        return averages;
    }

    SpectralLibrary SpectralLibrary::create_decoy_spectral_library(MS2ScoringModel model, vector<string> ions_to_extract, string allIons){
        unsigned int seed = 0;
        srand(seed);

        vector<string> library_peptides;
        for(int library_idx = 0; library_idx < specs.size(); library_idx++){
            string annotation = specs[library_idx].psmList.front()->m_annotation;
            annotation = create_annotation_ends(annotation);
            string stripped_annotation = annotation.substr(2, annotation.length() - 4);
            stripped_annotation += ('0' + specs[library_idx].parentCharge);
            library_peptides.push_back(stripped_annotation);
        }

        sort(library_peptides.begin(), library_peptides.end());

        SpectralLibrary decoy_library;
        for(int library_idx = 0; library_idx < specs.size(); library_idx++){
            string annotation = specs[library_idx].psmList.front()->m_annotation;

            //Randomize the annotation
            string stripped_annotation = annotation.substr(2, annotation.length() - 4);
            string annotation_orig_stripped = stripped_annotation;

            string random_annotation;

            bool valid_decoy_found = false;

            for(int random_retries = 0 ; random_retries < 50; random_retries++){
                random_annotation = create_decoy_peptide(annotation, specs[library_idx].parentCharge);

                string search_random_annotation = random_annotation;
                search_random_annotation += ('0'  + specs[library_idx].parentCharge);

                if(binary_search(library_peptides.begin(), library_peptides.end(), search_random_annotation)){
                    continue;   //try again, already in library
                }
                else{
                    valid_decoy_found = true;
                    break;
                }
            }

            if(!valid_decoy_found){
                cout<<"No Valid Decoy Found"<<endl;
                continue;
            }

            cout<<decoy_library.size()<<"\t"<<annotation_orig_stripped<<"\t"<<random_annotation<<endl;

            int peptide_length = getpeptideLength(annotation_orig_stripped);

            //Extracting the ions
            vector<pair <float, float> > ion_mass_intensity_pair_vector;
            specs[library_idx].psmList.front()->annotate(annotation,allIons,model,0,0,0.45);
            extractIons(specs[library_idx].psmList.front(), peptide_length, model, ions_to_extract, ion_mass_intensity_pair_vector, 0, 0);


            //Collecting the peaks that are not annotated
            vector<pair <float, float> > unannotated_peaks;
            for(int peak_idx = 0; peak_idx < specs[library_idx].size(); peak_idx++){
                bool annotated = false;
                for(int annotated_idx = 0; annotated_idx < ion_mass_intensity_pair_vector.size(); annotated_idx++){
                    if(specs[library_idx][peak_idx][0] == ion_mass_intensity_pair_vector[annotated_idx].first &&
                        specs[library_idx][peak_idx][1] == ion_mass_intensity_pair_vector[annotated_idx].second){
                        annotated = true;
                        break;
                    }
                }

                if(!annotated){
                    pair<float, float> unannotated_peak;
                    unannotated_peak.first = specs[library_idx][peak_idx][0];
                    unannotated_peak.second = specs[library_idx][peak_idx][1];
                    unannotated_peaks.push_back(unannotated_peak);
                }
            }


            vector<string> original_prefix_array;
            vector<string> original_suffix_array;
            vector<string> random_prefix_array;
            vector<string> random_suffix_array;
            generate_prefix_suffix_peptide(annotation_orig_stripped, original_prefix_array, original_suffix_array);
            generate_prefix_suffix_peptide(random_annotation, random_prefix_array, random_suffix_array);


            AAJumps aajumps(1);

            for(int i = 0; i < ion_mass_intensity_pair_vector.size(); i++){
                if( ions_to_extract[i/(peptide_length)].find('b') != -1 ||
                    ions_to_extract[i/(peptide_length)].find('a') != -1){

                    int charge = 2;
                    vector<float> masses;


                    string orig_prefix = original_prefix_array[i % (peptide_length)];
                    string random_prefix = random_prefix_array[i % (peptide_length)];

                    aajumps.getPRMMasses(create_annotation_ends(orig_prefix).c_str(), masses);
                    float original_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                    aajumps.getPRMMasses(create_annotation_ends(random_prefix).c_str(), masses);
                    float random_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                    float mass_difference = original_mass - random_mass;


                    if(ions_to_extract[i/(peptide_length)].find("++") != -1){
                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                    }
                    else{
                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference*2;
                    }
                }



                if(ions_to_extract[i/(peptide_length)].find('y') != -1){
                    int charge = 2;
                    vector<float> masses;


                    string orig_suffix = original_suffix_array[i % (peptide_length)];
                    string random_suffix = random_suffix_array[i % (peptide_length)];

                    aajumps.getPRMMasses(create_annotation_ends(orig_suffix).c_str(), masses);
                    float original_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                    aajumps.getPRMMasses(create_annotation_ends(random_suffix).c_str(), masses);
                    float random_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                    float mass_difference = original_mass - random_mass;


                    if(ions_to_extract[i/(peptide_length)].find("++") != -1){

                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                    }
                    else{
                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference*2;
                    }
                }

                if(ions_to_extract[i/(peptide_length)].find('P') != -1){
                    ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first;
                }
            }



            Spectrum new_synth_spec;

            //Adding back in the unannotated peaks
            //ion_mass_intensity_pair_vector.clear();
            ion_mass_intensity_pair_vector.insert(ion_mass_intensity_pair_vector.end(), unannotated_peaks.begin(), unannotated_peaks.end());

            //Adding these peaks to a spectrum
            new_synth_spec.resize(0);
            sort( ion_mass_intensity_pair_vector.begin(), ion_mass_intensity_pair_vector.end(), mass_intensity_pair_mass_comp);
            int new_peaklist_size = 0;
            for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                if(ion_mass_intensity_pair_vector[p].second > 0.f){
                    new_peaklist_size++;
                    //cout<<ion_mass_intensity_pair_vector[p].first<<" "<<ion_mass_intensity_pair_vector[p].second<<endl;
                }
            }


            new_synth_spec.resize(new_peaklist_size);
            int new_peaklist_idx = 0;
            for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                if(ion_mass_intensity_pair_vector[p].second > 0.f){
                    new_synth_spec[new_peaklist_idx][0] = ion_mass_intensity_pair_vector[p].first;
                    new_synth_spec[new_peaklist_idx][1] = ion_mass_intensity_pair_vector[p].second;
                    new_peaklist_idx++;
                }
            }


            //adding projection to this
            psmPtr decoy_psm(new PeptideSpectrumMatch);
            decoy_psm->m_annotation = create_annotation_ends(random_annotation);
            decoy_psm->m_charge = specs[library_idx].psmList.front()->m_charge;
            //decoy_psm->m_spectrum = &decoy_library.specs[decoy_library.specs.size()-1];

            new_synth_spec.psmList.push_back(decoy_psm);
            new_synth_spec.parentCharge = specs[library_idx].parentCharge;
            new_synth_spec.parentMass = specs[library_idx].parentMass;
            new_synth_spec.parentMZ = specs[library_idx].parentMZ;
            decoy_library.specs.push_back(new_synth_spec);
            decoy_library.specs[decoy_library.specs.size()-1].psmList.front()->m_spectrum = &decoy_library.specs[decoy_library.specs.size()-1];

        }

        return decoy_library;
    }


    string SpectralLibrary::create_decoy_peptide(string peptide, int charge){
        vector<string> deliminated_aminoacids;
        string stripped_peptide = remove_annotation_ends(peptide);
        int seed = hashpeptide(peptide, charge);
        //cout<<peptide<<"\tSEED\t"<<seed<<endl;
        srand(seed);

        //cout<<stripped_peptide<<endl;

        for(int pepidx = 0; pepidx < stripped_peptide.length(); pepidx++){
            if(stripped_peptide[pepidx] != '(' && stripped_peptide[pepidx] != '['){
                string temp_str = "";
                temp_str += stripped_peptide[pepidx];
                deliminated_aminoacids.push_back(temp_str);
                continue;
            }
            if(stripped_peptide[pepidx] == '('){
                int pepidx_right = pepidx+1;
                bool found_parenthesis = false;
                while(pepidx_right < stripped_peptide.length() ){
                    if(stripped_peptide[pepidx_right] == ')'){
                        found_parenthesis = true;
                        break;
                    }
                    pepidx_right++;
                }

                if(!found_parenthesis){
                    cout<<"Bad Annotation"<<endl;
                    cout<<peptide<<endl;
                    exit(1);
                }

                string temp_str = "";
                for(int i = pepidx; i <= pepidx_right; i++){
                    temp_str += stripped_peptide[i];
                }
                deliminated_aminoacids.push_back(temp_str);
                pepidx = pepidx_right;

                continue;
            }

            if(stripped_peptide[pepidx] == '['){
                int pepidx_right = pepidx+1;
                bool found_bracket = false;
                while(pepidx_right < stripped_peptide.length() ){
                    if(stripped_peptide[pepidx_right] == ']'){
                        found_bracket = true;
                        break;
                    }
                    pepidx_right++;
                }

                if(!found_bracket){
                    cout<<"Bad Annotation"<<endl;
                    cout<<peptide<<endl;
                    exit(1);
                }

                string temp_str = "";
                for(int i = pepidx; i <= pepidx_right; i++){
                    temp_str += stripped_peptide[i];
                }
                deliminated_aminoacids.push_back(temp_str);
                pepidx = pepidx_right;

                continue;
            }
        }

        vector<string> randomized_decoy_array;
        vector<string> remaining_aa_to_assign;
        vector<int> empty_indices;

        //Transferring over immobile amino acids
        for(int i = 0 ; i < deliminated_aminoacids.size(); i++){
            if( deliminated_aminoacids[i].find("R") != string::npos ||
                deliminated_aminoacids[i].find("K") != string::npos ||
                deliminated_aminoacids[i].find("P") != string::npos){

                randomized_decoy_array.push_back(deliminated_aminoacids[i]);
                continue;
            }
            else{
                string temp = "";
                randomized_decoy_array.push_back(temp);
            }
        }


        //Now we clear out R, K, P from original list
        for(int i = 0; i < deliminated_aminoacids.size(); i++){
            if(!( deliminated_aminoacids[i].find("R") != string::npos ||
                deliminated_aminoacids[i].find("K") != string::npos ||
                deliminated_aminoacids[i].find("P") != string::npos)){
                remaining_aa_to_assign.push_back(deliminated_aminoacids[i]);
                empty_indices.push_back(i);
            }
        }

        while(empty_indices.size() > 0){
            int rand_idx = rand()%remaining_aa_to_assign.size();
            randomized_decoy_array[empty_indices[0]] = remaining_aa_to_assign[rand_idx];
            empty_indices.erase(empty_indices.begin());
            remaining_aa_to_assign.erase(remaining_aa_to_assign.begin() + rand_idx);
        }

        string randomized = "";
        for(int i = 0; i < randomized_decoy_array.size(); i++){
            randomized +=randomized_decoy_array[i];
        }

        return randomized;
    }

    void SpectralLibrary::generate_prefix_suffix_peptide(string peptide, vector<string> &prefix, vector<string> &suffix){
        vector<string> deliminated_aminoacids;
        string stripped_peptide = remove_annotation_ends(peptide);

        //cout<<stripped_peptide<<endl;

        for(int pepidx = 0; pepidx < stripped_peptide.length(); pepidx++){
            if(stripped_peptide[pepidx] != '(' && stripped_peptide[pepidx] != '['){
                string temp_str = "";
                temp_str += stripped_peptide[pepidx];
                deliminated_aminoacids.push_back(temp_str);
                continue;
            }
            if(stripped_peptide[pepidx] == '('){
                int pepidx_right = pepidx+1;
                bool found_parenthesis = false;
                while(pepidx_right < stripped_peptide.length() ){
                    if(stripped_peptide[pepidx_right] == ')'){
                        found_parenthesis = true;
                        break;
                    }
                    pepidx_right++;
                }

                if(!found_parenthesis){
                    cout<<"Bad Annotation"<<endl;
                    cout<<peptide<<endl;
                    exit(1);
                }

                string temp_str = "";
                for(int i = pepidx; i <= pepidx_right; i++){
                    temp_str += stripped_peptide[i];
                }
                deliminated_aminoacids.push_back(temp_str);
                pepidx = pepidx_right;

                continue;
            }

            if(stripped_peptide[pepidx] == '['){
                int pepidx_right = pepidx+1;
                bool found_bracket = false;
                while(pepidx_right < stripped_peptide.length() ){
                    if(stripped_peptide[pepidx_right] == ']'){
                        found_bracket = true;
                        break;
                    }
                    pepidx_right++;
                }

                if(!found_bracket){
                    cout<<"Bad Annotation"<<endl;
                    cout<<peptide<<endl;
                    exit(1);
                }

                string temp_str = "";
                for(int i = pepidx; i <= pepidx_right; i++){
                    temp_str += stripped_peptide[i];
                }
                deliminated_aminoacids.push_back(temp_str);
                pepidx = pepidx_right;

                continue;
            }
        }

        prefix.clear();
        suffix.clear();
        for(int i = 0; i < deliminated_aminoacids.size(); i++){
            string prefix_annotation = "";
            string suffix_annotation = "";
            for(int j = 0; j <= i; j++){
                prefix_annotation += deliminated_aminoacids[j];
                suffix_annotation += deliminated_aminoacids[deliminated_aminoacids.size() - i + j - 1];
            }
            prefix.push_back(prefix_annotation);
            suffix.push_back(suffix_annotation);
        }
    }


    int SpectralLibrary::add_update_spectrum_to_Library(Spectrum & addition_spectrum){
        if(addition_spectrum.psmList.size() == 0){
            DEBUG_MSG("psm list is empty");
            return -1;
        }
        if(addition_spectrum.psmList.front()->m_organism.size() == 0){
            DEBUG_VAR(addition_spectrum.psmList.front()->m_organism);
            DEBUG_MSG("organism name is empty");
            return -1;
        }

        bool updated = false;
        for(int i = 0; i < this->specs.size(); i++){
            if(this->specs[i].fileName == addition_spectrum.fileName && this->specs[i].scan == addition_spectrum.scan){
                if(this->specs[i].psmList.front()->m_organism.size() > 0){
                    int most_recent_idx_lib = this->specs[i].psmList.front()->m_organism.size()-1;
                    int most_recent_idx_add = addition_spectrum.psmList.front()->m_organism.size()-1;
                    if( this->specs[i].psmList.front()->m_organism[most_recent_idx_lib] == addition_spectrum.psmList.front()->m_organism[most_recent_idx_add] &&
                        this->specs[i].psmList.front()->m_compound_name[most_recent_idx_lib] == addition_spectrum.psmList.front()->m_compound_name[most_recent_idx_add] &&
                        this->specs[i].psmList.front()->m_smiles[most_recent_idx_lib] == addition_spectrum.psmList.front()->m_smiles[most_recent_idx_add] &&
                        this->specs[i].psmList.front()->m_InChI[most_recent_idx_lib] == addition_spectrum.psmList.front()->m_InChI[most_recent_idx_add] &&
                        this->specs[i].psmList.front()->m_InChI_Aux[most_recent_idx_lib] == addition_spectrum.psmList.front()->m_InChI_Aux[most_recent_idx_add] &&
                        this->specs[i].psmList.front()->m_notes == addition_spectrum.psmList.front()->m_notes){

                        //return i+1;
                        return -2;
                    }
                    else{
                        this->specs[i].psmList.front()->m_submission_metadata.push_back(addition_spectrum.psmList.front()->m_submission_metadata[most_recent_idx_add]);
                        this->specs[i].psmList.front()->m_organism.push_back(addition_spectrum.psmList.front()->m_organism[most_recent_idx_add]);
                        this->specs[i].psmList.front()->m_compound_name.push_back(addition_spectrum.psmList.front()->m_compound_name[most_recent_idx_add]);
                        this->specs[i].psmList.front()->m_smiles.push_back(addition_spectrum.psmList.front()->m_smiles[most_recent_idx_add]);
                        this->specs[i].psmList.front()->m_InChI.push_back(addition_spectrum.psmList.front()->m_InChI[most_recent_idx_add]);
                        this->specs[i].psmList.front()->m_InChI_Aux.push_back(addition_spectrum.psmList.front()->m_InChI_Aux[most_recent_idx_add]);
                        updated = true;
                        return i+1;
                    }
                }
                else{
                    return i+1;
                }
            }
        }

        if (!updated){
            this->specs.push_back(addition_spectrum);
            return specs.size();
        }
        return 0;
    }



     void SpectralLibraryGFSearch::train_distribution(SpectralLibrary &library, MS2ScoringModel model, vector<string> ions_to_extract, string allIons){
        DEBUG_MSG("START Training");

        float annotation_tolerance = 0.1; //Eventually will be a param value

        //Initializing the bands
        int number_bands = 10;
        vector<vector< float > > band_deltas;
        vector< float > ion_deletion;
        vector< float > ion_insertion;
        vector< float > ion_nondeletion;
        vector< float > ion_noninsertion;
        for(int band_idx = 0; band_idx < number_bands; band_idx++){
            vector<float> band_delta;
            band_deltas.push_back(band_delta);
        }

        //Clustering all the same annotation spectra
        map<string, vector<Spectrum *> > same_annotation_clusters;
        for(int spectrum_idx = 0; spectrum_idx < library.size(); spectrum_idx++){
            stringstream ss (stringstream::in | stringstream::out);
            ss<<library[spectrum_idx].parentCharge;
            string annotation = library[spectrum_idx].psmList.front()->m_annotation;
            annotation += ss.str();
            same_annotation_clusters[annotation].push_back(&library[spectrum_idx]);
        }

        map<string, vector<Spectrum *> >::iterator it;
        for ( it=same_annotation_clusters.begin() ; it != same_annotation_clusters.end(); it++ ){
            DEBUG_MSG("TRAININGON\t"<<(*it).first<<"\t"<<(*it).second.size());
            string cluster_annotation = (*it).first;
            int cluster_charge = (*it).second[0]->parentCharge;

            if(((*it).second.size()) < 10) continue;

            AAJumps jumps(1);
            psmPtr training_library_consensus_psm;

            bool found_library_spectrum = false;
            for(int library_idx = 0; library_idx < this->size(); library_idx++){
                if(this->specs[library_idx].psmList.front()->m_annotation == (*it).second[0]->psmList.front()->m_annotation &&
                    this->specs[library_idx].parentCharge == cluster_charge){
                    DEBUG_MSG("FOUND LIBRARY");
                    training_library_consensus_psm = this->specs[library_idx].psmList.front();
                    training_library_consensus_psm->m_spectrum = &this->specs[library_idx];
                    found_library_spectrum = true;
                }
            }

            if(!found_library_spectrum) continue;

            string library_annotation = training_library_consensus_psm->m_annotation;
            training_library_consensus_psm->annotate(training_library_consensus_psm->m_annotation, allIons, model, 0,0, jumps);
            vector<float> library_consensus_ions;

            extractIons(training_library_consensus_psm, create_deliminated_aminoacids(library_annotation).size() ,model,ions_to_extract,library_consensus_ions, 0, 0);
            sqrt_vector(library_consensus_ions);
            normalize_extracted_ions(library_consensus_ions);

            for(int cluster_spec_idx = 0; cluster_spec_idx < (*it).second.size(); cluster_spec_idx++){
                vector<float> library_training_ions;

                float abundance = (*it).second[cluster_spec_idx]->getTotalIonCurrent();

                //if(abundance < 12000.f) continue;
                if(abundance > 12000.f) continue;

                psmPtr training_library_psm = (*it).second[cluster_spec_idx]->psmList.front();
                training_library_psm->annotate(training_library_psm->m_annotation, allIons, model, 0,0, jumps);
                extractIons(training_library_psm, create_deliminated_aminoacids(library_annotation).size() ,model,ions_to_extract,library_training_ions, 0, 0);
                sqrt_vector(library_training_ions);
                normalize_extracted_ions(library_training_ions);

                float dot_product = spectrum_similarity(training_library_psm, training_library_consensus_psm, create_deliminated_aminoacids(library_annotation).size(), model, ions_to_extract, allIons);
                float explained_intensity;
                float library_side_dot_product = spectrum_similarity_sqrt_librarypeaks(training_library_consensus_psm, training_library_psm, create_deliminated_aminoacids(library_annotation).size(), model, ions_to_extract, allIons, explained_intensity);
                float full_spectrum_dot = full_spectrum_similarity(*training_library_psm->m_spectrum, *training_library_consensus_psm->m_spectrum);
                DEBUG_MSG("SIM\t"<<dot_product<<"\tfull spec\t"<<full_spectrum_dot<<"\tLibrarySideSim\t"<<library_side_dot_product<<"\t"<<training_library_psm->m_spectrumFile<<"\t"<<training_library_psm->m_spectrum->scan<<"\t"<<(*training_library_psm->m_spectrum).getTotalIonCurrent()<<"\t"<<cluster_annotation<<"\t"<<training_library_consensus_psm->m_spectrum->parentMZ<<"\t"<<training_library_psm->m_spectrum->parentMZ);


                for(int ion_idx = 0; ion_idx < library_training_ions.size(); ion_idx++){

                    float consensus_ion = library_consensus_ions[ion_idx];
                    float library_ion = library_training_ions[ion_idx];
                    float max_library_ion = max_ion_intensity(library_consensus_ions);

                    if(consensus_ion > 0.0001 && library_ion < 0.0001){
                        //Deletion Probability
                        float percent_intensity = min(library_ion/max_library_ion, 1.f);
                        ion_deletion.push_back(percent_intensity);
                        continue;
                    }
                    if(consensus_ion < 0.0001 && library_ion < 0.0001){
                        //Insertion, not concerned
                        continue;
                    }
                    if(consensus_ion < 0.0001 && library_ion < 0.0001){
                        //Nothing
                        continue;
                    }


                    float ratio = library_ion/consensus_ion;
                    //DEBUG_MSG("library\t"<<library_ion<<"\t"<<consensus_ion);

                    float percent_intensity = min(consensus_ion/max_library_ion, 1.f);
                    ion_nondeletion.push_back(percent_intensity);
                    ion_noninsertion.push_back(percent_intensity);


                    float delta_ion_intensity = (log(library_ion/consensus_ion))/log(2);

                    //Both are present, lets find a bucket for it
                    int band_idx = min(consensus_ion/max_library_ion, 0.999999f)/(1.0/number_bands);
                    if(band_idx >= number_bands){
                        cout<<"BAND IDX ERROR"<<"\t"<<consensus_ion<<"\t"<<max_library_ion<<"\t"<<band_idx<<endl;
                        continue;
                    }

                    if(delta_ion_intensity > 1000) continue;

                    band_deltas[band_idx].push_back(delta_ion_intensity);

                }
            }

            //Skip Quadratic
            continue;

        }


        DEBUG_MSG("CREATING HISTOGRAMS");
        vector<float> band_averages;
        vector<float> band_variances;
        //Training Cosine Distributions
        histograms.clear();
        int buckets = 200;
        float start = -4;
        float end = 4;
        for(int band_idx = 0; band_idx < number_bands; band_idx++){
            vector<pair<float, float> > histogram = create_histogram(buckets, start, end, band_deltas[band_idx], true);
            histograms.push_back(histogram);
        }

        for(int i = 0; i < buckets; i++){
            cout<<histograms[0][i].first<<"\t";
            for(int band_idx = 0; band_idx < number_bands; band_idx++){
                cout<<histograms[band_idx][i].second<<"\t";
            }
            cout<<endl;
        }

        //Calculating Average and Variance
        for(int band_idx = 0; band_idx < number_bands; band_idx++){
            double sum = 0.f;
            for(int index = 0; index < band_deltas[band_idx].size(); index++){
                sum += band_deltas[band_idx][index];
                if(band_deltas[band_idx][index] > 100){
                    cout<<band_deltas[band_idx][index]<<endl;
                }
            }
            double mean = sum/band_deltas[band_idx].size();
            band_averages.push_back(mean);
            cout<<"Band\t"<<band_idx<<"\tAverage\t"<<sum/band_deltas[band_idx].size()<<"\t"<<band_deltas[band_idx].size()<<"\t"<<sum<<endl;
        }


        for(int band_idx = 0; band_idx < number_bands; band_idx++){
            double square_sum = 0.f;
            for(int index = 0; index < band_deltas[band_idx].size(); index++){
                square_sum += band_deltas[band_idx][index] * band_deltas[band_idx][index];
            }
            double variance = square_sum/band_deltas[band_idx].size() - band_averages[band_idx] * band_averages[band_idx];
            band_variances.push_back(variance);
            cout<<"Band\t"<<band_idx<<"\tVariance\t"<<variance<<endl;
        }

        //Training insertion/deletion distributions
        int insertion_deletion_buckets = 11;
        float insertion_delection_start = -0.15f;
        float insertion_delection_end = 0.95;
        vector < pair < float, float > > ion_deletion_histogram = create_histogram(insertion_deletion_buckets, insertion_delection_start, insertion_delection_end, ion_deletion, false);
        vector < pair < float, float > > ion_nondeletion_histogram = create_histogram(insertion_deletion_buckets, insertion_delection_start, insertion_delection_end, ion_nondeletion, false);


        ion_deletion_histogram.erase(ion_deletion_histogram.begin());
        ion_nondeletion_histogram.erase(ion_nondeletion_histogram.begin());

        deletion_probability.clear();

        for(int i = 0; i < ion_deletion_histogram.size(); i++){
            //DEBUG_MSG("DELETION\t"<<ion_deletion_histogram[i].second<<"\t"<<ion_nondeletion_histogram[i].second);

            pair< float, float > deletion_pair;
            deletion_pair.first = ion_deletion_histogram[i].first;
            deletion_pair.second = (ion_deletion_histogram[i].second/(ion_nondeletion_histogram[i].second + ion_deletion_histogram[i].second));
            deletion_probability.push_back(deletion_pair);

            cout<<"INSERTION_DELETION\t"<<deletion_pair.first<<"\t"<<deletion_pair.second<<endl;
        }

        DEBUG_MSG("DONE Training");

    }


    void SpectralLibraryGFSearch::generate_distributions(MS2ScoringModel model, vector<string> ions_to_extract, string allIons){
        DEBUG_MSG("START SLGF Generation");
        SpectralLibrary temp;
        train_library_dynamicprogramming(model, ions_to_extract, allIons, temp);

        return;

    }

    void SpectralLibraryGFSearch::train_library_dynamicprogramming(MS2ScoringModel model, vector<string> ions_to_extract, string allIons, SpectralLibrary &library){
        vector<float> band_variance;
        vector<float> band_mean;


        //Calculating the variance and mean from the band histograms
        for(int band_idx = 0; band_idx < histograms.size(); band_idx++){
            float sum_mass = 0.f;
            float sum_weighted_mass = 0.f;
            for(int idx = 0; idx < histograms[band_idx].size(); idx++){
                sum_mass += histograms[band_idx][idx].second;
                sum_weighted_mass += histograms[band_idx][idx].second * histograms[band_idx][idx].first;
            }
            float mean = sum_weighted_mass/sum_mass;
            band_mean.push_back(mean);
            DEBUG_MSG("BAND\t"<<band_idx<<"\t"<<sum_mass<<"\t"<<sum_weighted_mass);

            float sum_variance = 0.f;
            for(int idx = 0; idx < histograms[band_idx].size(); idx++){
                sum_variance += histograms[band_idx][idx].second * (histograms[band_idx][idx].first - mean)*(histograms[band_idx][idx].first - mean);
            }
            float variance = sum_variance/sum_mass;
            band_variance.push_back(variance);
            DEBUG_MSG("VAR\t"<<band_idx<<"\t"<<variance);
        }





        for(int spec_idx = 0; spec_idx < this->specs.size(); spec_idx++){
        //for(int spec_idx = 0; spec_idx < 100; spec_idx++){
            //int thread_id = omp_get_thread_num();

            //if(thread_id == 0){
            //    cout<<"Spec IDX\t"<<spec_idx*6<<" of "<<this->specs.size()<<endl;
            //}


            vector<pair<float, float> > insertion_probability;
            dynamicprogramming_spectrum(this->specs[spec_idx],
                                        model,
                                        ions_to_extract,
                                        allIons,
                                        band_mean,
                                        band_variance,
                                        insertion_probability,
                                        deletion_probability,
                                        this->specs[spec_idx].psmList.front()->SLGF_distribution,
                                        histograms, 4);
        }
    }

    void SpectralLibraryGFSearch::save_distributions(string output_file_path){
        string output_histogram_path = output_file_path + ".histogram";
        string output_deletion_path = output_file_path + ".deletionprob";
        save_cosine_distribution(histograms, output_histogram_path);
        save_deletion_distribution(deletion_probability, output_deletion_path);
    }

    void SpectralLibraryGFSearch::load_distributions(string input_path){
        string input_histogram_path = input_path + ".histogram";
        string input_deletion_path = input_path + ".deletionprob";
        load_cosine_distribution(histograms, input_histogram_path);
        load_deletion_distribution(deletion_probability, input_deletion_path);
    }

    int SpectralLibrary::search_target_decoy_specset_SLGF(SpectralLibrary &decoy,
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
                                                                         PeptideSpectrumMatchSet &output_psm_set){

        vector<int> accepted_fragmentation;
        accepted_fragmentation.push_back(Spectrum::FragType_CID);

        //Preextracting ions for target and decoy libraries
        for(int lib_idx = 0; lib_idx < target_library_ptr.size(); lib_idx++){
            preprocess_library_ion_extraction(target_library_ptr[lib_idx]->psmList.front(),
                          create_deliminated_aminoacids(target_library_ptr[lib_idx]->psmList.front()->m_annotation).size(),
                          model,
                          ionsToExtract,
                          allIons);
        }

        for(int lib_idx = 0; lib_idx < decoy_library_ptr.size(); lib_idx++){
            preprocess_library_ion_extraction(decoy_library_ptr[lib_idx]->psmList.front(),
                          create_deliminated_aminoacids(decoy_library_ptr[lib_idx]->psmList.front()->m_annotation).size(),
                          model,
                          ionsToExtract,
                          allIons);
        }

        int searched_count = 0;

        #pragma omp parallel for num_threads(8) schedule(guided)
        for(int query_idx = 0; query_idx < searchable_spectra.size(); query_idx++){
            //if(searchable_spectra[query_idx].scan != 9730) continue;
            //if(searchable_spectra[query_idx].scan != 3506) continue;

            //Filtering in acceptable fragmentation types
            bool valid_fragmentation = false;
            for(int fragmentation_idx = 0; fragmentation_idx < accepted_fragmentation.size(); fragmentation_idx++){
                if(searchable_spectra[query_idx].msFragType == accepted_fragmentation[fragmentation_idx]){
                    valid_fragmentation = true;
                    break;
                }
            }

            if(!valid_fragmentation) continue;



            //cout<<"Searching Scan:\t"<<searchable_spectra[query_idx].scan<<"\t";
            //cout<<"mslevel\t"<<searchable_spectra[query_idx].msLevel
            //cout<<searchable_spectra[query_idx].parentMass<<"\t"<<searchable_spectra[query_idx].parentMZ<<"\t"<<searchable_spectra[query_idx].parentCharge<<"\t";
            psmPtr targetdecoy_psm(new PeptideSpectrumMatch);
            int target_decoy_search = this->search_target_decoy_SLGF(decoy,
                                                                searchable_spectra[query_idx],
                                                                targetdecoy_psm,
                                                                parentmz_tolerance,
                                                                target_library_ptr,
                                                                decoy_library_ptr,
                                                                scoring_method,
                                                                model,
                                                                ionsToExtract,
                                                                allIons);


            if(target_decoy_search == 0){
                targetdecoy_psm->m_scanNum = searchable_spectra[query_idx].scan;
                //cout<<"Scan\t"<<targetdecoy_psm->m_scanNum<<"\t"<<"Library IDX\t"<<targetdecoy_psm->m_dbIndex<<"\t"<<targetdecoy_psm->m_annotation<<endl;
                float match_score = targetdecoy_psm->m_score;
                //cout<<"ISDECOY:\t"<<targetdecoy_psm->m_isDecoy<<"\t"<<targetdecoy_psm->m_annotation<<"\t"<<targetdecoy_psm->m_score<<"\t"<<targetdecoy_psm->m_spectrum->scan<<"\t";
                if( ((do_score_threshold == 0)) || match_score > score_threshold){
                    //cout<<match_score<<"\t"<<m_do_score_threshold<<endl;
                    #pragma omp critical
                    {
                        searched_count++;
                        cout<<"Searching Scan:\t"<<searchable_spectra[query_idx].scan<<"\tAND Index\t"<<searched_count<<"\tof\t"<<searchable_spectra.size()<<endl;
                        output_psm_set.push_back(targetdecoy_psm);
                    }
                }

            }
            //cout<<endl;

            continue;
        }

        return 0;
    }


    int SpectralLibrary::search_target_decoy_specset_SLGF(SpectralLibrary &decoy,
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
                                                                         int end_search_idx){



        vector<int> accepted_fragmentation;
        accepted_fragmentation.push_back(Spectrum::FragType_CID);

        //Preextracting ions for target and decoy libraries
        for(int lib_idx = 0; lib_idx < target_library_ptr.size(); lib_idx++){
           //target_library_ptr[lib_idx]->setResolution(1.0005, true);
        }

        for(int lib_idx = 0; lib_idx < decoy_library_ptr.size(); lib_idx++){
            //decoy_library_ptr[lib_idx]->setResolution(1.0005, true);
        }

        //Preextracting ions for target and decoy libraries iso combined
        for(int lib_idx = 0; lib_idx < target_isocombined_library_ptr.size(); lib_idx++){
            //target_isocombined_library_ptr[lib_idx]->setResolution(1.0005, true);
        }

        for(int lib_idx = 0; lib_idx < decoy_isocombined_library_ptr.size(); lib_idx++){
            //decoy_isocombined_library_ptr[lib_idx]->setResolution(1.0005, true);
        }




        //Preextracting ions for target and decoy libraries
        for(int lib_idx = 0; lib_idx < target_library_ptr.size(); lib_idx++){
            preprocess_library_ion_extraction(target_library_ptr[lib_idx]->psmList.front(),
                          create_deliminated_aminoacids(target_library_ptr[lib_idx]->psmList.front()->m_annotation).size(),
                          model,
                          ionsToExtract,
                          allIons);
        }

        for(int lib_idx = 0; lib_idx < decoy_library_ptr.size(); lib_idx++){
            preprocess_library_ion_extraction(decoy_library_ptr[lib_idx]->psmList.front(),
                          create_deliminated_aminoacids(decoy_library_ptr[lib_idx]->psmList.front()->m_annotation).size(),
                          model,
                          ionsToExtract,
                          allIons);
        }

        //Preextracting ions for target and decoy libraries iso combined
        for(int lib_idx = 0; lib_idx < target_isocombined_library_ptr.size(); lib_idx++){
            preprocess_library_ion_extraction_isocombine(target_isocombined_library_ptr[lib_idx]->psmList.front(),
                          create_deliminated_aminoacids(target_isocombined_library_ptr[lib_idx]->psmList.front()->m_annotation).size(),
                          model,
                          ionsToExtract,
                          allIons);
        }

        for(int lib_idx = 0; lib_idx < decoy_isocombined_library_ptr.size(); lib_idx++){
            preprocess_library_ion_extraction_isocombine(decoy_isocombined_library_ptr[lib_idx]->psmList.front(),
                          create_deliminated_aminoacids(decoy_isocombined_library_ptr[lib_idx]->psmList.front()->m_annotation).size(),
                          model,
                          ionsToExtract,
                          allIons);
        }



        int searched_count = 0;

        #pragma omp parallel for num_threads(1) schedule(guided)
        for(int query_idx = start_seach_idx; query_idx < end_search_idx; query_idx++){
            //Filtering on abundance
            int TIC = searchable_spectra[query_idx].getTotalIonCurrent();
            if(abundance == 1){//Low Abundance
                if(TIC > 12000){
                    continue;
                }
            }
            if(abundance == 2){//High Abundance
                if(TIC <= 12000){
                    continue;
                }
            }

            //Filtering in acceptable fragmentation types
            bool valid_fragmentation = false;
            for(int fragmentation_idx = 0; fragmentation_idx < accepted_fragmentation.size(); fragmentation_idx++){
                if(searchable_spectra[query_idx].msFragType == accepted_fragmentation[fragmentation_idx]){
                    valid_fragmentation = true;
                    break;
                }
            }

            if(!valid_fragmentation) continue;


            //cout<<"Searching Scan:\t"<<searchable_spectra[query_idx].scan<<"\t";
            //cout<<"mslevel\t"<<searchable_spectra[query_idx].msLevel
            //cout<<searchable_spectra[query_idx].parentMass<<"\t"<<searchable_spectra[query_idx].parentMZ<<"\t"<<searchable_spectra[query_idx].parentCharge<<"\t";

            /*
            psmPtr targetdecoy_psm(new PeptideSpectrumMatch);



            int target_decoy_search = this->search_target_decoy_SLGF(decoy,
                                                                searchable_spectra[query_idx],
                                                                targetdecoy_psm,
                                                                parentmz_tolerance,
                                                                target_library_ptr,
                                                                decoy_library_ptr,
                                                                target_isocombined,
                                                                decoy_isocombined,
                                                                target_isocombined_library_ptr,
                                                                decoy_isocombined_library_ptr,
                                                                scoring_method,
                                                                model,
                                                                ionsToExtract,
                                                                allIons);




            if(target_decoy_search == 0){
                targetdecoy_psm->m_scanNum = searchable_spectra[query_idx].scan;
                //cout<<"Scan\t"<<targetdecoy_psm->m_scanNum<<"\t"<<"Library IDX\t"<<targetdecoy_psm->m_dbIndex<<"\t"<<targetdecoy_psm->m_annotation<<endl;
                float match_score = targetdecoy_psm->m_score;
                //cout<<"ISDECOY:\t"<<targetdecoy_psm->m_isDecoy<<"\t"<<targetdecoy_psm->m_annotation<<"\t"<<targetdecoy_psm->m_score<<"\t"<<targetdecoy_psm->m_spectrum->scan<<"\t";
                if( ((do_score_threshold == 0)) || match_score > score_threshold){
                    //cout<<match_score<<"\t"<<m_do_score_threshold<<endl;
                    #pragma omp critical
                    {
                        searched_count++;
                        cout<<"Searching Scan:\t"<<searchable_spectra[query_idx].scan<<"\tAND Index\t"<<searched_count<<"\tof\t"<<searchable_spectra.size()<<endl;
                        output_psm_set.push_back(targetdecoy_psm);
                    }
                }

            }*/


           vector<psmPtr> targetdecoy_psms;
            int target_decoy_search = this->search_target_decoy_SLGFNew(decoy,
                                                                searchable_spectra[query_idx],
                                                                targetdecoy_psms,
                                                                output_psm_count,
                                                                parentmz_tolerance,
                                                                target_library_ptr,
                                                                decoy_library_ptr,
                                                                target_isocombined,
                                                                decoy_isocombined,
                                                                target_isocombined_library_ptr,
                                                                decoy_isocombined_library_ptr,
                                                                scoring_method,
                                                                model,
                                                                ionsToExtract,
                                                                allIons,
                                                                abundance);




            if(target_decoy_search == 0){
                for(int i = 0; i < targetdecoy_psms.size(); i++){
                    targetdecoy_psms[i]->m_scanNum = searchable_spectra[query_idx].scan;
                    float match_score = targetdecoy_psms[i]->m_score;
                    if( ((do_score_threshold == 0)) || match_score > score_threshold){
                        #pragma omp critical
                        {
                            searched_count++;
                            cout<<"Searching Scan:\t"<<searchable_spectra[query_idx].scan<<"\tAND Index\t"<<searched_count<<"\tof\t"<<searchable_spectra.size()<<endl;
                            output_psm_set.push_back(targetdecoy_psms[i]);
                        }
                    }
                }
            }



            continue;
        }

        return 0;
    }






    int SpectralLibrary::search_target_decoy_SLGF(SpectralLibrary &decoy,
                                             Spectrum query_spec,
                                             psmPtr output_psm,
                                             float parentmz_tolerance,
                                             vector<Spectrum *> target_library_ptr,
                                             vector<Spectrum *> decoy_library_ptr,
                                             int scoring_method,
                                             MS2ScoringModel &model,
                                             vector<string> &ionsToExtract,
                                             string allIons){
        AAJumps aajumps(1);
        vector<float> masses;
        int charge;

        float query_intensity = query_spec.getTotalIonCurrent();

        PeptideSpectrumMatchSet search_results_decoy;
        vector<score_results_tuple> scores_tuple_decoy;
        //psmPtr psm(new PeptideSpectrumMatch);

        int decoy_start_search_idx;
        int decoy_end_search_idx;
        spectrum_ptr_startend(decoy_library_ptr, query_spec.parentMZ, parentmz_tolerance, decoy_start_search_idx, decoy_end_search_idx);

        //cout<<decoy_start_search_idx<<"\t"<<decoy_end_search_idx<<endl;
        for(int library_idx = decoy_start_search_idx; library_idx <= decoy_end_search_idx; library_idx++){

            float library_mass = decoy_library_ptr[library_idx]->parentMZ;
            charge = decoy_library_ptr[library_idx]->parentCharge;

            if(abs(query_spec.parentMZ - library_mass) > parentmz_tolerance) continue;
            if(query_spec.parentCharge != 0 && query_spec.parentCharge != charge) continue;

            float sim;// = full_spectrum_similarity(*decoy_library_ptr[library_idx], query_spec);
            float dot_bias;// = full_spectrum_dotbias(*decoy_library_ptr[library_idx], query_spec, sim);
            float orig_dot = 0.f;
            float percent_intensity = 0.f;
            float full_cos = 0.f;



            //For rescoring the similarity as a pvalue
            if(scoring_method == SpectralLibrary::MatchScoreType_DotProduct_SLGF){
                psmPtr query_psm(new PeptideSpectrumMatch());
                query_psm->m_spectrum = & query_spec;
                query_psm->m_annotation = decoy_library_ptr[library_idx]->psmList.front()->m_annotation;
                vector<string> deliminated_aminoacids = create_deliminated_aminoacids(decoy_library_ptr[library_idx]->psmList.front()->m_annotation);

                float explained_intensity;
                sim = spectrum_similarity_sqrt_librarypeaks(decoy_library_ptr[library_idx]->psmList.front(), query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons, explained_intensity);

                //sim = sim * 1.01;

                float rescored_sim = SLGF_rescore(decoy_library_ptr[library_idx]->psmList.front()->SLGF_distribution, sim);
                //explained_intensity = percent_explained_intensity(query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons);

                full_cos = sim*explained_intensity;
                orig_dot = sim;
                //sim = rescored_sim*explained_intensity*explained_intensity;
                sim = rescored_sim*explained_intensity;
                //sim = rescored_sim;
                full_cos = sim;


                percent_intensity = explained_intensity;

            }

            score_results_tuple similarity_tuple;
            decoy.specs[library_idx].scan = decoy_library_ptr[library_idx]->scan;
            tr1::get<0>(similarity_tuple) = decoy_library_ptr[library_idx];
            tr1::get<1>(similarity_tuple) = full_cos;
            tr1::get<2>(similarity_tuple) = decoy_library_ptr[library_idx]->psmList.front()->m_dbIndex;
            tr1::get<3>(similarity_tuple) = (string)decoy_library_ptr[library_idx]->psmList.front()->m_annotation;
            tr1::get<4>(similarity_tuple) = dot_bias;
            tr1::get<5>(similarity_tuple) = orig_dot;
            tr1::get<6>(similarity_tuple) = percent_intensity;
            tr1::get<7>(similarity_tuple) = decoy_library_ptr[library_idx]->psmList.front();
            tr1::get<8>(similarity_tuple) = sim;
            scores_tuple_decoy.push_back(similarity_tuple);
        }



        //Finding the best
        //sort(scores.begin(), scores.end(), search_results_comparator);
        sort(scores_tuple_decoy.begin(), scores_tuple_decoy.end(), search_results_comparator);
        for(int i = 0; i < scores_tuple_decoy.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum->scan = psm->m_spectrum->scan;
            psm->m_spectrum = tr1::get<0>(scores_tuple_decoy[i]);
            psm->m_score = tr1::get<1>(scores_tuple_decoy[i]);
            psm->m_dbIndex = tr1::get<2>(scores_tuple_decoy[i]) + 1;
            psm->m_annotation = tr1::get<3>(scores_tuple_decoy[i]);
            psm->m_strict_envelope_score = tr1::get<5>(scores_tuple_decoy[i]);
            psm->m_unstrict_envelope_score = tr1::get<6>(scores_tuple_decoy[i]);
            psm->m_pValue = query_intensity;
            psm->m_fdr = tr1::get<8>(scores_tuple_decoy[i]);;
            search_results_decoy.push_back(psm);
        }




        PeptideSpectrumMatchSet search_results;
        vector<score_results_tuple> scores_tuple;

        int target_start_search_idx;
        int target_end_search_idx;
        spectrum_ptr_startend(target_library_ptr, query_spec.parentMZ, parentmz_tolerance, target_start_search_idx, target_end_search_idx);


        for(int library_idx = target_start_search_idx; library_idx <= target_end_search_idx; library_idx++){

            float library_mass = target_library_ptr[library_idx]->parentMZ;
            charge = target_library_ptr[library_idx]->parentCharge;

            if(abs(query_spec.parentMZ - library_mass) > parentmz_tolerance) continue;
            if(query_spec.parentCharge != 0 && query_spec.parentCharge != charge) continue;
            float sim;// = full_spectrum_similarity(*target_library_ptr[library_idx], query_spec);
            float dot_bias;// = full_spectrum_dotbias(*target_library_ptr[library_idx], query_spec, sim);
            float orig_dot = 0.f;
            float percent_intensity = 0.f;
            float full_cos = 0.f;


            //For rescoring the similarity as a pvalue
            if(scoring_method == SpectralLibrary::MatchScoreType_DotProduct_SLGF){
                psmPtr query_psm(new PeptideSpectrumMatch());
                query_psm->m_spectrum = & query_spec;
                query_psm->m_annotation = target_library_ptr[library_idx]->psmList.front()->m_annotation;
                vector<string> deliminated_aminoacids = create_deliminated_aminoacids(target_library_ptr[library_idx]->psmList.front()->m_annotation);

                float explained_intensity;
                sim = spectrum_similarity_sqrt_librarypeaks(target_library_ptr[library_idx]->psmList.front(), query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons, explained_intensity);

                //sim = sim * 1.01;

                float rescored_sim = SLGF_rescore(target_library_ptr[library_idx]->psmList.front()->SLGF_distribution, sim);
                //explained_intensity = percent_explained_intensity(query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons);

                //DEBUG_MSG("LIBRARY\t"<<target_library_ptr[library_idx]->psmList.front()->m_annotation<<"\t"<<sim<<"\t"<<explained_intensity<<"\t"<<rescored_sim);

                //cout<<rescored_sim<<"\t"<<explained_intensity<<"\t"<<rescored_sim*explained_intensity<<endl;

                orig_dot = sim;
                full_cos = sim*explained_intensity;
                //sim = rescored_sim*explained_intensity*explained_intensity;
                sim = rescored_sim*explained_intensity;
                //sim = rescored_sim;
                full_cos = sim;

                percent_intensity = explained_intensity;


            }
            //DEBUG_MSG(target_library_ptr[library_idx]->psmList.front()->m_dbIndex<<"\t"<<abs(query_spec.parentMZ - library_mass)<<"\t"<<(string)target_library_ptr[library_idx]->psmList.front()->m_annotation<<"\t"<<orig_dot<<"\t"<<percent_intensity<<"\t"<<sim);

            score_results_tuple similarity_tuple;

            tr1::get<0>(similarity_tuple) = target_library_ptr[library_idx];
            tr1::get<1>(similarity_tuple) = full_cos;
            tr1::get<2>(similarity_tuple) = target_library_ptr[library_idx]->psmList.front()->m_dbIndex;
            tr1::get<3>(similarity_tuple) = (string)target_library_ptr[library_idx]->psmList.front()->m_annotation;
            tr1::get<4>(similarity_tuple) = dot_bias;
            tr1::get<5>(similarity_tuple) = orig_dot;
            tr1::get<6>(similarity_tuple) = percent_intensity;
            tr1::get<8>(similarity_tuple) = sim;
            scores_tuple.push_back(similarity_tuple);
        }


        //Finding the best target
        sort(scores_tuple.begin(), scores_tuple.end(), search_results_comparator);
        for(int i = 0; i < scores_tuple.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum = tr1::get<0>(scores_tuple[i]);
            psm->m_score = tr1::get<1>(scores_tuple[i]);
            psm->m_dbIndex = tr1::get<2>(scores_tuple[i]) + 1;
            psm->m_annotation = tr1::get<3>(scores_tuple[i]);
            psm->m_strict_envelope_score = tr1::get<5>(scores_tuple[i]);
            psm->m_unstrict_envelope_score = tr1::get<6>(scores_tuple[i]);
            psm->m_pValue = query_intensity;
            psm->m_fdr = tr1::get<8>(scores_tuple[i]);;
            search_results.push_back(psm);

            //DEBUG_MSG(psm->m_dbIndex<<"\t"<<psm->m_score<<"\t"<<psm->m_annotation);
        }

        float target_top_scoring = 0.f;
        float target_second_scoring = 0.f;
        float decoy_top_scoring = 0.f;
        float decoy_second_scoring = 0.f;

        if(search_results.size() > 0){
            target_top_scoring = search_results[0]->m_score;
        }

        if(search_results.size() > 1){
            target_second_scoring = search_results[1]->m_score;
        }

        if(search_results_decoy.size() > 0){
            decoy_top_scoring = search_results_decoy[0]->m_score;
        }

        if(search_results_decoy.size() > 1){
            decoy_second_scoring = search_results_decoy[1]->m_score;
        }

        //DEBUG_MSG("TARGET\t"<<target_top_scoring<<"\t"<<search_results[0]->m_strict_envelope_score<<"\t"<<search_results[0]->m_unstrict_envelope_score<<"\t"<<search_results[0]->m_annotation);
        //DEBUG_MSG("DECOY\t"<<decoy_top_scoring<<"\t"<<search_results_decoy[0]->m_strict_envelope_score<<"\t"<<search_results_decoy[0]->m_unstrict_envelope_score<<"\t"<<search_results_decoy[0]->m_annotation);

        if(target_top_scoring >= decoy_top_scoring && search_results.size() > 0){
            float dot_bias = tr1::get<4>(scores_tuple[0]);
            float dot_product = target_top_scoring;
            float deltaD = target_top_scoring - max(target_second_scoring, decoy_top_scoring);
            float match_score = target_top_scoring * 0.6 + deltaD * 0.4 - dot_bias;
            float orig_cos = search_results[0]->m_strict_envelope_score;
            float explained_int = search_results[0]->m_unstrict_envelope_score;

            //cout<<"Dot Bias\t"<<dot_bias<<"\t"<<search_results[0]->m_annotation<<"\tdot\t"<<dot_product<<"\tdeltaD\t"<<deltaD<<"\t"<<orig_cos<<endl;


            switch(scoring_method){
                case SpectralLibrary::MatchScoreType_DotProduct:
                    match_score = dot_product;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD:
                    match_score = dot_product * 0.6 + deltaD * 0.4;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD_DotBias:
                    match_score = dot_product * 0.6 + deltaD * 0.4 - dot_bias;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_SLGF:
                    match_score = dot_product;
                    break;
                default:
                    break;
            }


            output_psm->m_spectrum      = search_results[0]->m_spectrum;
            output_psm->m_score         = search_results[0]->m_score;
            output_psm->m_dbIndex       = search_results[0]->m_dbIndex;
            output_psm->m_annotation    = search_results[0]->m_annotation;
            output_psm->m_score         = match_score;
            output_psm->m_strict_envelope_score =  orig_cos;
            output_psm->m_unstrict_envelope_score =  explained_int;
            output_psm->m_pValue = query_intensity;
            output_psm->m_isDecoy = false;
            output_psm->m_fdr = search_results[0]->m_fdr;
            output_psm->m_startMass = deltaD;
            output_psm->m_charge = query_spec.parentCharge;
            output_psm->m_mz = query_spec.parentMZ;

            return 0;
        }
        if(target_top_scoring < decoy_top_scoring && search_results_decoy.size() > 0){
            float dot_bias = tr1::get<4>(scores_tuple_decoy[0]);
            float dot_product = decoy_top_scoring;
            float deltaD = decoy_top_scoring - max(target_top_scoring, decoy_second_scoring);
            float match_score = decoy_top_scoring * 0.6 + deltaD * 0.4 - dot_bias;
            float orig_cos = search_results_decoy[0]->m_strict_envelope_score;
            float explained_int = search_results_decoy[0]->m_unstrict_envelope_score;

            switch(scoring_method){
                case SpectralLibrary::MatchScoreType_DotProduct:
                    match_score = dot_product;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD:
                    match_score = dot_product * 0.6 + deltaD * 0.4;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD_DotBias:
                    match_score = dot_product * 0.6 + deltaD * 0.4 - dot_bias;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_SLGF:
                    match_score = dot_product;
                    break;
                default:
                    break;
            }

            output_psm->m_spectrum      = search_results_decoy[0]->m_spectrum;
            output_psm->m_score         = search_results_decoy[0]->m_score;
            output_psm->m_dbIndex       = search_results_decoy[0]->m_dbIndex;
            output_psm->m_annotation    = search_results_decoy[0]->m_annotation;
            output_psm->m_score         = match_score;
            output_psm->m_strict_envelope_score =  orig_cos;
            output_psm->m_unstrict_envelope_score =  explained_int;
            output_psm->m_pValue = query_intensity;
            output_psm->m_isDecoy = true;
            output_psm->m_fdr = search_results_decoy[0]->m_fdr;
            output_psm->m_startMass = deltaD;
            output_psm->m_charge = query_spec.parentCharge;
            output_psm->m_mz = query_spec.parentMZ;

            return 0;
        }

        return -1;
    }

    //Searching with two models
    int SpectralLibrary::search_target_decoy_SLGF(SpectralLibrary &decoy,
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
                                             string allIons){
        AAJumps aajumps(1);
        vector<float> masses;
        int charge;

        float query_intensity = query_spec.getTotalIonCurrent();

        PeptideSpectrumMatchSet search_results_decoy;
        vector<score_results_tuple> scores_tuple_decoy;
        //psmPtr psm(new PeptideSpectrumMatch);

        float isotopic_mz_error_threshold = 1.2;

        int decoy_start_search_idx;
        int decoy_end_search_idx;
        spectrum_ptr_startend(decoy_library_ptr, query_spec.parentMZ, parentmz_tolerance, decoy_start_search_idx, decoy_end_search_idx);

        //cout<<decoy_start_search_idx<<"\t"<<decoy_end_search_idx<<endl;
        for(int library_idx = decoy_start_search_idx; library_idx <= decoy_end_search_idx; library_idx++){
            //DEBUG_MSG(decoy_library_ptr[library_idx]->parentMZ<<"\t"<<decoy_isocombined_library_ptr[library_idx]->parentMZ);

            float library_mass = decoy_library_ptr[library_idx]->parentMZ;
            charge = decoy_library_ptr[library_idx]->parentCharge;
            float mz_error = abs(query_spec.parentMZ - library_mass);

            if(abs(query_spec.parentMZ - library_mass) > parentmz_tolerance) continue;
            if(query_spec.parentCharge != 0 && query_spec.parentCharge != charge) continue;

            float sim;// = full_spectrum_similarity(*decoy_library_ptr[library_idx], query_spec);
            float dot_bias;// = full_spectrum_dotbias(*decoy_library_ptr[library_idx], query_spec, sim);
            float orig_dot = 0.f;
            float percent_intensity = 0.f;
            float full_cos = 0.f;

            //For rescoring the similarity as a pvalue
            if(scoring_method == SpectralLibrary::MatchScoreType_DotProduct_SLGF){
                psmPtr query_psm(new PeptideSpectrumMatch());
                query_psm->m_spectrum = & query_spec;
                query_psm->m_annotation = decoy_library_ptr[library_idx]->psmList.front()->m_annotation;
                vector<string> deliminated_aminoacids = create_deliminated_aminoacids(decoy_library_ptr[library_idx]->psmList.front()->m_annotation);



                float explained_intensity;
                sim = spectrum_similarity_sqrt_librarypeaks(decoy_library_ptr[library_idx]->psmList.front(), query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons,explained_intensity);


                //sim = sim * 1.01;

                float rescored_sim = SLGF_rescore(decoy_library_ptr[library_idx]->psmList.front()->SLGF_distribution, sim);
                if(mz_error > isotopic_mz_error_threshold){
                    sim = spectrum_similarity_sqrt_librarypeaks_isocombine(decoy_isocombined_library_ptr[library_idx]->psmList.front(), query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons, explained_intensity);
                    rescored_sim = SLGF_rescore(decoy_isocombined_library_ptr[library_idx]->psmList.front()->SLGF_distribution, sim);
                }



                //explained_intensity = percent_explained_intensity(query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons);
                full_cos = sim*explained_intensity;
                orig_dot = sim;
                //sim = rescored_sim*explained_intensity*explained_intensity;
                sim = rescored_sim*explained_intensity;
                sim = rescored_sim;
                full_cos = sim;

                percent_intensity = explained_intensity;

            }

            score_results_tuple similarity_tuple;
            decoy.specs[library_idx].scan = decoy_library_ptr[library_idx]->scan;
            tr1::get<0>(similarity_tuple) = decoy_library_ptr[library_idx];
            tr1::get<1>(similarity_tuple) = full_cos;
            tr1::get<2>(similarity_tuple) = decoy_library_ptr[library_idx]->psmList.front()->m_dbIndex;
            tr1::get<3>(similarity_tuple) = (string)decoy_library_ptr[library_idx]->psmList.front()->m_annotation;
            tr1::get<4>(similarity_tuple) = dot_bias;
            tr1::get<5>(similarity_tuple) = orig_dot;
            tr1::get<6>(similarity_tuple) = percent_intensity;
            tr1::get<7>(similarity_tuple) = decoy_library_ptr[library_idx]->psmList.front();
            tr1::get<8>(similarity_tuple) = sim;
            scores_tuple_decoy.push_back(similarity_tuple);
        }



        //Finding the best
        //sort(scores.begin(), scores.end(), search_results_comparator);
        sort(scores_tuple_decoy.begin(), scores_tuple_decoy.end(), search_results_comparator);
        for(int i = 0; i < scores_tuple_decoy.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum->scan = psm->m_spectrum->scan;
            psm->m_spectrum = tr1::get<0>(scores_tuple_decoy[i]);
            psm->m_score = tr1::get<1>(scores_tuple_decoy[i]);
            psm->m_dbIndex = tr1::get<2>(scores_tuple_decoy[i]) + 1;
            psm->m_annotation = tr1::get<3>(scores_tuple_decoy[i]);
            psm->m_strict_envelope_score = tr1::get<5>(scores_tuple_decoy[i]);
            psm->m_unstrict_envelope_score = tr1::get<6>(scores_tuple_decoy[i]);
            psm->m_pValue = query_intensity;
            psm->m_fdr = tr1::get<8>(scores_tuple_decoy[i]);;
            search_results_decoy.push_back(psm);
        }




        PeptideSpectrumMatchSet search_results;
        vector<score_results_tuple> scores_tuple;

        int target_start_search_idx;
        int target_end_search_idx;
        spectrum_ptr_startend(target_library_ptr, query_spec.parentMZ, parentmz_tolerance, target_start_search_idx, target_end_search_idx);



        for(int library_idx = target_start_search_idx; library_idx <= target_end_search_idx; library_idx++){

            float library_mass = target_library_ptr[library_idx]->parentMZ;
            charge = target_library_ptr[library_idx]->parentCharge;
            float mz_error = abs(query_spec.parentMZ - library_mass);

            if(abs(query_spec.parentMZ - library_mass) > parentmz_tolerance) continue;
            if(query_spec.parentCharge != 0 && query_spec.parentCharge != charge) continue;
            float sim;// = full_spectrum_similarity(*target_library_ptr[library_idx], query_spec);
            float dot_bias;// = full_spectrum_dotbias(*target_library_ptr[library_idx], query_spec, sim);
            float orig_dot = 0.f;
            float percent_intensity = 0.f;
            float full_cos = 0.f;

            //For rescoring the similarity as a pvalue
            if(scoring_method == SpectralLibrary::MatchScoreType_DotProduct_SLGF){
                psmPtr query_psm(new PeptideSpectrumMatch());
                query_psm->m_spectrum = & query_spec;
                query_psm->m_annotation = target_library_ptr[library_idx]->psmList.front()->m_annotation;
                vector<string> deliminated_aminoacids = create_deliminated_aminoacids(target_library_ptr[library_idx]->psmList.front()->m_annotation);

                float explained_intensity;
                sim = spectrum_similarity_sqrt_librarypeaks(target_library_ptr[library_idx]->psmList.front(), query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons, explained_intensity);

                //sim = sim * 1.01;

                float rescored_sim = SLGF_rescore(target_library_ptr[library_idx]->psmList.front()->SLGF_distribution, sim);
                if(mz_error > isotopic_mz_error_threshold){
                    sim = spectrum_similarity_sqrt_librarypeaks_isocombine(target_isocombined_library_ptr[library_idx]->psmList.front(), query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons, explained_intensity);
                    rescored_sim = SLGF_rescore(target_isocombined_library_ptr[library_idx]->psmList.front()->SLGF_distribution, sim);
                }

                //explained_intensity = percent_explained_intensity(query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons);
                //cout<<rescored_sim<<"\t"<<explained_intensity<<"\t"<<rescored_sim*explained_intensity<<endl;

                orig_dot = sim;
                full_cos = sim*explained_intensity;
                //sim = rescored_sim*explained_intensity*explained_intensity;
                sim = rescored_sim*explained_intensity;
                sim = rescored_sim;
                full_cos = sim;

                percent_intensity = explained_intensity;


            }
            //DEBUG_MSG(target_library_ptr[library_idx]->psmList.front()->m_dbIndex<<"\t"<<abs(query_spec.parentMZ - library_mass)<<"\t"<<(string)target_library_ptr[library_idx]->psmList.front()->m_annotation<<"\t"<<orig_dot<<"\t"<<percent_intensity<<"\t"<<sim);

            score_results_tuple similarity_tuple;

            tr1::get<0>(similarity_tuple) = target_library_ptr[library_idx];
            tr1::get<1>(similarity_tuple) = full_cos;
            tr1::get<2>(similarity_tuple) = target_library_ptr[library_idx]->psmList.front()->m_dbIndex;
            tr1::get<3>(similarity_tuple) = (string)target_library_ptr[library_idx]->psmList.front()->m_annotation;
            tr1::get<4>(similarity_tuple) = dot_bias;
            tr1::get<5>(similarity_tuple) = orig_dot;
            tr1::get<6>(similarity_tuple) = percent_intensity;
            tr1::get<8>(similarity_tuple) = sim;
            scores_tuple.push_back(similarity_tuple);
        }


        //Finding the best target
        sort(scores_tuple.begin(), scores_tuple.end(), search_results_comparator);
        for(int i = 0; i < scores_tuple.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum = tr1::get<0>(scores_tuple[i]);
            psm->m_score = tr1::get<1>(scores_tuple[i]);
            psm->m_dbIndex = tr1::get<2>(scores_tuple[i]) + 1;
            psm->m_annotation = tr1::get<3>(scores_tuple[i]);
            psm->m_strict_envelope_score = tr1::get<5>(scores_tuple[i]);
            psm->m_unstrict_envelope_score = tr1::get<6>(scores_tuple[i]);
            psm->m_pValue = query_intensity;
            psm->m_fdr = tr1::get<8>(scores_tuple[i]);;
            search_results.push_back(psm);

            //DEBUG_MSG(psm->m_dbIndex<<"\t"<<psm->m_score<<"\t"<<psm->m_annotation);
        }

        float target_top_scoring = 0.f;
        float target_second_scoring = 0.f;
        float decoy_top_scoring = 0.f;
        float decoy_second_scoring = 0.f;

        if(search_results.size() > 0){
            target_top_scoring = search_results[0]->m_score;
        }

        if(search_results.size() > 1){
            target_second_scoring = search_results[1]->m_score;
        }

        if(search_results_decoy.size() > 0){
            decoy_top_scoring = search_results_decoy[0]->m_score;
        }

        if(search_results_decoy.size() > 1){
            decoy_second_scoring = search_results_decoy[1]->m_score;
        }

        //DEBUG_MSG("TARGET\t"<<target_top_scoring<<"\t"<<search_results[0]->m_strict_envelope_score<<"\t"<<search_results[0]->m_unstrict_envelope_score<<"\t"<<search_results[0]->m_annotation);
        //DEBUG_MSG("DECOY\t"<<decoy_top_scoring<<"\t"<<search_results_decoy[0]->m_strict_envelope_score<<"\t"<<search_results_decoy[0]->m_unstrict_envelope_score<<"\t"<<search_results_decoy[0]->m_annotation);

        if(target_top_scoring >= decoy_top_scoring && search_results.size() > 0){
            float dot_bias = tr1::get<4>(scores_tuple[0]);
            float dot_product = target_top_scoring;
            float deltaD = target_top_scoring - max(target_second_scoring, decoy_top_scoring);
            float match_score = target_top_scoring * 0.6 + deltaD * 0.4 - dot_bias;
            float orig_cos = search_results[0]->m_strict_envelope_score;
            float explained_int = search_results[0]->m_unstrict_envelope_score;

            //cout<<"Dot Bias\t"<<dot_bias<<"\t"<<search_results[0]->m_annotation<<"\tdot\t"<<dot_product<<"\tdeltaD\t"<<deltaD<<"\t"<<orig_cos<<endl;

            switch(scoring_method){
                case SpectralLibrary::MatchScoreType_DotProduct:
                    match_score = dot_product;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD:
                    match_score = dot_product * 0.6 + deltaD * 0.4;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD_DotBias:
                    match_score = dot_product * 0.6 + deltaD * 0.4 - dot_bias;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_SLGF:
                    match_score = dot_product;
                    break;
                default:
                    break;
            }


            output_psm->m_spectrum      = search_results[0]->m_spectrum;
            output_psm->m_score         = search_results[0]->m_score;
            output_psm->m_dbIndex       = search_results[0]->m_dbIndex;
            output_psm->m_annotation    = search_results[0]->m_annotation;
            output_psm->m_score         = match_score;
            output_psm->m_strict_envelope_score =  orig_cos;
            output_psm->m_unstrict_envelope_score =  explained_int;
            output_psm->m_pValue = query_intensity;
            output_psm->m_isDecoy = false;
            output_psm->m_fdr = search_results[0]->m_fdr;
            output_psm->m_startMass = deltaD;

            return 0;
        }
        if(target_top_scoring < decoy_top_scoring && search_results_decoy.size() > 0){
            float dot_bias = tr1::get<4>(scores_tuple_decoy[0]);
            float dot_product = decoy_top_scoring;
            float deltaD = decoy_top_scoring - max(target_top_scoring, decoy_second_scoring);
            float match_score = decoy_top_scoring * 0.6 + deltaD * 0.4 - dot_bias;
            float orig_cos = search_results_decoy[0]->m_strict_envelope_score;
            float explained_int = search_results_decoy[0]->m_unstrict_envelope_score;

            switch(scoring_method){
                case SpectralLibrary::MatchScoreType_DotProduct:
                    match_score = dot_product;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD:
                    match_score = dot_product * 0.6 + deltaD * 0.4;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD_DotBias:
                    match_score = dot_product * 0.6 + deltaD * 0.4 - dot_bias;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_SLGF:
                    match_score = dot_product;
                    break;
                default:
                    break;
            }

            output_psm->m_spectrum      = search_results_decoy[0]->m_spectrum;
            output_psm->m_score         = search_results_decoy[0]->m_score;
            output_psm->m_dbIndex       = search_results_decoy[0]->m_dbIndex;
            output_psm->m_annotation    = search_results_decoy[0]->m_annotation;
            output_psm->m_score         = match_score;
            output_psm->m_strict_envelope_score =  orig_cos;
            output_psm->m_unstrict_envelope_score =  explained_int;
            output_psm->m_pValue = query_intensity;
            output_psm->m_isDecoy = true;
            output_psm->m_fdr = search_results_decoy[0]->m_fdr;
            output_psm->m_startMass = deltaD;

            return 0;
        }

        return -1;
    }


}





int SpectralLibrary::search_target_decoy_SLGFNew(SpectralLibrary &decoy,
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
                                             int abundance){

    DEBUG_MSG("SEARCHING SCAN "<<query_spec.scan);

    vector<psmPtr> search_results;

    float ISOTOPIC_MZ_ERROR_THRESHOLD = 1.2;

    float query_intensity = query_spec.getTotalIonCurrent();
    psmPtr query_psm(new PeptideSpectrumMatch());
    query_psm->m_spectrum = & query_spec;

    for(int library_idx = 0; library_idx < this->size(); library_idx++){
        float library_mass = specs[library_idx].parentMZ;
        float query_mass = query_spec.parentMZ;

        int charge = specs[library_idx].parentCharge;
        int query_charge = query_spec.parentCharge;

        if(charge != query_charge) continue;


        float sim;
        float percent_intensity = 0.f;

        query_psm->m_annotation = specs[library_idx].psmList.front()->m_annotation;

        float mass_difference = (query_mass - library_mass) * charge;
        float mz_difference = (query_mass - library_mass);
        if(abs(mz_difference) > parentmz_tolerance) continue;

        vector<string> deliminated_peptide = create_deliminated_aminoacids(query_psm->m_annotation);

        float explained_intensity;
        float rescored_sim;

        string library_name = "";

        if(abundance == 1 || (abundance == 2 && mz_difference < ISOTOPIC_MZ_ERROR_THRESHOLD) || abundance == 0){
            sim = spectrum_similarity_sqrt_librarypeaks(specs[library_idx].psmList.front(), query_psm, deliminated_peptide.size(), model, ionsToExtract, allIons,explained_intensity);
            rescored_sim = SLGF_rescore(specs[library_idx].psmList.front()->SLGF_distribution, sim);
            //library_name = specs[library_idx].fileName;
            library_name = specs[library_idx].psmList.front()->m_spectrumFile;
        }

        if(abundance == 2 && mz_difference > ISOTOPIC_MZ_ERROR_THRESHOLD){
            sim = spectrum_similarity_sqrt_librarypeaks_isocombine(target_isocombined[library_idx].psmList.front(), query_psm, deliminated_peptide.size(), model, ionsToExtract, allIons, explained_intensity);
            rescored_sim = SLGF_rescore(target_isocombined[library_idx].psmList.front()->SLGF_distribution, sim);
            //library_name = target_isocombined[library_idx].fileName;
            library_name = target_isocombined[library_idx].psmList.front()->m_spectrumFile;
        }

        float final_score = rescored_sim*explained_intensity;


        psmPtr search_result(new PeptideSpectrumMatch());
        search_result->m_annotation = query_psm->m_annotation;
        search_result->m_spectrum = & query_spec;
        search_result->m_isDecoy = false;
        search_result->m_score = final_score;
        search_result->m_pValue = sim;
        search_result->m_spectrumFile = query_spec.fileName;
        search_result->m_dbIndex = library_idx + 1;
        search_result->m_library_name = get_only_filename(library_name);



        search_results.push_back(search_result);

        //DEBUG_MSG("TARGET\t"<<final_score);
    }

    for(int decoy_idx = 0; decoy_idx < decoy.size(); decoy_idx++){
        float library_mass = decoy[decoy_idx].parentMZ;
        float query_mass = query_spec.parentMZ;

        int charge = decoy[decoy_idx].parentCharge;
        int query_charge = query_spec.parentCharge;

        if(charge != query_charge) continue;


        float sim;
        float percent_intensity = 0.f;


        query_psm->m_annotation = decoy[decoy_idx].psmList.front()->m_annotation;

        float mass_difference = (query_mass - library_mass) * charge;
        float mz_difference = (query_mass - library_mass);
        if(abs(mz_difference) > parentmz_tolerance) continue;

        vector<string> deliminated_peptide = create_deliminated_aminoacids(query_psm->m_annotation);

        float explained_intensity;
        float rescored_sim;
        rescored_sim = SLGF_rescore(decoy[decoy_idx].psmList.front()->SLGF_distribution, sim);

        string library_name = "";

        if(abundance == 1 || (abundance == 2 && mz_difference < ISOTOPIC_MZ_ERROR_THRESHOLD) || abundance == 0){
            sim = spectrum_similarity_sqrt_librarypeaks(decoy[decoy_idx].psmList.front(), query_psm, deliminated_peptide.size(), model, ionsToExtract, allIons,explained_intensity);
            rescored_sim = SLGF_rescore(decoy[decoy_idx].psmList.front()->SLGF_distribution, sim);
            //library_name = decoy[decoy_idx].fileName;
            library_name = decoy[decoy_idx].psmList.front()->m_spectrumFile;
        }

        if(abundance == 2 && mz_difference > ISOTOPIC_MZ_ERROR_THRESHOLD){
            sim = spectrum_similarity_sqrt_librarypeaks_isocombine(decoy_isocombined[decoy_idx].psmList.front(), query_psm, deliminated_peptide.size(), model, ionsToExtract, allIons, explained_intensity);
            rescored_sim = SLGF_rescore(decoy_isocombined[decoy_idx].psmList.front()->SLGF_distribution, sim);
            //library_name = decoy_isocombined[decoy_idx].fileName;
            library_name = decoy_isocombined[decoy_idx].psmList.front()->m_spectrumFile;
        }

        float final_score = rescored_sim*explained_intensity;

        psmPtr search_result(new PeptideSpectrumMatch());
        search_result->m_annotation = query_psm->m_annotation;
        search_result->m_spectrum = & query_spec;
        search_result->m_isDecoy = true;
        search_result->m_score = final_score;
        search_result->m_pValue = sim;
        search_result->m_spectrumFile = query_spec.fileName;
        search_result->m_dbIndex = decoy_idx + 1;
        search_result->m_library_name = get_only_filename(library_name);


        search_results.push_back(search_result);

        //DEBUG_MSG("DECOY\t"<<final_score);
    }

    sort(search_results.begin(), search_results.end(), search_results_comparator_psmPtr);
    if(search_results.size() > 0){
        for(int i = 0; i < top_psm_number; i++){
            if(search_results.size() <= i) break;
            output_psms.push_back(search_results[i]);
        }
        return 0;
    }


    return -1;
}
