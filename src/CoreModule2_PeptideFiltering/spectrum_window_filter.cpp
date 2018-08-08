#include "spectrum_window_filter.h"


bool mass_intensity_pair_comp(pair < float, float > mass_intensity_1, pair < float, float > mass_intensity_2){
    return (mass_intensity_1.second > mass_intensity_2.second);
}

bool mass_intensity_pair_mass_comp(pair < float, float > mass_intensity_1, pair < float, float > mass_intensity_2){
    return (mass_intensity_1.first < mass_intensity_2.first);
}

int filter_window(Spectrum * s, int window, int max_peaks){
    float max_mass = 0.f;
    for(int peak_idx = 0; peak_idx < s->size(); peak_idx++){
        if( (*s)[peak_idx][0] > max_mass){
            max_mass = (*s)[peak_idx][0];
        }
    }
    //cout<<"Window Peaks: "<<s->peakList.size()<<endl;
    
    vector< vector < pair < float, float > > > buckets;
    int number_buckets = (max_mass)/window + 1;
    //cout<<"number buckets: "<<number_buckets<<endl;
    
    for(int bucket_idx = 0; bucket_idx < number_buckets; bucket_idx++){
        vector < pair < float, float > > temp_bucket;
        buckets.push_back(temp_bucket);
    }
    
    //cout<<"added buckets"<<endl;
    
    for(int peak_idx = 0; peak_idx < s->size(); peak_idx++){
        int bucket = (*s)[peak_idx][0]/window;
        pair <float, float> mass_intensity_pair;
        mass_intensity_pair.first = (*s)[peak_idx][0];
        mass_intensity_pair.second = (*s)[peak_idx][1];
        buckets[bucket].push_back(mass_intensity_pair);
    }
    
    //cout<<"put stuff into buckets"<<endl;
    
    for(int bucket_idx = 0; bucket_idx < number_buckets; bucket_idx++){
        sort(buckets[bucket_idx].begin(), buckets[bucket_idx].end(), mass_intensity_pair_comp);
    }
    
    //cout<<"sorted buckets"<<endl;
    vector< pair < float, float > > new_spec_list;
    
    for(int bucket_idx = 0; bucket_idx < number_buckets; bucket_idx++){
        //for(int i = 0; i < min(max_peaks,(int)buckets[bucket_idx].size()); i++){
        for(int i = 0; i < buckets[bucket_idx].size(); i++){
            //cout<<"Mass: "<<buckets[bucket_idx][i].first<<" Intensity: "<<buckets[bucket_idx][i].second<<endl;
            pair < float, float > mass_intensity_pair;
            mass_intensity_pair.first = buckets[bucket_idx][i].first;
            mass_intensity_pair.second = buckets[bucket_idx][i].second;
            if(i >= max_peaks){
                mass_intensity_pair.second /= 1000;
            }
            else{
                new_spec_list.push_back(mass_intensity_pair);
            }
        }
    }
    
    sort( new_spec_list.begin(), new_spec_list.end(), mass_intensity_pair_mass_comp);
    
    s->resize(new_spec_list.size());
    for(int peak_idx = 0; peak_idx < new_spec_list.size(); peak_idx++){
        (*s)[peak_idx][0] = new_spec_list[peak_idx].first;
        (*s)[peak_idx][1] = new_spec_list[peak_idx].second;
    }
    
    //cout<<"Window Peaks: "<<s->peakList.size()<<endl;
    
    return 0;
}
