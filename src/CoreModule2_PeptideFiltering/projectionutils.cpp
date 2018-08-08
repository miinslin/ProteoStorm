#include "projectionutils.h"

float spectrum_similarity(psmPtr psm1, string annotation1,  psmPtr psm2, string annotation2,
                          int peptideLength, MS2ScoringModel &model,
                          vector<string> &ionsToExtract, string allIons){
    
    psmPtr psm1temp(new PeptideSpectrumMatch);
    psmPtr psm2temp(new PeptideSpectrumMatch);
    
    Spectrum temp_spectrum1 = *(psm1->m_spectrum);
    Spectrum temp_spectrum2 = *(psm2->m_spectrum);
    
    Spectrum * old_spectrum1 = psm1->m_spectrum;
    Spectrum * old_spectrum2 = psm2->m_spectrum;
    
    psm1temp->m_spectrum = &temp_spectrum1;
    psm2temp->m_spectrum = &temp_spectrum2;
    psm1temp->m_annotation = psm1->m_annotation;
    psm2temp->m_annotation = psm2->m_annotation;
    
    //Making max 1000
    preprocess_spectrum_intensities_max_intensity(psm1temp->m_spectrum, 1000.f);
    preprocess_spectrum_intensities_max_intensity(psm2temp->m_spectrum, 1000.f);
    
    //Applying SQRT
    preprocess_spectrum_intensities(psm1temp->m_spectrum, 0, 1);
    preprocess_spectrum_intensities(psm2temp->m_spectrum, 0, 1);
    
    psm1temp->annotate(psm1temp->m_annotation,allIons,model,0,0,.45);
    psm2temp->annotate(psm2temp->m_annotation,allIons,model,0,0,.45);
    
    vector<float> ion_intensities1;
    vector<float> ion_intensities2;
    
    extractIons(psm1temp,peptideLength,model,ionsToExtract,ion_intensities1, 0, 0);
    extractIons(psm2temp,peptideLength,model,ionsToExtract,ion_intensities2, 0, 0);

    normalize_extracted_ions(ion_intensities1);
    normalize_extracted_ions(ion_intensities2);
    
    float cosine_val = cosine(ion_intensities1, ion_intensities2);
    
    psm1->m_spectrum = old_spectrum1;
    psm2->m_spectrum = old_spectrum2;
    
    return cosine_val;
}



float spectrum_similarity(psmPtr psm1, psmPtr psm2, int peptideLength, MS2ScoringModel &model,
        vector<string> &ionsToExtract, string allIons){

    psmPtr psm1temp(new PeptideSpectrumMatch);
    psmPtr psm2temp(new PeptideSpectrumMatch);
    
    Spectrum temp_spectrum1 = *(psm1->m_spectrum);
    Spectrum temp_spectrum2 = *(psm2->m_spectrum);
    
    Spectrum * old_spectrum1 = psm1->m_spectrum;
    Spectrum * old_spectrum2 = psm2->m_spectrum;
    
    psm1temp->m_spectrum = &temp_spectrum1;
    psm2temp->m_spectrum = &temp_spectrum2;
    psm1temp->m_annotation = psm1->m_annotation;
    psm2temp->m_annotation = psm2->m_annotation;
    
    //Making max 1000
    preprocess_spectrum_intensities_max_intensity(psm1temp->m_spectrum, 1000.f);
    preprocess_spectrum_intensities_max_intensity(psm2temp->m_spectrum, 1000.f);
    
    //Applying SQRT
    preprocess_spectrum_intensities(psm1temp->m_spectrum, 0, 1);
    preprocess_spectrum_intensities(psm2temp->m_spectrum, 0, 1);
    
    psm1temp->annotate(psm1temp->m_annotation,allIons,model,0,0,.45);
    psm2temp->annotate(psm2temp->m_annotation,allIons,model,0,0,.45);
    
    vector<float> ion_intensities1;
    vector<float> ion_intensities2;
    
    extractIons(psm1temp,peptideLength,model,ionsToExtract,ion_intensities1, 0, 0);
    extractIons(psm2temp,peptideLength,model,ionsToExtract,ion_intensities2, 0, 0);

    normalize_extracted_ions(ion_intensities1);
    normalize_extracted_ions(ion_intensities2);
    
    float cosine_val = cosine(ion_intensities1, ion_intensities2);
    
    psm1->m_spectrum = old_spectrum1;
    psm2->m_spectrum = old_spectrum2;
    
    return cosine_val;
}


float spectrum_similarity(psmPtr psm1, psmPtr psm2, int peptideLength, MS2ScoringModel &model,
        vector<string> &ionsToExtract, string allIons, vector<vector<float> > average_ions){
    
    psmPtr psm1temp(new PeptideSpectrumMatch);
    psmPtr psm2temp(new PeptideSpectrumMatch);
    
    Spectrum temp_spectrum1 = *(psm1->m_spectrum);
    Spectrum temp_spectrum2 = *(psm2->m_spectrum);
    
    Spectrum * old_spectrum1 = psm1->m_spectrum;
    Spectrum * old_spectrum2 = psm2->m_spectrum;
    
    psm1temp->m_spectrum = &temp_spectrum1;
    psm2temp->m_spectrum = &temp_spectrum2;
    psm1temp->m_annotation = psm1->m_annotation;
    psm2temp->m_annotation = psm2->m_annotation;
    
    //Making max 1000
    preprocess_spectrum_intensities_max_intensity(psm1temp->m_spectrum, 1000.f);
    preprocess_spectrum_intensities_max_intensity(psm2temp->m_spectrum, 1000.f);
    
    //Applying SQRT
    preprocess_spectrum_intensities(psm1temp->m_spectrum, 0, 1);
    preprocess_spectrum_intensities(psm2temp->m_spectrum, 0, 1);
    
    psm1temp->annotate(psm1temp->m_annotation,allIons,model,0,0,.45);
    psm2temp->annotate(psm2temp->m_annotation,allIons,model,0,0,.45);
    
    vector<float> ion_intensities1;
    vector<float> ion_intensities2;
    
    extractIons(psm1temp,peptideLength,model,ionsToExtract,ion_intensities1, 0, 0);
    extractIons(psm2temp,peptideLength,model,ionsToExtract,ion_intensities2, 0, 0);

    normalize_extracted_ions(ion_intensities1);
    normalize_extracted_ions(ion_intensities2);
    
    //Subtracting off average
    for(int i = 0; i < ion_intensities1.size(); i++){
        ion_intensities1[i] = ion_intensities1[i] - average_ions[peptideLength][i];
        ion_intensities2[i] = ion_intensities2[i] - average_ions[peptideLength][i];
    }
    
    normalize_extracted_ions(ion_intensities1);
    normalize_extracted_ions(ion_intensities2);
    
    float cosine_val = cosine(ion_intensities1, ion_intensities2);
    
    psm1->m_spectrum = old_spectrum1;
    psm2->m_spectrum = old_spectrum2;
    
    return cosine_val;
}


float full_spectrum_similarity(Spectrum spec1, Spectrum spec2){
    float max_mass = 0.f;
    
    Spectrum spec1_temp = spec1;
    Spectrum spec2_temp = spec2;
    spec1_temp.setResolution(1.0005, 1);
    spec2_temp.setResolution(1.0005, 1);

    preprocess_spectrum_intensities(&spec1_temp, 0,1);
    preprocess_spectrum_intensities(&spec2_temp, 0,1);
    
    //Finding the maximum mass in a spectrum
    for(int spec1_idx = 0; spec1_idx < spec1_temp.size(); spec1_idx++){
        if(max_mass < spec1_temp[spec1_idx][0])
            max_mass = spec1_temp[spec1_idx][0];
    }
    
    for(int spec2_idx = 0; spec2_idx < spec2_temp.size(); spec2_idx++){
        if(max_mass < spec2_temp[spec2_idx][0])
            max_mass = spec2_temp[spec2_idx][0];
    }
    
    
    vector<float> spec1_peak_vector( (int)max_mass + 1);
    vector<float> spec2_peak_vector( (int)max_mass + 1);
    for(int i = 0; i < spec1_temp.size(); i++)
        spec1_peak_vector[(int)spec1_temp[i][0]] = spec1_temp[i][1];
    for(int i = 0; i < spec2_temp.size(); i++)
        spec2_peak_vector[(int)spec2_temp[i][0]] = spec2_temp[i][1];
    
    normalize_extracted_ions(spec1_peak_vector);
    normalize_extracted_ions(spec2_peak_vector);
    
    
    float cosine_val = cosine(spec1_peak_vector, spec2_peak_vector);
    
    
    return cosine_val;
}


float full_spectrum_similarity(Spectrum spec1, Spectrum spec2, int & shared_peaks){
    float max_mass = 0.f;
    
    Spectrum spec1_temp = spec1;
    Spectrum spec2_temp = spec2;
    spec1_temp.setResolution(1.0005, 1);
    spec2_temp.setResolution(1.0005, 1);

    preprocess_spectrum_intensities(&spec1_temp, 0,1);
    preprocess_spectrum_intensities(&spec2_temp, 0,1);
    
    //Finding the maximum mass in a spectrum
    for(int spec1_idx = 0; spec1_idx < spec1_temp.size(); spec1_idx++){
        if(max_mass < spec1_temp[spec1_idx][0])
            max_mass = spec1_temp[spec1_idx][0];
    }
    
    for(int spec2_idx = 0; spec2_idx < spec2_temp.size(); spec2_idx++){
        if(max_mass < spec2_temp[spec2_idx][0])
            max_mass = spec2_temp[spec2_idx][0];
    }
    
    
    vector<float> spec1_peak_vector( (int)max_mass + 1);
    vector<float> spec2_peak_vector( (int)max_mass + 1);
    for(int i = 0; i < spec1_temp.size(); i++)
        spec1_peak_vector[(int)spec1_temp[i][0]] = spec1_temp[i][1];
    for(int i = 0; i < spec2_temp.size(); i++)
        spec2_peak_vector[(int)spec2_temp[i][0]] = spec2_temp[i][1];
    
    normalize_extracted_ions(spec1_peak_vector);
    normalize_extracted_ions(spec2_peak_vector);
    
    
    float cosine_val = cosine(spec1_peak_vector, spec2_peak_vector);
    shared_peaks = shared_peak_count(spec1_peak_vector, spec2_peak_vector);
    
    return cosine_val;
}


// Computes to cosine between two normalized vectors of the same size
int shared_peak_count(vector<float> &u, vector<float> &v) {
    if(u.size()!=v.size()) return 0;
    int count = 0;
    for(unsigned int i=0; i<u.size(); i++){
        if(u[i] > 0 && v[i] > 0){
            count++;
        }
    }
    return count;
}


//Calculating amino acid mass
float getMass(string peptide, map<char, float> amino_acid_map){
    float total_mass = 0.f;
    for(int i = 0; i < peptide.length(); i++){
        char amino_acid = tolower(peptide[i]);
        map<char,float>::iterator it;
        it = amino_acid_map.find(amino_acid);
        if(it != amino_acid_map.end()){
            total_mass += amino_acid_map[amino_acid];
        }
    }
    return total_mass;
}

float getSubstitutionCosineDepression(char a, char b, map<string, float> amino_acid_transformation_map){
	a = toupper(a);
	b = toupper(b);
	
	string lookup1 = "";
	lookup1 += a;
	lookup1 += b;
	string lookup2 = "";
	lookup2 += b;
	lookup2 += a;
	
	float lookup1_score = 0.f;
	float lookup2_score = 0.f;
	
	if(amino_acid_transformation_map.find(lookup1) != amino_acid_transformation_map.end()){
		lookup1_score = amino_acid_transformation_map[lookup1];
	}
	
	if(amino_acid_transformation_map.find(lookup2) != amino_acid_transformation_map.end()){
		lookup2_score = amino_acid_transformation_map[lookup2];
	}
	
	return max(lookup2_score, lookup1_score);
}

//Generate lookup table for amino acid substitutions and cosine depression
map<string, float> getAminoAcidSubstitutionLookup(){
    map<string, float> amino_acid_transform_map;
	amino_acid_transform_map["WT"] = 0.912276;
	amino_acid_transform_map["AS"] = 0.908659;
	amino_acid_transform_map["YA"] = 0.900532;
	amino_acid_transform_map["DI"] = 0.920199;
	amino_acid_transform_map["DL"] = 0.919106;
	amino_acid_transform_map["FE"] = 0.936028;
	amino_acid_transform_map["EV"] = 0.922741;
	amino_acid_transform_map["SG"] = 0.886207;
	amino_acid_transform_map["LN"] = 0.90118;
	amino_acid_transform_map["QV"] = 0.91915;
	return amino_acid_transform_map;
}

//Generate lookup table for amino acid masses
map<char, float> getAminoAcidLookup(){
    map<char, float> amino_acid_map;
    amino_acid_map['g'] = 57.02146;
    amino_acid_map['a'] = 71.03711f;
    amino_acid_map['s'] = 87.03203f;
    amino_acid_map['p'] = 97.05276f;
    amino_acid_map['v'] = 99.06841f;
    amino_acid_map['t'] = 101.04768f;
    amino_acid_map['c'] = 103.00919f;
    amino_acid_map['l'] = 113.08406f;
    amino_acid_map['i'] = 113.08406f;
    amino_acid_map['n'] = 114.04293f;
    amino_acid_map['d'] = 115.02694f;
    amino_acid_map['q'] = 128.05858f;
    amino_acid_map['k'] = 128.09496f;
    amino_acid_map['e'] = 129.04259f;
    amino_acid_map['m'] = 131.04049f;
    amino_acid_map['h'] = 137.05891f;
    amino_acid_map['f'] = 147.06841f;
    amino_acid_map['r'] = 156.10111f;
    amino_acid_map['y'] = 163.06333f;
    amino_acid_map['w'] = 186.07931f;
    return amino_acid_map;
}

string stripString(string peptide){
	string stripped_peptide = peptide;
	//Removing * and . at ends
	if(stripped_peptide[0] == '*'){
		stripped_peptide.erase(stripped_peptide.begin() + 0);
	}
	
	if(stripped_peptide[stripped_peptide.length()-1] == '*'){
		stripped_peptide.erase(stripped_peptide.begin() + stripped_peptide.length()-1);
	}
	
	if(stripped_peptide[0] == '.'){
		stripped_peptide.erase(stripped_peptide.begin() + 0);
	}
	
	if(stripped_peptide[stripped_peptide.length()-1] == '.'){
		stripped_peptide.erase(stripped_peptide.begin() + stripped_peptide.length()-1);
	}
	
	return stripped_peptide;
}

string getAdditionalMassString(int parent_mass_difference){
	//inserting the additional mass
	stringstream mass_addition_stream;
	if(parent_mass_difference > 0){
		mass_addition_stream<<"["<<parent_mass_difference<<"]";
	}
	else if(parent_mass_difference < 0){
		mass_addition_stream<<"["<<parent_mass_difference<<"]";
	}
	//otherwise equal, so do nothing
	
	return mass_addition_stream.str();
}


string cleanAnnotationEnds(string annotation){
    if(annotation[0] != '*'){
            annotation[0] = '*';
            //cout<<new_peptides[i]<<"\t"<<annotation<<endl;
    }
    if(annotation[annotation.length() - 1] != '*'){
        annotation[annotation.length()-1] = '*';
        //cout<<new_peptides[i]<<"\t"<<annotation<<endl;
    }
    return annotation;
}


string create_annotation_ends(string annotation){
    int firstdotloc = annotation.find_first_of(".");
    int lastdotloc = annotation.find_last_of(".");
    if(lastdotloc != (annotation.length() - 2)) annotation.insert(annotation.length(), ".*");
    if(firstdotloc != 1) annotation.insert(0, "*.");
    
    return annotation;    
}

string remove_annotation_ends(string annotation){
    int firstdotloc = annotation.find_first_of(".");
    int lastdotloc = annotation.find_last_of(".");
    if(lastdotloc == (annotation.length() - 2)) annotation = annotation.substr(0, annotation.length() - 2);
    if(firstdotloc == 1) annotation = annotation.substr(2, annotation.length() - 2);
    
    return annotation;
}



int getDifferenceIndex(string in1, string in2){
    if(in1.length() != in2.length())
        return -1;
    int difference = 0;
    for(int i = 0; i < in1.length(); i++){
        if(in1[i] != in2[i]){
            return i;
        }
    }
    return -1;
}

int getStringDifference(string in1, string in2){
    if(in1.length() != in2.length())
        return -1;
    int difference = 0;
    for(int i = 0; i < in1.length(); i++){
        if(in1[i] != in2[i]){
            difference++;
        }
    }
    return difference;
}

// Computes to cosine between two normalized vectors of the same size
float cosine(vector<float> &u, vector<float> &v) {
    if(u.size()!=v.size()) return 0;
    float c=0;
    for(unsigned int i=0; i<u.size(); i++) c+=u[i]*v[i];
    return c;
}

void preprocess_spectrum_intensities(Spectrum * s, int do_early_normalize, int do_filter){
    if(do_filter == 1){
        for(int peakIdx = 0; peakIdx < s->size(); peakIdx++){
            (*s)[peakIdx][1] =  sqrt((*s)[peakIdx][1]);
        }
    }
    if(do_filter == 2){
        for(int peakIdx = 0; peakIdx < s->size(); peakIdx++){
            (*s)[peakIdx][1] =  log((*s)[peakIdx][1]);
        }
    }
    
    //normalization
    //cout<<"preprocesses"<<endl;
    float total_intensity = 0.f;
    for(int peakIdx = 0; peakIdx < s->size(); peakIdx++){
        total_intensity += (*s)[peakIdx][1] * (*s)[peakIdx][1];
    }
    total_intensity = sqrt(total_intensity);
    
    for(int peakIdx = 0; peakIdx < s->size(); peakIdx++){
        (*s)[peakIdx][1] = (*s)[peakIdx][1]/total_intensity;
    }
}


void preprocess_spectrum_intensities_max_intensity(Spectrum * s, float max){
    float current_max_intesity = 0.f;
    for(int i = 0; i < s->size(); i++){
        if(current_max_intesity < (*s)[i][1]) current_max_intesity = (*s)[i][1];
    }
    
    for(int i = 0; i < s->size(); i++){
        (*s)[i][1] = (*s)[i][1] * max / current_max_intesity;
    }
}

void normalize_extracted_ions(vector<float> &v) {
    //if(do_early_normalize == 1)
    //    return;
    
    float n=0;
    for(unsigned int i=0; i<v.size(); i++) n+=v[i]*v[i];
        if(n>0) {
        n=sqrt(n);
        for(unsigned int i=0; i<v.size(); i++) v[i]/=n;
    }
}

void extractIons(psmPtr psm, int peptideLength, MS2ScoringModel &model,
        vector<string> &ionsToExtract, vector<float> &ions, int do_early_normalize, int do_filter) {
	ions.resize(peptideLength*ionsToExtract.size());
	for(unsigned int i=0; i<ions.size(); i++) ions[i]=0;
	
	if(psm->m_peakAnnotations.size()!=psm->m_spectrum->size()) return;
	
	//Normalize before extraction, so that we take into account noise
    //normalize_spectrum(s, do_early_normalize, do_filter);

	// Select target ions
	unsigned int baseIdx=0;
	for(unsigned int ionType=0; ionType<ionsToExtract.size(); ionType++) {
        for(unsigned int peakIdx=0; peakIdx<psm->m_peakAnnotations.size(); peakIdx++)
            if(psm->m_peakAnnotations[peakIdx].first and psm->m_peakAnnotations[peakIdx].first->name==ionsToExtract[ionType]){
				ions[baseIdx+psm->m_peakAnnotations[peakIdx].second-1]=(*psm->m_spectrum)[peakIdx][1];
			}
		baseIdx+=peptideLength;
	}
}

void extractIons(psmPtr psm, int peptideLength, MS2ScoringModel &model,
        vector<string> &ionsToExtract, vector<pair <float, float> > &ions, int do_early_normalize, int do_filter) {
	ions.resize(peptideLength*ionsToExtract.size());
	for(unsigned int i=0; i<ions.size(); i++){
		ions[i].first=0;
		ions[i].second=0;
	}
	if(psm->m_peakAnnotations.size()!=psm->m_spectrum->size()) return;

	//Normalize before extraction, so that we take into account noise
    //normalize_spectrum(s, do_early_normalize, do_filter);

	// Select target ions
	unsigned int baseIdx=0;
	for(unsigned int ionType=0; ionType<ionsToExtract.size(); ionType++) {
		for(unsigned int peakIdx=0; peakIdx<psm->m_peakAnnotations.size(); peakIdx++)
			if(psm->m_peakAnnotations[peakIdx].first and psm->m_peakAnnotations[peakIdx].first->name==ionsToExtract[ionType]){
				ions[baseIdx+psm->m_peakAnnotations[peakIdx].second-1].second=(*psm->m_spectrum)[peakIdx][1];
				ions[baseIdx+psm->m_peakAnnotations[peakIdx].second-1].first=(*psm->m_spectrum)[peakIdx][0];
			}
		baseIdx+=peptideLength;
	}
}


void norm_vector(vector<float> & input){
    float sum = 0.f;
    for(int i = 0; i < input.size(); i++){
        sum += input[i]  * input[i];
    }
    
    if(sum == 0.f) return;
    
    float euclidian_norm = sqrt(sum);
    for(int i = 0; i < input.size(); i++){
        input[i] = input[i] / euclidian_norm;
    }
}

//Comparater function for sorting
bool search_results_comparator (score_results_tuple i, score_results_tuple j){
    return tr1::get<1>(i) > tr1::get<1>(j);
}

bool search_results_comparator_psmPtr (psmPtr i, psmPtr j){
    return i->m_score > j->m_score;
}


string reverse_string(string input){
    string output = "";
    for(int i = 0; i < input.length(); i++){
        output += input[input.length() - 1 - i];
    }
    return output;
}

float full_spectrum_dotbias(Spectrum spec1, Spectrum spec2, float spec_sim){
    float max_mass = 0.f;
    
    Spectrum spec1_temp = spec1;
    Spectrum spec2_temp = spec2;
    spec1_temp.setResolution(1.0005, 1);
    spec2_temp.setResolution(1.0005, 1);

    preprocess_spectrum_intensities(&spec1_temp, 0,1);
    preprocess_spectrum_intensities(&spec2_temp, 0,1);
    
    //Finding the maximum mass in a spectrum
    for(int spec1_idx = 0; spec1_idx < spec1_temp.size(); spec1_idx++){
        if(max_mass < spec1_temp[spec1_idx][0])
            max_mass = spec1_temp[spec1_idx][0];
    }
    
    for(int spec2_idx = 0; spec2_idx < spec2_temp.size(); spec2_idx++){
        if(max_mass < spec2_temp[spec2_idx][0])
            max_mass = spec2_temp[spec2_idx][0];
    }
    
    for(int i = 0; i < spec1_temp.size(); i++){
        //cout<<spec1_temp.peakList[i][0]<<"\t"<<spec1_temp.peakList[i][1]<<"\t"<<spec1.peakList[i][0]<<"\t"<<spec1.peakList[i][1]<<endl;
    }
    
    vector<float> spec1_peak_vector( (int)max_mass + 1);
    vector<float> spec2_peak_vector( (int)max_mass + 1);
    for(int i = 0; i < spec1_temp.size(); i++)
        spec1_peak_vector[(int)spec1_temp[i][0]] = spec1_temp[i][1]*spec1_temp[i][1];
    for(int i = 0; i < spec2_temp.size(); i++)
        spec2_peak_vector[(int)spec2_temp[i][0]] = spec2_temp[i][1]*spec2_temp[i][1];
    
    normalize_extracted_ions(spec1_peak_vector);
    normalize_extracted_ions(spec2_peak_vector);
    
    for(int i = 0; i < spec1_peak_vector.size(); i++)
        spec1_peak_vector[i] = spec1_peak_vector[i] * spec1_peak_vector[i];
    for(int i = 0; i < spec2_peak_vector.size(); i++)
        spec2_peak_vector[i] = spec2_peak_vector[i] * spec2_peak_vector[i];
    
    
    
    
    float cosine_val = cosine(spec1_peak_vector, spec2_peak_vector);
    
    float sqrt_cosine_val = sqrt(cosine_val);
    
    float DB = sqrt_cosine_val/spec_sim;
    
    //return DB;
    
    if(DB < 0.1)
        return 0.12;
    if(DB <= 0.4 && DB > 0.35)
        return .12;
    if(DB > 0.4 && DB <= 0.45)
        return 0.18;
    if(DB > 0.45)
        return 0.24;
    
    return 0;
}

void sorted_vector_library(vector<Spectrum *> &library_ptr, SpectralLibrary & real_library){
    for(int i = 0; i < real_library.size(); i++){
        library_ptr.push_back(&real_library[i]);
    }
    
    sort(library_ptr.begin(), library_ptr.end(), spectrum_ptr_compare);
    
}

bool spectrum_ptr_compare (Spectrum* first, Spectrum* second){ 
    return (first->parentMZ<second->parentMZ);
}

int spectrum_ptr_startend(vector<Spectrum *> &library_ptr, float parentMZ, float tolerance, int &start_idx, int &end_idx){
    if(library_ptr.size() == 0){
        start_idx = 0;
        end_idx = -1;
        return 0;
    }
    
    float looking_mass = parentMZ - tolerance;
    
    int start_range = 0;
    int end_range = library_ptr.size()-1;
    while(end_range - start_range > 1){
        int middle = ( start_range + end_range )/2;
        if(looking_mass > library_ptr[middle]->parentMZ){
            start_range = middle + 1;
            continue;
        }
        if(looking_mass < library_ptr[middle]->parentMZ){
            end_range = middle - 1;
            continue;
        }
        if(looking_mass == library_ptr[middle]->parentMZ){
            start_range = middle;
            end_range = middle;
            break;
        }
    }
    
    start_idx = start_range;
    
    looking_mass= parentMZ + tolerance;
    start_range = 0;
    end_range = library_ptr.size()-1;
    while(end_range - start_range > 1){
        int middle = ( start_range + end_range )/2;
        if(looking_mass > library_ptr[middle]->parentMZ){
            start_range = middle + 1;
            continue;
        }
        if(looking_mass < library_ptr[middle]->parentMZ){
            end_range = middle - 1;
            continue;
        }
        if(looking_mass == library_ptr[middle]->parentMZ){
            start_range = middle;
            end_range = middle;
            break;
        }
    }
    
    end_idx = end_range;
    
    return 0;
}

int getpeptideLength(string annotation){
    int length = 0;
    for(int i = 0; i < annotation.length(); i++){
        if(annotation[i] == 'A' || 
            annotation[i] == 'R' ||
            annotation[i] == 'N' ||
            annotation[i] == 'D' ||
            annotation[i] == 'C' ||
            annotation[i] == 'E' ||
            annotation[i] == 'Q' ||
            annotation[i] == 'G' ||
            annotation[i] == 'H' ||
            annotation[i] == 'I' ||
            annotation[i] == 'L' ||
            annotation[i] == 'K' ||
            annotation[i] == 'M' ||
            annotation[i] == 'F' ||
            annotation[i] == 'P' ||
            annotation[i] == 'S' ||
            annotation[i] == 'T' ||
            annotation[i] == 'W' ||
            annotation[i] == 'Y' ||
            annotation[i] == 'V' ||
            annotation[i] == '[')
            length++;    
    }
    return length;
}

int hashpeptide(string input, int charge){
    input += 'a' + charge;
    return hashstring (input);
}

int hashstring(string input){
    int hash = 0;
    int offset = 'A' - 1;
    for(string::const_iterator it=input.begin(); it!=input.end(); ++it) {
        hash = hash + (*it);
    }
    return hash;
    
}

string get_only_filename(string input_path){
    int last_slash = input_path.find_last_of("/");
    if(last_slash != string::npos){
        return input_path.substr(last_slash+1);
    }
    return input_path;
}

string get_only_path(string input_path){
    int last_slash = input_path.find_last_of("/");
    if(last_slash != string::npos){
        return input_path.substr(0, last_slash);
    }
    return input_path;
}

string remove_char(string input, char removal){
    for(int i = 0; i < input.length(); i++){
        if(input[i] == removal){
            input.erase(input.begin() + i);
            i--;
        }
    }
    return input;
}

vector<string> create_deliminated_aminoacids(string peptide){
    string stripped_peptide = remove_annotation_ends(peptide);
    vector<string> deliminated_aminoacids;
    
    for(int pepidx = 0; pepidx < stripped_peptide.length(); pepidx++){
        if(stripped_peptide[pepidx] == '*' || stripped_peptide[pepidx] == '.'){
            continue;
        }
        
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
    
    return deliminated_aminoacids;
}

void sqrt_vector(vector<float> & input){
    for(int i = 0; i < input.size(); i++){
        input[i] = sqrt(input[i]);
    }
}

float spectrum_similarity_sqrt_librarypeaks_isocombine(psmPtr library,
                          psmPtr query,
                          int peptideLength, 
                          MS2ScoringModel &model,
                          vector<string> &ionsToExtract, 
                          string allIons, 
                          float &percent_explained_intensity){
    
    return spectrum_similarity_sqrt_librarypeaks_isocombine_speed_optimization(library, query, peptideLength, model, ionsToExtract, allIons, percent_explained_intensity);
}

float spectrum_similarity_sqrt_librarypeaks_isocombine_speed_optimization(psmPtr library,
                          psmPtr query,
                          int peptideLength, 
                          MS2ScoringModel &model,
                          vector<string> &ionsToExtract, 
                          string allIons, 
                          float &percent_explained_intensity){
    
    //return 0.f;
    
    
    if(library->m_ion_extraction.size() == 0){
        preprocess_library_ion_extraction_isocombine(library, peptideLength, model, ionsToExtract, allIons);
    }
    
    psmPtr psm2temp(new PeptideSpectrumMatch);
    
    Spectrum temp_spectrum2 = *(query->m_spectrum);
    
    Spectrum * old_spectrum2 = query->m_spectrum;
    
    psm2temp->m_spectrum = &temp_spectrum2;
    psm2temp->m_annotation = query->m_annotation;
    
    
    
    float explained_intensity = preprocess_library_ion_extraction_isocombine(psm2temp, peptideLength, model, ionsToExtract, allIons);
    
    
    for(int i = 0; i < library->m_ion_extraction.size(); i++){
        if(library->m_ion_extraction[i] < 0.0001){
            psm2temp->m_ion_extraction[i] = 0.f;
        }
    }
    
    
    normalize_extracted_ions(psm2temp->m_ion_extraction);
    
    float cosine_val = cosine(library->m_ion_extraction, psm2temp->m_ion_extraction);
    
    query->m_spectrum = old_spectrum2;
    
   
    return cosine_val;
}

float spectrum_similarity_sqrt_librarypeaks(psmPtr library,
                          psmPtr query,
                          int peptideLength, 
                          MS2ScoringModel &model,
                          vector<string> &ionsToExtract, 
                          string allIons,
                          float &percent_explained_intensity){
    
    return spectrum_similarity_sqrt_librarypeaks_speed_optimization(library, query, peptideLength, model, ionsToExtract, allIons, percent_explained_intensity);
}

float spectrum_similarity_sqrt_librarypeaks_speed_optimization(psmPtr library,
                          psmPtr query,
                          int peptideLength, 
                          MS2ScoringModel &model,
                          vector<string> &ionsToExtract, 
                          string allIons, 
                          float &percent_explained_intensity){
    

    if(library->m_ion_extraction.size() == 0){
        DEBUG_MSG("EXTRACTING LIBRARY");
        DEBUG_MSG(library->m_annotation<<"\t"<<library->m_ion_extraction.size());
        //Removing Precursor Peaks
        /*std::vector<std::string> ionsToExtract_precursor;
        ionsToExtract_precursor.push_back("P++");
        ionsToExtract_precursor.push_back("P++-H2O");
        ionsToExtract_precursor.push_back("P++-NH3");
        ionsToExtract_precursor.push_back("P++-H2O-H2O");
        
        library->annotate(library->m_annotation,allIons,model,0,0,.45);
        deleteIons(library,peptideLength,model,ionsToExtract_precursor);
        
        preprocess_spectrum_intensities_max_intensity(library->m_spectrum, 1000.f);
        
        preprocess_spectrum_intensities(library->m_spectrum, 0, 1);
        
        library->annotate(library->m_annotation,allIons,model,0,0,.45);
        extractIons(library,peptideLength,model,ionsToExtract,library->m_ion_extraction, 0, 0);
        normalize_extracted_ions(library->m_ion_extraction);*/
        preprocess_library_ion_extraction(library, peptideLength, model, ionsToExtract, allIons);
    }
    
    
    psmPtr psm2temp(new PeptideSpectrumMatch);
    
    Spectrum temp_spectrum2 = *(query->m_spectrum);
    
    Spectrum * old_spectrum2 = query->m_spectrum;
    
    psm2temp->m_spectrum = &temp_spectrum2;
    psm2temp->m_annotation = query->m_annotation;
    
    /*
    //Removing Precursor Peaks
    std::vector<std::string> ionsToExtract_precursor;
    ionsToExtract_precursor.push_back("P++");
    ionsToExtract_precursor.push_back("P++-H2O");
    ionsToExtract_precursor.push_back("P++-NH3");
    ionsToExtract_precursor.push_back("P++-H2O-H2O");
    
    psm2temp->annotate(psm2temp->m_annotation,allIons,model,0,0,.45);
    
    deleteIons(psm2temp,peptideLength,model,ionsToExtract_precursor);
    
    
    //Making max 1000
    preprocess_spectrum_intensities_max_intensity(psm2temp->m_spectrum, 1000.f);
    
    //Applying SQRT and euclidian norm
    preprocess_spectrum_intensities(psm2temp->m_spectrum, 0, 1);
    
    psm2temp->annotate(psm2temp->m_annotation,allIons,model,0,0,.45);
    
    vector<float> ion_intensities2;
    
    
    
    extractIons(psm2temp,peptideLength,model,ionsToExtract,ion_intensities2, 0, 0);*/
    
    


    float explained_intensity = preprocess_library_ion_extraction(psm2temp, peptideLength, model, ionsToExtract, allIons);

    percent_explained_intensity = explained_intensity;
    //combine_isotopic(ionsToExtract, ion_intensities1, peptideLength);
    //combine_isotopic(ionsToExtract, ion_intensities2, peptideLength);

    for(int i = 0; i < library->m_ion_extraction.size(); i++){
        if(library->m_ion_extraction[i] < 0.0001){
            psm2temp->m_ion_extraction[i] = 0.f;
        }
    }
    
    //DEBUG_MSG("NEW EXPLAINED INTENSITY\t"<<explained_intensity);

    normalize_extracted_ions(psm2temp->m_ion_extraction);
    
    
    
    float cosine_val = cosine(library->m_ion_extraction, psm2temp->m_ion_extraction);
    
    //if(library->m_annotation == "*.AIIVLSTSGTTPR.*" && query->m_spectrum->scan == 8322)
    //if(library->m_annotation == "*.ANGTTVLVGMPAGAK.*")
    //if(library->m_annotation == "*.ALNEEAEAR.*")
    //if(library->m_annotation == "*.VTPSFVAFTPEER.*")
    if(false)
    if(cosine_val*explained_intensity > 0.4){
        #pragma omp critical
        {
            query->annotate(library->m_annotation,allIons,model,0,0,.50);
            
            DEBUG_MSG("SCAN: "<<query->m_spectrum->scan<<"\t"<<library->m_annotation);
            //printDoubleIon(ionsToExtract, library->m_ion_extraction, psm2temp->m_ion_extraction);
            DEBUG_MSG("ANNOTATION SIZE\t"<<library->m_peakAnnotations.size());
            
            for(int i = 0; i < query->m_peakAnnotations.size(); i++){
                cout<<"QUERY\t"<<(*query->m_spectrum)[i][0] << "\t"<< (*query->m_spectrum)[i][1] << "\t";
                
                stringstream ss (stringstream::in | stringstream::out);
                //ss<<(*library->m_spectrum)[i][0] << "\t"<< (*library->m_spectrum)[i][1] << "\t";
                if(query->m_peakAnnotations[i].first){
                    cout<<query->m_peakAnnotations[i].first->name<<"\t"<<query->m_peakAnnotations[i].second;
                    //ss<<library->m_peakAnnotations[i].first->name;
                }
                cout<<endl;
                //DEBUG_MSG(ss.str());
            }
            
            for(int i = 0; i < library->m_peakAnnotations.size(); i++){
                cout<<"LIBRARY\t"<<(*library->m_spectrum)[i][0] << "\t"<< (*library->m_spectrum)[i][1] << "\t";
                
                //ss<<(*library->m_spectrum)[i][0] << "\t"<< (*library->m_spectrum)[i][1] << "\t";
                if(library->m_peakAnnotations[i].first){
                    cout<<library->m_peakAnnotations[i].first->name<<library->m_peakAnnotations[i].second;
                    //ss<<library->m_peakAnnotations[i].first->name;
                }
                cout<<endl;
                //DEBUG_MSG(ss.str());
            }
        }
    }
    
    query->m_spectrum = old_spectrum2;
    
    /*
    for(int i = 0; i < ion_intensities2.size(); i++){
        DEBUG_MSG(ion_intensities2[i]<<"\t"<<library->m_ion_extraction[i]);
    }
    
    DEBUG_MSG("COSINE SPEED\t"<<cosine_val);
    */
    return cosine_val;
}

//Returns percent explained intensity
float preprocess_library_ion_extraction(psmPtr library,
                          int peptideLength, 
                          MS2ScoringModel &model,
                          vector<string> &ionsToExtract, 
                          string allIons){

    //DEBUG_MSG("EXTRACTING LIBRARY PREPROCESS");
    //Removing Precursor Peaks
    std::vector<std::string> ionsToExtract_precursor;
    ionsToExtract_precursor.push_back("P++");
    ionsToExtract_precursor.push_back("P++-H2O");
    ionsToExtract_precursor.push_back("P++-NH3");
    ionsToExtract_precursor.push_back("P++-H2O-H2O");
    
    //library->annotate(library->m_annotation,allIons,model,0,0,.45);
    //deleteIons(library,peptideLength,model,ionsToExtract_precursor);
    
    preprocess_spectrum_intensities_max_intensity(library->m_spectrum, 1000.f);
    
    preprocess_spectrum_intensities(library->m_spectrum, 0, 1);
    
    library->annotate(library->m_annotation,allIons,model,0,0,.50);
    extractIons(library,peptideLength,model,ionsToExtract,library->m_ion_extraction, 0, 0);
    
    float percent_explained_intensity = 0.f;
    
    for(int i = 0; i < library->m_ion_extraction.size(); i++){
        //if(library->SLGF_distribution.size() == 0)
        //    DEBUG_MSG(library->m_ion_extraction[i]);
        percent_explained_intensity += library->m_ion_extraction[i] * library->m_ion_extraction[i];
    }
    
    percent_explained_intensity = sqrt(percent_explained_intensity);
    
    normalize_extracted_ions(library->m_ion_extraction);
    
    return percent_explained_intensity;


}

//Returns percent explained intensity
float preprocess_library_ion_extraction_isocombine(psmPtr library,
                          int peptideLength, 
                          MS2ScoringModel &model,
                          vector<string> &ionsToExtract, 
                          string allIons){

    //DEBUG_MSG("EXTRACTING LIBRARY PREPROCESS");
    //Removing Precursor Peaks
    std::vector<std::string> ionsToExtract_precursor;
    ionsToExtract_precursor.push_back("P++");
    ionsToExtract_precursor.push_back("P++-H2O");
    ionsToExtract_precursor.push_back("P++-NH3");
    ionsToExtract_precursor.push_back("P++-H2O-H2O");
    
    //library->annotate(library->m_annotation,allIons,model,0,0,.45);
    //deleteIons(library,peptideLength,model,ionsToExtract_precursor);
    
    
    preprocess_spectrum_intensities_max_intensity(library->m_spectrum, 1000.f);
    
    preprocess_spectrum_intensities(library->m_spectrum, 0, 1);
    
    library->annotate(library->m_annotation,allIons,model,0,0,.50);
    extractIons(library,peptideLength,model,ionsToExtract,library->m_ion_extraction, 0, 0);
    
    vector<float> temp_ions = library->m_ion_extraction;
    
    combine_isotopic(ionsToExtract, library->m_ion_extraction, peptideLength);
    
    float percent_explained_intensity = 0.f;
    
    for(int i = 0; i < library->m_ion_extraction.size(); i++){
        //percent_explained_intensity += library->m_ion_extraction[i] * library->m_ion_extraction[i];
        percent_explained_intensity += temp_ions[i] * temp_ions[i];
    }
    
    percent_explained_intensity = sqrt(percent_explained_intensity);
    
    /*
    if("*.YSSFTR.*" == library->m_annotation){
        DEBUG_MSG(peptideLength);
        DEBUG_MSG(library->m_annotation);
        DEBUG_MSG(library->m_ion_extraction.size());
        DEBUG_MSG(ionsToExtract.size());
        for(int i = 0; i < library->m_ion_extraction.size(); i++){
            DEBUG_MSG(i<<"\t"<<temp_ions[i]);
        }
    }*/
    
    //if(percent_explained_intensity > 1.0){
        //DEBUG_MSG(library->m_annotation<<"\t"<<percent_explained_intensity);
        //for(int i = 0; i < library->m_ion_extraction.size(); i++){
        //    DEBUG_MSG(temp_ions[i]<<"\t"<<library->m_ion_extraction[i]);
        //}
    //}
    
    normalize_extracted_ions(library->m_ion_extraction);
    
    return percent_explained_intensity;

}


float max_ion_intensity(vector<float> ion_intensity){
    float max = 0.f;
    for(int i = 0; i < ion_intensity.size(); i++){
        if ( max < ion_intensity[i] ){
            max = ion_intensity[i];
        }
    }
    return max;
}

vector<pair<float, float> > create_histogram(int buckets, float start_range, float end_range, vector<float> input, bool normalize){
    vector<pair<float, float> > histogram;
    for(int i = 0 ; i < buckets; i++){
        pair<float, float> histogram_cell;
        histogram_cell.first = start_range + (end_range - start_range)/buckets * i + (end_range - start_range)/buckets * 0.5;
        histogram_cell.second = 0;
        histogram.push_back(histogram_cell);
    }
    
    int inbounds_count = 0;
    for(int i = 0 ; i < input.size(); i++){
        int histogram_idx = (input[i]-start_range)/(end_range - start_range) * buckets;
        if(histogram_idx >= buckets || histogram_idx < 0) continue;
        
        histogram[histogram_idx].second++;
        inbounds_count++;
        //cout<<input[i]<<"\t"<<histogram[histogram_idx].first<<endl;
    }
    
    if(normalize){
        for(int i = 0 ; i < buckets; i++){
            histogram[i].second = histogram[i].second/inbounds_count;
        }
        if(inbounds_count != input.size()){
            cout<<input.size() - inbounds_count<<" Fell outside range"<<endl;
        }
    }
    
    return histogram;
    
}


void combine_isotopic(vector<string> ions, vector<float> &extracted_ions, int length){
    for(int i = 0; i < ions.size(); i++){
        for(int j = 0; j < ions.size(); j++){
            string first_ion = ions[i];
            string second_ion = ions[j];
            string first_ion_iso = first_ion + "-iso";
            if(first_ion_iso == second_ion){
                int first_idx = i*length;
                int second_idx = j*length;
                for(int index = 0; index < length; index++){
                    extracted_ions[first_idx+index] += extracted_ions[second_idx+index];
                    extracted_ions[second_idx+index] = 0.f;
                }
                //cout<<first_ion<<"\t"<<second_ion<<endl;
            }
        }
    }
}


void dynamicprogramming_spectrum(Spectrum &spectrum, MS2ScoringModel model, vector<string> ions_to_extract, string allIons, vector<float> ion_mean, vector<float> ion_variance, vector< pair < float, float > > insertion_probability , vector< pair < float, float > > deletion_probability, vector< pair< float, float> > &output_distribution, vector<vector<pair<float, float> > > histograms, int training_model){
    
    stringstream output_stream (stringstream::in | stringstream::out);;
    
    //Extracting Ions
    psmPtr psm = spectrum.psmList.front();
    psm->annotate(psm->m_annotation,allIons,model,0,0,.45);
    vector<float> ion_mass_intensity_pair_vector;
    
    //Removing Precursor Peaks
    std::vector<std::string> ionsToExtract_precursor;
    ionsToExtract_precursor.push_back("P++");
    ionsToExtract_precursor.push_back("P++-H2O");
    ionsToExtract_precursor.push_back("P++-NH3");
    ionsToExtract_precursor.push_back("P++-H2O-H2O");
    
    deleteIons(psm, create_deliminated_aminoacids(psm->m_annotation).size(), model,ionsToExtract_precursor);

    extractIons(psm,create_deliminated_aminoacids(psm->m_annotation).size(),model,ions_to_extract,ion_mass_intensity_pair_vector, 0, 0);
    
    
    sqrt_vector(ion_mass_intensity_pair_vector);
    normalize_extracted_ions(ion_mass_intensity_pair_vector);
    
    
    dynamicprogramming_ionvector(ion_mass_intensity_pair_vector,
                                  ion_mean, 
                                  ion_variance, 
                                  insertion_probability, 
                                  deletion_probability, 
                                  output_distribution, 
                                  histograms, 
                                  training_model,
                                  800,
                                  400);
}


void dynamicprogramming_ionvector(vector<float> ion_intensities,
                                  vector<float> ion_mean, 
                                  vector<float> ion_variance, 
                                  vector< pair < float, float > > insertion_probability , 
                                  vector< pair < float, float > > deletion_probability, 
                                  vector< pair< float, float> > &output_distribution, 
                                  vector<vector<pair<float, float> > > histograms, 
                                  int training_model,
                                  int cosine_dimension,
                                  int intensity_dimension){
    
    stringstream output_stream (stringstream::in | stringstream::out);;
    

    vector<float> ion_mass_intensity_pair_vector;
    //sqrt_vector(ion_mass_intensity_pair_vector);
    normalize_extracted_ions(ion_intensities);
    
    //Sorting in descedning ions
    ion_mass_intensity_pair_vector = sort_descending(ion_intensities);
    
    for(int i = 0; i < ion_mass_intensity_pair_vector.size(); i++){
        if(ion_mass_intensity_pair_vector[i] < 0.00001){
            ion_mass_intensity_pair_vector.erase(ion_mass_intensity_pair_vector.begin()+i);
            i--;
        }
    }
    
    ion_mass_intensity_pair_vector.insert(ion_mass_intensity_pair_vector.begin(), 1, 0.f);
    
    output_stream<<"NONZERO Size:\t"<<ion_mass_intensity_pair_vector.size()<<endl;

    int intensity_bins = intensity_dimension;
    int cosine_bins = cosine_dimension;
    int ion_size = ion_mass_intensity_pair_vector.size();
    
    int dynamic_prog_table_size = intensity_bins * cosine_bins * ion_size;
    
    output_stream<<"Memory footprint\t"<<dynamic_prog_table_size*sizeof(double)<<endl;
    
    double * programming_matrix = new double[dynamic_prog_table_size];
    
    for(int intensity_idx = 0; intensity_idx < intensity_bins; intensity_idx++){
        for(int cosine_idx = 0; cosine_idx< cosine_bins; cosine_idx++){
            for(int ion_idx = 0; ion_idx < ion_size; ion_idx++){
                set_3d_array(intensity_bins, cosine_bins, ion_size, intensity_idx, cosine_idx, ion_idx, programming_matrix, 0.f);
            }
        }
    }
    
    //Initializing table
    set_3d_array(intensity_bins, cosine_bins, ion_size, 0, 0, 0, programming_matrix, 1.f);
    
    
    double max_intensity = max_ion_intensity(ion_mass_intensity_pair_vector);
    
    
    //======================================================================
    //Precomputation, we're really tied to using squared intensity as the index
    vector<double> intensity_difference_to_expected_ion;
    for(int i = 0; i < intensity_bins; i++){
        double expected_ion = sqrt(((double)(i))/((double)intensity_bins));
        intensity_difference_to_expected_ion.push_back(expected_ion);
    }
    
   //=======================================================================
    
    for(int ion_idx = 0; ion_idx < ion_size-1; ion_idx++){
        double ion_intensity = ion_mass_intensity_pair_vector[ion_idx+1];
        float percent_intensity = min(ion_intensity/max_intensity, 0.999999);
        int band_idx = min( ((int)(min(ion_intensity/max_intensity, 0.999999)/(1.0/ion_variance.size()))), ((int)(ion_variance.size() - 1)));
        DEBUG_MSG("Progress\t"<<ion_idx<<" of "<<ion_size);

        
        for(int intensity_idx = 0; intensity_idx < intensity_bins; intensity_idx++){
            for(int cosine_idx = 0; cosine_idx< cosine_bins; cosine_idx++){
                
                double expected_cosine = ((float)cosine_idx)/((float)cosine_bins);
                //double intensity_used = sqrt(((float)intensity_idx)/((float)intensity_bins));
                //double intensity_used = (((float)intensity_idx)/((float)intensity_bins));
                
                double current_probability  = get_3d_array(intensity_bins, cosine_bins, ion_size, intensity_idx, cosine_idx, ion_idx, programming_matrix);
                
                if(current_probability < 0.0000000001) continue;
                
                
                //Deletion Probability
                float specific_deletion_probability = 0.f;
                for(int deletion_idx = 0; deletion_idx < deletion_probability.size(); deletion_idx++){
                    if(percent_intensity > deletion_probability[deletion_idx].first){
                        specific_deletion_probability = deletion_probability[deletion_idx].second;
                    }
                    else{
                        break;
                    }
                }
                
                double next_same_cell_probability = get_3d_array(intensity_bins, cosine_bins, ion_size, intensity_idx, cosine_idx, ion_idx+1, programming_matrix);
                
                double new_probability_for_deletion = current_probability * specific_deletion_probability + next_same_cell_probability;
                
                set_3d_array(intensity_bins, cosine_bins, ion_size, intensity_idx, cosine_idx, ion_idx+1, programming_matrix, new_probability_for_deletion);
                
                vector<double> points_in_gaussain;
                //Now we start looking back at each column 
                for(int forward_intensity_idx = intensity_idx; forward_intensity_idx < intensity_bins; forward_intensity_idx++){
                    //double forward_intensity_used = (((float)forward_intensity_idx+1)/((float)intensity_bins));
                    //double forward_intensity_used = (((float)forward_intensity_idx)/((float)intensity_bins));
                    //double forward_expection_ion = sqrt(forward_intensity_used*forward_intensity_used - intensity_used*intensity_used);
                    
                    //double forward_expection_ion = sqrt(forward_intensity_used - intensity_used);
                    double forward_expection_ion = intensity_difference_to_expected_ion[forward_intensity_idx-intensity_idx];
                    
                    double log_ratio = log(forward_expection_ion/ion_intensity)/log(2);
                    
                    double forward_cosine_value = max(expected_cosine + forward_expection_ion * ion_intensity, 0.0);
                    
                    //cout<<ion_idx<<"\t"<<intensity_idx<<"\t"<<cosine_idx<<"\t"<<forward_intensity_idx<<"\t"<<forward_intensity_used<<"\t"<<forward_expection_ion<<"\t"<<log_ratio<<"\t"<<forward_cosine_value<<endl;
                    
                    points_in_gaussain.push_back(log_ratio);
                }
                
                double max_gaussian_lookup = points_in_gaussain[points_in_gaussain.size()-1];
                
                for(int forward_intensity_idx = intensity_idx; forward_intensity_idx < intensity_bins; forward_intensity_idx++){
                    //double forward_intensity_used = (((float)forward_intensity_idx+1)/((float)intensity_bins));
                    //double forward_intensity_used = (((float)forward_intensity_idx)/((float)intensity_bins));
                    //double forward_expection_ion = sqrt(forward_intensity_used*forward_intensity_used - intensity_used*intensity_used);
                    //double forward_expection_ion = sqrt(forward_intensity_used - intensity_used);
                    double forward_expection_ion = intensity_difference_to_expected_ion[forward_intensity_idx-intensity_idx];
                    
                    double forward_cosine_value = min(expected_cosine + forward_expection_ion * ion_intensity, 0.9999999999);
                    
                    int forward_cosine_idx = forward_cosine_value*cosine_bins;
                    
                    double log_ratio = log(forward_expection_ion/ion_intensity)/log(2);
                    
                    double variance = ion_variance[band_idx];
                    double mean = ion_mean[band_idx];
                    
                    double current_gaussian_lookup = points_in_gaussain[forward_intensity_idx - (intensity_idx)];
                    double previous_gaussian_lookup = 0.f;
                    if(forward_intensity_idx - (intensity_idx)  == 0 ){
                        previous_gaussian_lookup = -999999999999;
                    }
                    else{
                        previous_gaussian_lookup = points_in_gaussain[forward_intensity_idx - (intensity_idx) - 1];
                    }
                    
                    
                    double normalizing_factor;
                    //double ion_variance_probability =  abs(Utils::gaussiancdf(previous_gaussian_lookup, 0.f, variance) - Utils::gaussiancdf(current_gaussian_lookup, 0.f, variance));
                    
                    double ion_variance_probability;
                    
                    //Switching on model of ion variance
                    switch(training_model){
                        case 1://Gaussian Zero'd Means
                            normalizing_factor = Utils::gaussiancdf(max_gaussian_lookup, 0.f, variance);
                            ion_variance_probability =  abs(Utils::gaussiancdf(previous_gaussian_lookup, 0.f, variance) - Utils::gaussiancdf(current_gaussian_lookup, 0.f, variance));
                            break;
                        case 2://Gaussian actual means
                            normalizing_factor = Utils::gaussiancdf(max_gaussian_lookup, mean, variance);
                            ion_variance_probability =  abs(Utils::gaussiancdf(previous_gaussian_lookup, mean, variance) - Utils::gaussiancdf(current_gaussian_lookup, mean, variance));
                            break;
                        case 3://Empirical with mean offset
                            normalizing_factor = ion_variance_cdf(histograms, band_idx, max_gaussian_lookup, mean);
                            ion_variance_probability =  abs(ion_variance_cdf(histograms, band_idx, previous_gaussian_lookup, mean) - ion_variance_cdf(histograms, band_idx, current_gaussian_lookup, mean));
                            break;
                        case 4://Empirical without mean offset
                            //normalizing_factor = ion_variance_cdf(histograms, band_idx, max_gaussian_lookup, 0.f);
                            //ion_variance_probability =  abs(ion_variance_cdf(histograms, band_idx, previous_gaussian_lookup, 0.f) - ion_variance_cdf(histograms, band_idx, current_gaussian_lookup, 0.f));
                            ion_variance_probability = ion_variance_integration(histograms, band_idx, previous_gaussian_lookup, current_gaussian_lookup, 0.f);
                            break;
                        default:
                            normalizing_factor = 0.f;
                            ion_variance_probability = 0.f;
                    }
                    
                    
                    
                    double ion_variance_probability_normalized = ion_variance_probability;
                    if(ion_variance_probability == 0.f) continue;
                        
                    /*if(normalizing_factor > 0.0001){
                        //ion_variance_probability_normalized = ion_variance_probability / normalizing_factor;
                        ion_variance_probability_normalized = ion_variance_probability;
                    }
                    else{
                        ion_variance_probability_normalized = 0.f;
                    }*/
                    
                    double specific_nondeleted_probability = (1.f - specific_deletion_probability);
                    
                    current_probability = get_3d_array(intensity_bins, cosine_bins, ion_size, intensity_idx, cosine_idx, ion_idx, programming_matrix);
                    
                    double propogating_probability = current_probability * specific_nondeleted_probability * ion_variance_probability_normalized;
                    
                    if(propogating_probability != propogating_probability){
                        DEBUG_MSG("NAN at ion\t"<<ion_idx<<"\tintensity\t"<<forward_intensity_idx<<"\tcosine\t"<<forward_cosine_idx<<"\t"<<normalizing_factor);
                    }
                    
                    double forward_probability = get_3d_array(intensity_bins, cosine_bins, ion_size, forward_intensity_idx, forward_cosine_idx, ion_idx+1, programming_matrix);
                    
                    propogating_probability += forward_probability;
                    
                    set_3d_array(intensity_bins, cosine_bins, ion_size, forward_intensity_idx, forward_cosine_idx, ion_idx+1, programming_matrix, propogating_probability);
                    
                    if(propogating_probability > 0.0){
                        
                     //   cout<<"DEBUG\t"<<ion_idx<<"\t"<<"intensity_idx"<<"\t"<<"cosine_idx"<<"\t"<<"forward_intensity_idx"<<"\t"<<"forward_intensity_used"<<"\t";
                     //   cout<<"current_probability"<<"\t"<<"forward_probability"<<"\t"<<"propogating_probability"<<"\t"<<"ion_variance_probability"<<"\t"<<"current_gaussian_lookup"<<"\t"<<"previous_gaussian_lookup"<<endl;
                     //   cout<<"DEBUG\t"<<ion_idx<<"\t"<<intensity_idx<<"\t"<<cosine_idx<<"\t"<<forward_intensity_idx<<"\t"<<forward_intensity_used<<"\t";
                     //   cout<<current_probability<<"\t"<<forward_probability<<"\t"<<propogating_probability<<"\t"<<ion_variance_probability<<"\t"<<current_gaussian_lookup<<"\t"<<ion_variance_cdf(histograms, band_idx, current_gaussian_lookup, 0.f)<<"\t"<<previous_gaussian_lookup<<"\t"<<ion_variance_cdf(histograms, band_idx, previous_gaussian_lookup, 0.f)<<endl;
                      //cout<<ion_idx<<"\t"<<intensity_idx<<"\t"<<cosine_idx<<"\t"<<forward_intensity_idx<<"\t"<<forward_intensity_used<<"\t";
                      //cout<<current_gaussian_lookup<<"\t"<<ion_variance_probability<<"\t"<<normalizing_factor<<"\t"<<ion_variance_probability_normalized<<"\t"<<current_probability<<"\t"<<propogating_probability<<"\t"<<forward_cosine_idx<<endl;
                      //cout<<"NondeleteProb\t"<<specific_nondeleted_probability<<"\tIONVARIANCEPROB\t"<<ion_variance_probability_normalized<<"\t"<<current_gaussian_lookup<<"\t"<<Utils::gaussiancdf(current_gaussian_lookup, mean, variance)<<"\t"<<previous_gaussian_lookup<<"\t"<<Utils::gaussiancdf(previous_gaussian_lookup, mean, variance)<<"\t"<<variance<<"\t"<<forward_expection_ion<<"\t"<<ion_intensity<<"\t"<<log_ratio<<endl;
                      //cout<<"BandIDX\t"<<band_idx<<endl;
                    }
                    
                }
                
            }
        }
    }
    
    bool display_output = false;
    
    //DEBUG
    for(int i = 0; i < ion_mass_intensity_pair_vector.size(); i++){
        //cout<<"IONINT\t"<<ion_mass_intensity_pair_vector[i]<<endl;
    }
    
    if(display_output){
        int resize_dim_cosine = 100;
        int resize_dim_intensity = 100;
        
        int resize_factor_cosine = cosine_bins/resize_dim_cosine;
        int resize_factor_intensity = intensity_bins/resize_dim_intensity;
        
        for(int ion_idx = 0; ion_idx < ion_size; ion_idx++){
            double * resize_matrix = new double[resize_dim_cosine*resize_dim_intensity*sizeof(double)];
            
            for(int cosine_idx = 0; cosine_idx < resize_dim_cosine; cosine_idx++){
                for(int intensity_idx = 0; intensity_idx < resize_dim_intensity; intensity_idx++){
                    set_3d_array(resize_dim_intensity, resize_dim_cosine, 1, intensity_idx, cosine_idx, 0, resize_matrix, 0.f);
                }
            }
            
            
            cout<<"ION IDX\t"<<"0.0000"<<ion_idx<<"\t"<<ion_mass_intensity_pair_vector[ion_idx]<<endl;
            for(int cosine_idx = 0; cosine_idx < cosine_bins; cosine_idx++){
                for(int intensity_idx = 0; intensity_idx < intensity_bins; intensity_idx++){
                    //cout<<get_3d_array(intensity_bins, cosine_bins, ion_size, intensity_idx, cosine_idx, ion_idx, programming_matrix)<<"\t";
                    
                    double current = get_3d_array(resize_dim_intensity, resize_dim_cosine, 1, intensity_idx/resize_factor_intensity, cosine_idx/resize_factor_cosine, 0, resize_matrix);
                    double additional = get_3d_array(intensity_bins, cosine_bins, ion_size, intensity_idx, cosine_idx, ion_idx, programming_matrix);
                    current += additional;
                    set_3d_array(resize_dim_intensity, resize_dim_cosine, 1, intensity_idx/resize_factor_intensity, cosine_idx/resize_factor_cosine, 0, resize_matrix, current);
                }
                //cout<<endl;
            }
            //cout<<endl;
            
            for(int cosine_idx = 0; cosine_idx < resize_dim_cosine; cosine_idx++){
                for(int intensity_idx = 0; intensity_idx < resize_dim_intensity; intensity_idx++){
                    cout<<get_3d_array(resize_dim_intensity, resize_dim_cosine, 1, intensity_idx, cosine_idx, 0, resize_matrix)<<"\t";
                }
                cout<<endl;
            }
            cout<<endl;
            
            
            
            delete resize_matrix;
        }
    }
    
    double total_prob = 0.f;
    for(int cosine_idx = 0; cosine_idx< cosine_bins; cosine_idx++){
        total_prob += get_3d_array(intensity_bins, cosine_bins, ion_size, intensity_bins-1, cosine_idx, ion_size-1, programming_matrix);
    }
    
    if(total_prob != 0.f){
        for(int cosine_idx = 0; cosine_idx< cosine_bins; cosine_idx++){
            float cosine_bin = cosine_idx * (1.f/cosine_bins);
            float normalized_bin_prob = get_3d_array(intensity_bins, cosine_bins, ion_size, intensity_bins-1, cosine_idx, ion_size-1, programming_matrix) / total_prob;
            pair<float, float> cosine_prob_pair;
            cosine_prob_pair.first = cosine_bin;
            cosine_prob_pair.second = normalized_bin_prob;
            
            output_distribution.push_back(cosine_prob_pair);
            output_stream<<cosine_bin<<"\t"<<normalized_bin_prob<<endl;
        }
    }
    
    
    delete programming_matrix;
    
}

double get_3d_array(int x_dim, int y_dim, int z_dim, int x_addr, int y_addr, int z_addr, double * ptr){
    int index = x_addr + y_addr * x_dim + z_addr * x_dim * y_dim;
    return ptr[index];
}

int set_3d_array(int x_dim, int y_dim, int z_dim, int x_addr, int y_addr, int z_addr, double * ptr, double value){
    int index = x_addr + y_addr * x_dim + z_addr * x_dim * y_dim;
    ptr[index] = value;
    return 0;
}

vector<float> sort_descending(vector<float> input){
    sort(input.begin(), input.end());
    return input;
}


double ion_variance_cdf(vector<vector<pair<float, float> > > &histograms, int bandidx, float x, float mean_adjustment){
    double cdf = 0.f;
    
    for(int i = 0; i < histograms[bandidx].size(); i++){
        if((x + mean_adjustment) > histograms[bandidx][i].first){
            cdf += histograms[bandidx][i].second;
        }
        else{
            return cdf;
        }
    }
    return cdf;
}

double ion_variance_integration(vector<vector<pair<float, float> > > &histograms, int bandidx, double left, double right, double mean_adjustment){
    double integration = 0.f;
    
    if(left > 4.0 && right > 4.0){
        return 0.f;
    }
    
    for(int i = 0; i < histograms[bandidx].size(); i++){
        if((right + mean_adjustment) > histograms[bandidx][i].first){
            if((left + mean_adjustment) <= histograms[bandidx][i].first){
                integration += histograms[bandidx][i].second;
            }
        }
        else{
            return integration;
        }
    }
    return integration;
}



void deleteIons(psmPtr psm, int peptideLength, MS2ScoringModel &model, vector<string> &ionsToExtract) {
    if(psm->m_peakAnnotations.size()!=psm->m_spectrum->size()) return;

    // Select target ions
    unsigned int baseIdx=0;
    for(unsigned int ionType=0; ionType<ionsToExtract.size(); ionType++) {
    for(unsigned int peakIdx=0; peakIdx<psm->m_peakAnnotations.size(); peakIdx++)
        if(psm->m_peakAnnotations[peakIdx].first and psm->m_peakAnnotations[peakIdx].first->name==ionsToExtract[ionType]){
            //ions[baseIdx+psm->m_peakAnnotations[peakIdx].second-1]=(*psm->m_spectrum)[peakIdx][1];
            (*psm->m_spectrum)[peakIdx][1] = 0.f;
        }
    }
}


void save_cosine_distribution(vector<vector<pair<float, float> > > histograms, string filename){
    ofstream outfile (filename.c_str());
    
    outfile<<"IntensityBands\t"<<histograms.size()<<endl;
    cout<<"IntensityBands\t"<<histograms.size()<<endl;
    for(int i = 0; i < histograms[0].size(); i++){
        outfile<<histograms[0][i].first<<"\t";
        cout<<histograms[0][i].first<<"\t";
        for(int j = 0; j < histograms.size(); j++){ 
            outfile<<histograms[j][i].second<<"\t";
            cout<<histograms[j][i].second<<"\t";
        }
        outfile<<endl;
        cout<<endl;
    }
}


void load_cosine_distribution(vector<vector<pair<float, float> > > & histograms, string filename){
    ifstream infile (filename.c_str());
    string temp;
    int bands;
    
    infile>>temp>>bands;
    
    histograms.clear();
    
    for(int i = 0; i < bands; i++){
        vector<pair<float, float> > histogram;
        histograms.push_back(histogram);
    }
    
    while (infile.good()){
        float histogram_first;
        infile>>histogram_first;
        for(int j = 0; j < bands; j++){
            float probability;
            infile>>probability;
            
            pair<float, float> histogram_pair;
            histogram_pair.first = histogram_first;
            histogram_pair.second = probability;
            
            histograms[j].push_back(histogram_pair);
        }
    }   
}

void save_deletion_distribution(vector < pair < float, float > > &deletion_prob, string filename){
    ofstream outfile (filename.c_str());
    
    for(int i = 0; i < deletion_prob.size(); i++){
        outfile<<deletion_prob[i].first<<"\t"<<deletion_prob[i].second<<endl;
    }
}


void load_deletion_distribution(vector < pair < float, float > > &deletion_prob, string filename){
    ifstream infile (filename.c_str());

    while (infile.good()){
        float delete_first;
        float delete_prob;

        infile>>delete_first>>delete_prob;

        cout<<delete_first<<"\t"<<delete_prob<<endl;

        pair < float, float > deletion_pair;
        deletion_pair.first = delete_first;
        deletion_pair.second = delete_prob;

        deletion_prob.push_back(deletion_pair);
    }
}


void save_SLGF(SpectralLibrary lib, string filename){
    vector<vector<float> > AllSLGF;
    int SLGF_size = 0;
    for(int lib_idx = 0; lib_idx < lib.size(); lib_idx++){
        vector<float> current_SLGF;
        for(int SLGF_idx = 0; SLGF_idx < lib[lib_idx].psmList.front()->SLGF_distribution.size(); SLGF_idx++){
            current_SLGF.push_back(lib[lib_idx].psmList.front()->SLGF_distribution[SLGF_idx].second);
        }
        SLGF_size = max(SLGF_size, (int)current_SLGF.size());
        
        if(current_SLGF.size() == 0){
            current_SLGF.resize(SLGF_size, 0.f);
        }
        
        DEBUG_MSG(lib_idx<<"\tCOUNT\t"<<current_SLGF.size()<<"\t"<<lib[lib_idx].psmList.front()->m_annotation);
        AllSLGF.push_back(current_SLGF);
    }
    
    Save_binArray( filename.c_str(), AllSLGF);
    DEBUG_MSG("SAVED_BIN\t" << filename);
}

void load_SLGF(SpectralLibrary &lib, string filename){
    vector<vector<float> > AllSLGF;
    Load_binArray( filename.c_str(), AllSLGF);
    //DEBUG_MSG("SLGF SIZE:\t"<<AllSLGF.size());
    for(int lib_idx = 0; lib_idx < lib.size(); lib_idx++){

        int SLGF_size = AllSLGF[lib_idx].size();
        //DEBUG_MSG("SLGF LIB SIZE\t"<<SLGF_size);
        float increment_amount = 1.0 / (float)(SLGF_size);
        
        lib[lib_idx].psmList.front()->SLGF_distribution.clear();
        for(int SLGF_idx = 0; SLGF_idx < SLGF_size; SLGF_idx++){
            pair<float, float> SLGF_pair ( SLGF_idx * increment_amount, AllSLGF[lib_idx][SLGF_idx]);
            lib[lib_idx].psmList.front()->SLGF_distribution.push_back(SLGF_pair);
            //DEBUG_MSG(SLGF_pair.first<<"\t"<<SLGF_pair.second);
        }
        
        //DEBUG_MSG(lib_idx<<"\tCOUNT\t"<<SLGF_size<<"\t"<<lib[lib_idx].psmList.front()->m_annotation);
        
    }
}

float SLGF_rescore(vector<pair < float, float> > &cosine_distribution, float cosine){
    //SLGF Correction
    //Finding mean
    
    float mean = 0.f;
    for(int i = 0; i < cosine_distribution.size(); i++){
        mean += cosine_distribution[i].first * cosine_distribution[i].second;
    }
    float distance_to_one = 1.f - mean;
    
    float variance = 0.f;
    for(int i = 0; i < cosine_distribution.size(); i++){
        variance += (cosine_distribution[i].first  - mean) * (cosine_distribution[i].first  - mean) * cosine_distribution[i].second;
    }
    
    
    float shift_to_mean = mean - cosine;
    float shift_factor = 800.0;
    
    //cosine += shift_to_mean * shift_factor * distance_to_one * distance_to_one * distance_to_one; 
    //DEBUG_MSG("MEAN: "<<mean<<"\tVAR: "<<variance);
    
    
    
    //Alternate offset
    //float offset = (pow((distance_to_one*100)/5.0, 1.7)/100.f);
    
    //float offset = distance_to_one/1.5;//seemed good for low and high abundance
    //float offset = distance_to_one/4; 
    
    //cosine += offset;
    //cosine += 0.04;
    
    float integrate_from_left = 0.f;
    for(int i = 0; i < cosine_distribution.size(); i++){
        if(cosine_distribution[i].first < cosine){
            integrate_from_left += cosine_distribution[i].second;
        }
    }
    
    /*if(integrate_from_left > 0.1){
        DEBUG_MSG("COS\t"<<cosine);
        for(int i = 0; i < cosine_distribution.size(); i++){
            DEBUG_MSG(cosine_distribution[i].first<<"\t"<<cosine_distribution[i].second);
        }
        DEBUG_MSG("SLGF\t"<<integrate_from_left);
    }*/
    
    
    
    return integrate_from_left;
}

string strip_extension(string input_string){
    int lastindex = input_string.find_last_of("."); 
    string rawname = input_string.substr(0, lastindex); 
    return rawname;
}



void qtof_filtering_spectrum(SpecSet &specs, const ParameterList & m_params){
    
    SpecSet * m_outputSpecset = & specs;
    float snr_thresh = m_params.getValueFloat("FILTER_SNR_PEAK_INT", 0.0);
    if (snr_thresh > 0.0){
        DEBUG_MSG("FILTERING SNR PEAK" << snr_thresh);
        for (int i = 0; i < specs.size(); i++){
            list<int> peaksToRemove;
            float noise_level = specs[i].getNoiseLevel();
            float intensity_threshold = noise_level * m_params.getValueFloat("FILTER_SNR_PEAK_INT", 0.0);
            for (int peak_idx = 0; peak_idx < specs[i].size();peak_idx++){
                if (specs[i][peak_idx][1] < intensity_threshold){
                    peaksToRemove.push_back(peak_idx);
                }
            }
            
            specs[i].removePeaks(peaksToRemove);
        }   
    }
    
    
    if (m_params.getValueInt("FILTER_PRECURSOR_WINDOW", 0) > 0)
    {
      DEBUG_MSG("Removing peaks near precursor");
      for (int i = 0; i < m_outputSpecset->size(); i++)
      {
        //Removing peaks around precursor
        float precursor_mz = (*m_outputSpecset)[i].parentMZ;
        list<int> peaksToRemove;
        for (int peak_idx = 0; peak_idx < (*m_outputSpecset)[i].size();
            peak_idx++)
        {
          if ((*m_outputSpecset)[i][peak_idx][0] > (precursor_mz - 17)
              && (*m_outputSpecset)[i][peak_idx][0] < (precursor_mz + 15))
          {
            //(*m_outputSpecset)[i][peak_idx][1] = 0.f;
            peaksToRemove.push_back(peak_idx);
          }
        }
        
        (*m_outputSpecset)[i].removePeaks(peaksToRemove);
      } 
    }
    
    
    if (m_params.getValueFloat("FILTER_STDDEV_PEAK_INT", 0.0) > 0.0)
    {
      DEBUG_MSG("Removing peaks below " << m_params.getValueFloat("FILTER_STDDEV_PEAK_INT", 0.0) <<" STD DEV ABOVE MEAN" );
      
      
      for (int i = 0; i < m_outputSpecset->size(); i++)
        {
        
        list<int> peaksToRemove;
        
        //Removing peaks lower than min_peak_intensity
        peaksToRemove.clear();
        
      
        //Calculating Mean
        float total_int = 0.0;
        int total_count = 0;
        vector<float> peak_intensity_list;
        
        for (int peak_idx = 0; peak_idx < (*m_outputSpecset)[i].size();
            peak_idx++)
        {
          peak_intensity_list.push_back((*m_outputSpecset)[i][peak_idx][1]);
        }
        
        //sorting and finding bottom 25%
        sort(peak_intensity_list.begin(), peak_intensity_list.end());
        if(peak_intensity_list.size() < 10) {
            continue;
        }
        
        //float bottom_threshold = peak_intensity_list[peak_intensity_list.size()/4];
        for (int peak_idx = 0; peak_idx < peak_intensity_list.size()/4; peak_idx++){
            total_int+= peak_intensity_list[peak_idx];
            total_count++;
        }
        
        float mean = total_int / total_count;
        float sum_variance = 0.0;
        for (int peak_idx = 0; peak_idx < peak_intensity_list.size()/4; peak_idx++){
            sum_variance += (mean - peak_intensity_list[peak_idx]) * (mean - peak_intensity_list[peak_idx]);
        }
        float variance = sum_variance / total_count;
        float std_dev = sqrt(variance);
        
        float upper_threshold = m_params.getValueFloat("FILTER_STDDEV_PEAK_INT", 0.0) * std_dev + mean;
      
        for (int peak_idx = 0; peak_idx < (*m_outputSpecset)[i].size();
            peak_idx++)
        {
            if ((*m_outputSpecset)[i][peak_idx][1] < upper_threshold){
                peaksToRemove.push_back(peak_idx);
            }
        }
        
        (*m_outputSpecset)[i].removePeaks(peaksToRemove);
        
      }
    }
    
    if (m_params.getValueFloat("MIN_PEAK_INT", 0.0) > 0.0)
    {
      DEBUG_MSG("Removing peaks below " << m_params.getValueFloat("MIN_PEAK_INT", 0.0) );
      for (int i = 0; i < m_outputSpecset->size(); i++)
        {
        
        list<int> peaksToRemove;
        
        //Removing peaks lower than min_peak_intensity
        peaksToRemove.clear();
        float min_peak_intensity = m_params.getValueFloat("MIN_PEAK_INT", 50.0);
        for (int peak_idx = 0; peak_idx < (*m_outputSpecset)[i].size();
            peak_idx++)
        {
          if ((*m_outputSpecset)[i][peak_idx][1] < min_peak_intensity)
          {
            //(*m_outputSpecset)[i][peak_idx][1] = 0.f;
            peaksToRemove.push_back(peak_idx);
          }
        }
        
        (*m_outputSpecset)[i].removePeaks(peaksToRemove);
        
      }
    }
    
    
    if (m_params.getValueInt("WINDOW_FILTER", 0) == 1){
        //Window filtering
        DEBUG_MSG("FILTER PEAKS");
        for(int i = 0; i< specs.size(); i++){
            specs[i].rankFilterPeaks(6,50);
        }
    }
    
}

string get_extension(string input_string){
    string spectra_file_name = input_string;
    int dotpos = spectra_file_name.find_last_of('.');
    string extension = spectra_file_name.substr(dotpos+1);
    return extension;
}

int diff_molecular_psm(psmPtr psm1, psmPtr psm2){
    if( psm1->m_annotation != psm2->m_annotation ||
        psm1->m_organism != psm2->m_organism ||
        psm1->m_compound_name != psm2->m_compound_name ||
        psm1->m_smiles != psm2->m_smiles ||
        psm1->m_InChI != psm2->m_InChI ||
        psm1->m_InChI_Aux != psm2->m_InChI_Aux ||
        psm1->m_notes != psm2->m_notes ||
        psm1->m_ionmode != psm2->m_ionmode ||
        psm1->m_exactmass != psm2->m_exactmass ||
        psm1->m_library_quality != psm2->m_library_quality
        ){
        
        cout<<psm1->m_annotation <<"\t"<< psm2->m_annotation <<endl;
        cout<<psm1->m_organism <<"\t"<< psm2->m_organism <<endl;
        cout<<psm1->m_compound_name <<"\t"<< psm2->m_compound_name <<endl;
        cout<<psm1->m_smiles <<"\t"<< psm2->m_smiles <<endl;
        cout<<psm1->m_InChI <<"\t"<< psm2->m_InChI <<endl;
        cout<<psm1->m_InChI_Aux <<"\t"<< psm2->m_InChI_Aux <<endl;
        cout<<psm1->m_notes <<"\t"<< psm2->m_notes <<endl;
        cout<<psm1->m_ionmode <<"\t"<< psm2->m_ionmode <<endl;
        cout<<psm1->m_exactmass <<"\t"<< psm2->m_exactmass <<endl;
        cout<<psm1->m_library_quality <<"\t"<< psm2->m_library_quality <<endl;
        
        return 1;
    }
    
    return 0;
}



string string_replace(std::string s,
                      std::string toReplace,
                      std::string replaceWith)
{
    return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
}

string msp_annotation_to_specnets(string input_annotation){
    
    string modified_annot = "";
    //Need to change the mods
    modified_annot = string_replace_all(input_annotation, "C[160]" , "C");      //Assume Capped Cysteines.
    modified_annot = string_replace_all(modified_annot, "n" , "");              //remove n term n
    modified_annot = string_replace_all(modified_annot, "M[147]" , "(M,16)");   //M Oxidation
    modified_annot = string_replace_all(modified_annot, "Q[111]" , "(Q,-17)");              //remove n term n
    modified_annot = string_replace_all(modified_annot, "C[339]" , "(C,236)");              //remove n term n
    
    return modified_annot;
}

string specnets_annotation_to_msgf(string input_annotation){
    string modified_annot = "";
    
    modified_annot = string_replace_all(input_annotation, "(" , "");
    modified_annot = string_replace_all(modified_annot, ")" , "");
    modified_annot = string_replace_all(modified_annot, "," , "+");
    modified_annot = string_replace_all(modified_annot, "+-" , "-");
    modified_annot = string_replace_all(modified_annot, "++" , "+");
    
    modified_annot = create_annotation_ends(modified_annot);
    
    modified_annot = string_replace_all(modified_annot, ".." , ".");
    
    return modified_annot;
}

