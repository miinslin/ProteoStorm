# -*- coding: utf-8 -*-

from Bio import SeqIO
import urllib, os, time
import xml.etree.ElementTree as ET
import shutil

def PanMicrobial_preprocessing(FASTAdir, MBsize, FASTAchunkdir):
    start_time = time.time()
    
    DBsizeBytes = 1024*1024*MBsize
    filenum= -1
    filesize = 0
    
    for fipart, fpart in enumerate([f for f in os.listdir(FASTAdir) if f.endswith('.fasta')]):
        for record in SeqIO.parse(os.path.join(FASTAdir,fpart), "fasta"):
            sequence= str(record.seq).replace('L','I')
            sequence_split = [sequence[i:i+60] for i in range(0, len(sequence), 60)]
            header = '>'+str(record.description)
            if filesize == 0:
                filenum += 1
                fastabin = os.open(os.path.join(FASTAchunkdir,'combinedDB_part_'+str(filenum)+'.fasta'),os.O_RDWR|os.O_CREAT)
            filesize+=os.write(fastabin, header+'\n'+'\n'.join(sequence_split)+'\n')
            
            if filesize>DBsizeBytes:
                os.close(fastabin)
                filesize = 0
    if filesize>0:
        os.close(fastabin)
    
    return time.time()-start_time

def downloadNCBI(ncbitaxdir, taxonomyID):
    if not os.path.exists(os.path.join(ncbitaxdir, taxonomyID + '.xml')):
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=' + taxonomyID + '&retmode=xml'
        urllib.urlretrieve(url, os.path.join(ncbitaxdir, taxonomyID + '.xml'))
        
        if not os.path.exists(os.path.join(ncbitaxdir, taxonomyID + '.xml')):
            print taxonomyID
            raise ValueError('download failed')
        
    superkingdom = ''
    genusname = ''
    kingdom = ''
    
    tree = ET.parse(os.path.join(ncbitaxdir, taxonomyID + '.xml'))
    root = tree.getroot()
    for child in root.findall('./Taxon/LineageEx/Taxon'):
        Ranks = child.find('Rank').text
        if Ranks == 'superkingdom':
            superkingdom = child.find('ScientificName').text
        if Ranks == 'kingdom':
            kingdom = child.find('ScientificName').text
        if Ranks == 'genus':
            genusname = child.find('ScientificName').text

    if superkingdom == 'Bacteria':
        if not genusname:
            os.remove(os.path.join(ncbitaxdir, taxonomyID + '.xml'))
    elif superkingdom == 'Eukaryota':
        if kingdom in ['fungi','Fungi']:
            if not genusname:
                os.remove(os.path.join(ncbitaxdir, taxonomyID + '.xml'))
    else:
        os.remove(os.path.join(ncbitaxdir, taxonomyID + '.xml'))

    return {'superkingdom':superkingdom,'kingdom':kingdom,'genusname':genusname}


def RefSeq_genera_preprocessing(S1_preprocessing_dir, ncbidir,
                                RefSeq_fasta_files,RefSeq_catalog, ProteoStormLOG):
    
    ProteoStormLOG.write('Processing RefSeq genera...'+'\n')
    refseq_taxon_ids_bac_fungi = os.path.join(S1_preprocessing_dir,'RefSeq_Bac_Fungi_taxonID_acessionIDs.txt')
    RefSeq_taxon_IDs = os.path.join(S1_preprocessing_dir,'RefSeq_Bac_Fungi_genera.txt')
    
    if not os.path.exists(refseq_taxon_ids_bac_fungi):
        bac_fungi_accessions = set()
        for fastafile in os.listdir(RefSeq_fasta_files):
            with open(os.path.join(RefSeq_fasta_files,fastafile),'r') as infile:
                for line in infile:
                    if line[0]=='>':
                        accessionID = line.strip().split()[0][1:]
                        bac_fungi_accessions.add(accessionID)
              
        with open(RefSeq_catalog,'r') as infile:
            refseq_taxids = {}
            for line in infile:
                accessionID = line.strip().split('\t')[2]
                if accessionID in bac_fungi_accessions:
                    taxid = line.strip().split('\t')[0]
                    if taxid in refseq_taxids:
                        refseq_taxids[taxid].add(accessionID)
                    elif taxid not in refseq_taxids:
                        refseq_taxids[taxid] = set([accessionID])
        
        del bac_fungi_accessions
        with open(refseq_taxon_ids_bac_fungi,'w') as outfile:
            for taxid in refseq_taxids:
                outfile.write(taxid+'\t'+';'.join(refseq_taxids[taxid])+'\n')
        
        ProteoStormLOG.write(str(len(refseq_taxids))+ ' bacteria and fungi ids...'+'\n')
        del refseq_taxids
    
    accession_genus = {}
    with open(refseq_taxon_ids_bac_fungi,'r') as infile,\
    open(RefSeq_taxon_IDs,'w') as outfile:
        for line in infile:
            ID = line.strip().split('\t')[0]
            acessionIDlist = line.strip().split('\t')[1]
            results = downloadNCBI(ncbidir, ID)
            superkingdom = results['superkingdom']
            genusname = results['genusname']
            kingdom = results['kingdom']
    
            if superkingdom == 'Bacteria':
                if not genusname:
                    continue
                outfile.write('\t'.join([ID, acessionIDlist, genusname, superkingdom])+'\n')
                acessionIDlist = acessionIDlist.split(';')
                for ac in acessionIDlist:
                    accession_genus[ac] = genusname
            elif superkingdom == 'Eukaryota':
                if kingdom in ['fungi','Fungi']:
                    if not genusname:
                        continue
                    outfile.write('\t'.join([ID, acessionIDlist, genusname, kingdom])+'\n')
                    acessionIDlist = acessionIDlist.split(';')
                    for ac in acessionIDlist:
                        accession_genus[ac] = genusname
    
    return accession_genus

def Genera_preprocessing(FASTAdir, S1_preprocessing_dir, \
                                FASTAchunkdir, MBsize,\
                                RefSeq_catalog, ProteoStormLOG):
    
    ProteoStormLOG.write('Processing UniProt genera...'+'\n')
    start_time = time.time()
    
    # check for uniprot and refseq directories
    uniprot_fasta_files = os.path.join(FASTAdir,'UniProt_fasta')
    RefSeq_fasta_files = os.path.join(FASTAdir,'RefSeq_fasta')
    
    temp_UPID_TAXONID_dir = os.path.join(S1_preprocessing_dir,'temp_UPIDtaxonID')
    ncbidir = os.path.join(S1_preprocessing_dir, 'NCBI_taxonomy_xml')
    os.makedirs(temp_UPID_TAXONID_dir)
    os.makedirs(ncbidir)
    
    Uniprot_UPID_taxonID_mapping = os.path.join(S1_preprocessing_dir,'UniProtID_taxonID_mapping.txt')
    Uniprot_taxonid_UPID_genus = os.path.join(S1_preprocessing_dir,'UniProt_taxonID_UPID_genus_mappings.txt')
    
    fastafiles = ['id%3A'+x[:-6].split('_')[0] for x in os.listdir(uniprot_fasta_files)]
    # download ncbi .xml files, parse for genus
    # query cannot be too long, so limiting to 300 queries per iteration
    f_idx = 0
    for i in range(0,len(fastafiles)+300,300):
        IDlist = fastafiles[i:i+300]
        if IDlist:
            f_idx+=1
            tab_upid_taxonid = "http://www.uniprot.org/proteomes/?compress=no&query="+"+OR+".join(IDlist)+\
            "&sort=score&force=no&format=tab&columns=id,organism-id"
            urllib.urlretrieve(tab_upid_taxonid, os.path.join(temp_UPID_TAXONID_dir, 'UPID_TAXONID_'+str(f_idx)+'.txt'))
    
    mappingtempfiles = os.listdir(temp_UPID_TAXONID_dir)
    
    with open(Uniprot_UPID_taxonID_mapping,'w') as outfile:
        outfile.write('Proteome ID'+'\t'+'Organism ID'+'\n')
        for f in mappingtempfiles:
            with open(os.path.join(temp_UPID_TAXONID_dir,f),'r') as infile:
                lines = [z.strip() for z in infile.readlines()[1:]]
            outfile.write('\n'.join(lines)+'\n')

    shutil.rmtree(temp_UPID_TAXONID_dir) 

    # download ncbi files
    UPID_genus_map = {}
    with open(Uniprot_UPID_taxonID_mapping,'r') as infile,\
    open(Uniprot_taxonid_UPID_genus,'w') as outfile:
        for index, line in enumerate(infile):
            if index==0:
                continue
            ID = line.strip().split('\t')[1]
            UPID = line.strip().split('\t')[0]
            
            results = downloadNCBI(ncbidir, ID)
            superkingdom = results['superkingdom']
            genusname = results['genusname']
            kingdom = results['kingdom']
    
            if superkingdom == 'Bacteria':
                if not genusname:
                    continue
                outfile.write('\t'.join([ID, UPID, genusname, superkingdom])+'\n')
                UPID_genus_map[UPID] = genusname
            elif superkingdom == 'Eukaryota':
                if kingdom in ['fungi','Fungi']:
                    if not genusname:
                        continue
                    outfile.write('\t'.join([ID, UPID, genusname, kingdom])+'\n')
                    UPID_genus_map[UPID] = genusname
    
    os.remove(Uniprot_UPID_taxonID_mapping)   
    
    DBsizeBytes = 1024*1024*MBsize
    filenum= -1
    filesize = 0
    
    for fastafile in os.listdir(uniprot_fasta_files):
        UPID = fastafile[:-6].split('_')[0]
        if UPID in UPID_genus_map:
            genus = UPID_genus_map[UPID]
            for record in SeqIO.parse(os.path.join(uniprot_fasta_files, fastafile), "fasta"):
                sequence= str(record.seq).replace('L','I')
                sequence_split = [sequence[i:i+60] for i in range(0, len(sequence), 60)]
                header = '>'+str(record.id)+'\t'+genus
                
                if filesize == 0:
                    filenum += 1
                    fastabin = os.open(os.path.join(FASTAchunkdir,'combinedDB_part_'+str(filenum)+'.fasta'),os.O_RDWR|os.O_CREAT)
                filesize+=os.write(fastabin, header+'\n'+'\n'.join(sequence_split)+'\n')
                
                if filesize>DBsizeBytes:
                    os.close(fastabin)
                    filesize = 0

    if RefSeq_catalog!='na':
        ref_results = RefSeq_genera_preprocessing(S1_preprocessing_dir, ncbidir,
                                    RefSeq_fasta_files,RefSeq_catalog, ProteoStormLOG)
        
        for fastafile in os.listdir(RefSeq_fasta_files):
            for record in SeqIO.parse(os.path.join(RefSeq_fasta_files,fastafile), "fasta"):
                acessionID = str(record.id)
                if acessionID not in ref_results:
                    continue
                genus = ref_results[acessionID]
                
                sequence= str(record.seq).replace('L','I')
                sequence_split = [sequence[i:i+60] for i in range(0, len(sequence), 60)]
                header = '>'+acessionID+'\t'+genus
    
                if filesize == 0:
                    filenum += 1
                    fastabin = os.open(os.path.join(FASTAchunkdir,'combinedDB_part_'+str(filenum)+'.fasta'),os.O_RDWR|os.O_CREAT)
                filesize+=os.write(fastabin, header+'\n'+'\n'.join(sequence_split)+'\n')
                
                if filesize>DBsizeBytes:
                    os.close(fastabin)
                    filesize = 0
        
    if filesize>0:
        os.close(fastabin)
    
    ProteoStormLOG.write('Created fasta chunks with genera information in header...'+'\n')
    return time.time()-start_time