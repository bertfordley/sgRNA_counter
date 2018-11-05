# -*- coding: utf-8 -*-
"""
By: Rob Moseley

Script for counting sgRNA in fastq files.
Uses SequenceMatcher of the difflib package to find matches between the 
sgRNA (and other parts) sequences and reads in the fastq file.
"""

import pymongo
import pprint
import json
import os
import gzip
import itertools
# !pip install --user fuzzywuzzy
# !pip install --user python-Levenshtein
import fuzzywuzzy.process
import pandas as pd
import argparse

def query_db_for_fastq_file(fastq_file_name):
    
    dbURI = 'mongodb://readonly:WNCPXXd8ccpjs73zxyBV@data-catalog.sd2e.org:27020/admin?readPreference=primary'
    client = pymongo.MongoClient(dbURI)
    db = client.catalog_dev
    science_table=db.science_table
    query={}
    query["experiment_reference"] = "YeastSTATES-gRNA-Seq-Diagnosis"
    query['lab']='Ginkgo'
    query['measurement_type'] = 'RNA_SEQ'
    query['file_type'] = 'FASTQ'
    all_sgrna = []
    for sample in science_table.find(query):
        all_sgrna.append(sample)
    
    fastq_metadata = {}
    for sample in all_sgrna:
        try:
            if sample['strain_input_state'] is not None and os.path.basename(sample["filename"]) == fastq_file_name:
                fastq_metadata = sample
        except Exception as e: continue   

    return fastq_metadata


def get_metadata(sample):
    # Get metadata for sample
    
    out_dict = {}
    
    out_dict["gate"] = sample["strain_circuit"]
    out_dict["left_signal"] = sample["strain_input_state"][0]
    out_dict["right_signal"] = sample["strain_input_state"][1]
    out_dict["R1_fastq_filename"] = os.path.basename(sample["filename"])
    out_dict["sample_id"] = sample["sample_id"]
    out_dict["reference_sample_id"] = sample["reference_sample_id"]
    out_dict["replicate"] = sample["replicate"]
    out_dict["strain_lab_id"] = sample["strain_lab_id"]
    out_dict["strain"] = sample["strain"]
        
    return out_dict


def create_logic_count_dict(meta_info_dict, part_seq_dict, circuit_dict):
    # extend metadata dictionary to include circuit information, e.g., input signal state,
    # which sgRNA are the inputs and which are the intermediates. Also add in key:value
    # pairs for when counting occurs
    
    gate_name = meta_info_dict["gate"]
    
    meta_info_dict["left_input"] = circuit_dict[gate_name]["input"][0]
    meta_info_dict["right_input"] = circuit_dict[gate_name]["input"][1]

    for i, part in enumerate(circuit_dict[gate_name]["intermediate"]):
        if part:
            meta_info_dict["intermediate" + str(i + 1)] = part
        else:
            meta_info_dict["intermediate" + str(i + 1)] = "-" 

    for part, seq in part_seq_dict.items():
        meta_info_dict[part] = 0

    return meta_info_dict


def return_gate_grna(meta_info_dict, part_seq_dict, circuit_dict):
    # make a list of appropriate sgRNA sequences for sample being analyzed to be
    # used for matching to reads in fastq file. Also make dictionary specifically for
    # the gate that includes sgRNA names and sequences for matching.
    
    gate_name = meta_info_dict["gate"]
    
    gRNA_list = []
    gate_part_dict = {}
    
    for part, seq in part_seq_dict.items():
        if part in circuit_dict[gate_name]["input"] or part in circuit_dict[gate_name]["intermediate"]:
            gRNA_list.append(seq)
            gate_part_dict[part] = seq
    
    return gRNA_list, gate_part_dict


def read_fastq_seqs(filepath):
    # unzip and read fastq file
    
    with gzip.open(filepath, 'rt', encoding='utf') as fh:
        for seq_header, seq, qual_header, qual in itertools.zip_longest(*[fh] * 4):
            
            yield seq.rstrip('\n')


def match_parts(line_1, line_2, gRNA_seqs, gate_part_dict, thres):
    # match sgRNAs in the cicuit to reads in the fastq R1 and R2 files
    # fuzzy matching is used with user picking a threshold for percent
    # match. Default is 80.
    
    line_1_rev = reverse_complement(line_1)
    line_2_rev = reverse_complement(line_2)

    goodmatches = []

    match_r1 = fuzzywuzzy.process.extractOne(line_1, gRNA_seqs)
    if match_r1[1] > thres:
        goodmatches.append(match_r1)
    
    match_r1_rev = fuzzywuzzy.process.extractOne(line_1_rev, gRNA_seqs)
    if match_r1_rev[1] > thres:
        goodmatches.append(match_r1_rev)
        
    match_r2 = fuzzywuzzy.process.extractOne(line_2, gRNA_seqs)
    if match_r2[1] > thres:
        goodmatches.append(match_r2)
    
    match_r2_rev = fuzzywuzzy.process.extractOne(line_2_rev, gRNA_seqs)
    if match_r2_rev[1] > thres:
        goodmatches.append(match_r2_rev)

    parts_list = []
    
    for part, seq in gate_part_dict.items():
        if goodmatches:
            for goodm in goodmatches:
                if seq == goodm[0]:
                    parts_list.append(part)

    return parts_list


def reverse_complement(seq): 
    # make reverese complement of read
    
    alt_map = {'ins': '0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    
    for k, v in alt_map.items():
        seq = seq.replace(k, v)
        
    bases = list(seq) 
    bases = reversed([complement.get(base, base) for base in bases])
    bases = ''.join(bases)
    
    for k, v in alt_map.items():
        bases = bases.replace(v, k)
        
    return bases


def count_file(sample, gate_seqs_dict, gate_info_dict, threshold):
    # process a single fastq file for counting sgRNA
    
    # get metadata dictionary
    meta_info_dict = get_metadata(sample)
    
    # get names for both R1 and R2 files for sample
    filename_R1 = sample["jupyter_path"]
#     filename_R1 = os.path.join(fastq_dir, os.path.basename(meta_info_dict["R1_fastq_filename"]))
    filename_R2 = filename_R1[:-10] + "2.fastq.gz"

    print("Gate {}, Input State {}{}: {} {}".format(meta_info_dict["gate"],
                                                    meta_info_dict["left_signal"],
                                                    meta_info_dict["right_signal"],
                                                    os.path.basename(filename_R1),
                                                    os.path.basename(filename_R2)))
    
    # expand metadata dictionary to include circuit info
    logic_patternCounts = create_logic_count_dict(meta_info_dict, gate_seqs_dict, gate_info_dict)

    # get specific gate sgRNA name and sequences
    gRNA_seq_list, gate_part_dict = return_gateg_rna(meta_info_dict, gate_seqs_dict, gate_info_dict)

    # read R1 and R2 fastq files and make iterators for both
    seq_gen_r1 = read_fastq_seqs(filename_R1)
    seq_gen_r2 = read_fastq_seqs(filename_R2)

    # for each line pair in R1 and R2 fastq files, match sgRNAs found in sample's gate
    for i, (line_r1, line_r2) in enumerate(zip(seq_gen_r1, seq_gen_r2)):
        parts_list = match_parts(line_r1, line_r2, gRNA_seq_list, gate_part_dict, threshold)
        if parts_list:
            for part in parts_list:
                print("Match in line {}: {}".format(i+1, part))
                logic_patternCounts[part] += 1

#     #return pd.DataFrame(logic_patternCounts, index=[1])
    return(logic_patternCounts)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(prog="sgRNA counter",
                                     description="Using fuzzy string matching to find sgRNA target sequences within a fastq file")
                            
    parser.add_argument("fastq_R1_file", 
                        action="store", 
                        help="R1 fastq file containing reads, don't include path")
                                  
    parser.add_argument("-o", "--output_dir",
                        type=str,
                        default=None,
                        required=True,
                        help="output directory to write results")                                  
                                     
    parser.add_argument("-t", "--threshold",
                        type = int,
                        default=80,
                        required=True,
                        help="threshold for string matching")  
                                     
    args = parser.parse_args()
    
    fastq_file = args.fastq_R1_file
    threshold = args.threshold
    output_dir = args.output_dir
    
    gate_dict = {"OR": {"input": ["r3", "r6"], "intermediate": ["r7", "", ""], 
                    "additional": ["HHRz5p", "HHRz3p", "sgRNA_handle", "HDVRz"]},
             "AND": {"input": ["r2", "r1"], "intermediate": ["r7", "r5", ""], 
                     "additional": ["HHRz5p", "HHRz3p", "sgRNA_handle", "HDVRz"]},
             "NAND": {"input": ["r10", "r2"], "intermediate": ["r1", "r3", "r5"], 
                      "additional": ["HHRz5p", "HHRz3p", "sgRNA_handle", "HDVRz"]},
             "XNOR": {"input": ["r2", "r10"], "intermediate": ["r7", "r3", "r1"], 
                      "additional": ["HHRz5p", "HHRz3p", "sgRNA_handle", "HDVRz"]},
             "XOR": {"input": ["r3", "r6"], "intermediate": ["r1", "r9", "r7"], 
                     "additional": ["HHRz5p", "HHRz3p", "sgRNA_handle", "HDVRz"]}}

    gate_part_seqs = {'r1': 'GGAACGTGATTGAATAACTT',
                      'r2': 'ACCAACGCAAAAAGATTTAG',
                      'r3': 'CATTGCCATACACCTTGAGG',
                      'r4': 'GAAAATCACAACTCTACTGA',
                      'r5': 'GAAGTCAGTTGACAGAGTCG',
                      'r6': 'GTGGTAACTTGCTCCATGTC',
                      'r7': 'CTTTACGTATAGGTTTAGAG',
                      'r8': 'CGCATTTCCTATTCAAACTT',
                      'r9': 'GCAACCCACAAATATCCAGT',
                      'r10': 'GTGACATAAACATTCGACTC',
                      'r11': 'GGGCAAAGAGACGCTTGTCG',
                      'r12': 'GAAGTCATCGCTTCTTGTCG',
                      'r13': 'GAGTTGACAAAGTATAACTT',
                      'r14': 'GAAGTTTCAGAATCTCGACG',
                      'r15': 'GGCTAGGATCCATCTGACTT',
                      'r16': 'GCAACCATAGACTCTCCAGG',
                      'r17': 'ACCACAACTGAGTCGAACCT',
                      'r18': 'GGGTAGCAACACTCGTACTT',
                      'r19': 'GTAAAAGATAACTCTGTTGC',
                      'r20': 'TCTACCCGAGACTCAAACGG',
                      'HHRz5p': 'GGATTCTAGAACTAGTGGATCTACAAA',
                      'HHRz3p': 'CTGATGAGTCCGTGAGGACGAAACGAGTAAGCTCGTC',
                      'sgRNA_handle': 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGG',
                      'HDVRz': 'CACCGAGTCGGTGCTTTTGGCCGGCATGGTCCCAGCCTCCTCGCTGGCGCCGGCTGGGCAACATGCTTCGGCATGGCGAATGGGACTGATACCGTCGACCTCGAGTC'}

    fastq_metadata = query_db_for_fastq_file(fastq_file)
            
    if fastq_metadata is None:
        print("Either {} isn't in database or no circuit information is available".format(fastq_file_name))
    else:   
        results = count_file(fastq_metadata, gate_part_seqs, gate_dict, threshold)
        print(pd.DataFrame(results, index=[1]).T)                           
