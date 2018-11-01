# -*- coding: utf-8 -*-
"""
By: Rob Moseley
Spyder Editor

Script for counting sgRNA in fastq files.
Uses SequenceMatcher of the difflib package to find matches between the 
sgRNA (and other parts) sequences and reads in the fastq file.
"""

import pandas as pd
import os
import gzip
import itertools
from joblib import Parallel, delayed
import fuzzywuzzy.process


def reverse_complement(seq): 
    
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


def getGatePartSeqs(gate_seqs_file):

    one_gate_seqs_dict = {}
    
    with open(gate_seqs_file) as f:
        seqs = [line.rstrip('\n') for line in f]
        for s in seqs:
            name, seq = s.split(":")
            one_gate_seqs_dict[name] = seq
        
    return one_gate_seqs_dict


def read_fastq_seqs(filepath):

    with gzip.open(filepath, 'rt', encoding='utf') as fh:
        for seq_header, seq, qual_header, qual in itertools.zip_longest(*[fh] * 4):
#             if any(line is None for line in (seq_header, seq, qual_header, qual)):
#                 raise Exception(
#                     "Number of lines in FASTQ file must be multiple of four "
#                     "(i.e., each record must be exactly four lines long).")
#             if not seq_header.startswith('@'):
#                 raise Exception("Invalid FASTQ sequence header: %r" % seq_header)
#             if qual_header != '+\n':
#                 raise Exception("Invalid FASTQ quality header: %r" % qual_header)
#             if qual == '\n':
#                 raise Exception("FASTQ record is missing quality scores.")

            yield seq.rstrip('\n')


def createLogicCountDict(filename, gate_name, seq_dict, circuit_dict):

    logic_count_dict = {}
    logic_count_dict["fastq_R1_filename"] = filename
    logic_count_dict["left_input"] = circuit_dict[gate_name]["input"][0]
    logic_count_dict["right_input"] = circuit_dict[gate_name]["input"][1]

    for i, part in enumerate(circuit_dict[gate_name]["intermediate"]):
        logic_count_dict["intermediate" + str(i + 1)] = part

    for part, seq in seq_dict.items():
        logic_count_dict[part] = 0

    return logic_count_dict


def returnGategRNA(gatename, seq_dict, gate_dict):
    
    gRNA_list = []
    gate_part_dict = {}
    
    for part, seq in seq_dict.items():
        if part in gate_dict[gatename]["input"] or part in gate_dict[gatename]["intermediate"]:
            gRNA_list.append(seq)
            gate_part_dict[part] = seq
    
    return gRNA_list, gate_part_dict


def matchParts(line_1, line_2, gRNA_seqs, gate_part_dict, thres):

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


def countFile(onefile, wdir, meta_df, gate_seqs_dict, gate_dict, idx, num_batches, batch_num, num_files, thres):
    
    if onefile == None: 
        return {}
    else:
        filename_R1 = os.path.join(wdir, onefile)
        filename_R2 = filename_R1[:-10] + "2.fastq.gz"
        print("File {} of {} in Batch {} of {}: {}, {}".format(idx + 1, num_files, batch_num + 1, num_batches, os.path.basename(filename_R1), os.path.basename(filename_R2)))
        
        gate = meta_df.loc[meta_df["fastq_R1_filename"] == onefile]["gate"].values[0]
    
        logic_patternCounts = createLogicCountDict(onefile, gate, gate_seqs_dict, gate_dict)
        
        gRNA_seq_list, gate_part_dict = returnGategRNA(gate, gate_seqs_dict, gate_dict)
        
        seq_gen_r1 = read_fastq_seqs(filename_R1)
        seq_gen_r2 = read_fastq_seqs(filename_R2)
        
        for (line_r1, line_r2) in zip(seq_gen_r1, seq_gen_r2):
            parts_list = matchParts(line_r1, line_r2, gRNA_seq_list, gate_part_dict, thres)
            if parts_list:
                for part in parts_list:
                    logic_patternCounts[part] += 1
    
        #return pd.DataFrame(logic_patternCounts, index=[1])
        return(logic_patternCounts)

def grouper(n, iterable, fillvalue=None):
    "grouper(3, ('cat','dog','fish','human'), 'x') --> ('cat','dog','fish'), ('human',None,None)"
    args = [iter(iterable)] * n
    return itertools.zip_longest(fillvalue=fillvalue, *args)


def getMeta(meta_file_name, wdir, media):
    
    meta_file = pd.read_csv(os.path.join(wdir, meta_file_name), index_col=0)
    meta_file_sliced = meta_file.loc[(meta_file["Media"] == media) & 
                               (meta_file["fastq_R1_filename"].str.contains("gz")) & 
                               (~meta_file["gate"].str.contains("wt", na=False))]
        
    return meta_file_sliced

    
if __name__ == "__main__":
    # Set input variables
    
    # Where meta file and gate seqs file are stored
    maindir = "/Users/haase-admin/Desktop/rmoseley/sgRNA_counter"
    
    # Where fastq files are stored
    filedir = "/Users/haase-admin/Desktop/rmoseley/sgRNA_counter/fastq_files"
    output_file_name = "merged_meta_counts_data"
    
    # test_file_names = ["9918080-R1.fastq.gz", "9918096-R1.fastq.gz"]
    # test_file_names = ["10004578-R1.fastq.gz", "10004562-R1.fastq.gz"]
     
    # testdir = "/Users/haase-admin/Desktop/rmoseley/test"
    
    media_selection = "standard"
    meta_file_name = "full_meta_merge.csv"
    gate_seq_file = "gate_seqs.txt"
    
    # Grab gate part sequences and meta data based on media of interest. 
    # Then make list of fastq files from meta data for sgRNA counting
    
    # gate part sequences (dictionary)
    gate_seqs_dict = getGatePartSeqs(os.path.join(maindir, gate_seq_file))
    
    # meta data file for 'standard' media (dataframe)
    meta_file_stad = getMeta(meta_file_name, maindir, media_selection)
    
    # list of fastq file names (list)
    full_file_list = meta_file_stad["fastq_R1_filename"].tolist()
    # full_file_list = ["10004550-R1.fastq.gz", "9918101-R1.fastq.gz", "9249509-R1.fastq.gz", "10003155-R1.fastq.gz"]

    # Form ciruits. First two indices in values are the inputs for the respective gate.
    # First index is left-side of truth table and second index is right-side of truth table
    # from Gander et al. 2017. For example, "OR" gate's input is "r3" and "r6".
    # For 0,1 in the truth table, "r3" is turned off and "r6" is turned on.
    
    #                   "input": [left, right]
    gate_dict = {"OR": {"input": ["r3", "r6"], 
                               "intermediate": ["r7", "", ""], 
                               "additional": ["HHRz5p", "HHRz3p", "sgRNA_handle", "HDVRz"]},
                       "AND": {"input": ["r2", "r1"], 
                               "intermediate": ["r7", "r5", ""], 
                               "additional": ["HHRz5p", "HHRz3p", "sgRNA_handle", "HDVRz"]},
                       "NAND": {"input": ["r10", "r2"], 
                                "intermediate": ["r1", "r3", "r5"], 
                                "additional": ["HHRz5p", "HHRz3p", "sgRNA_handle", "HDVRz"]},
                       "XNOR": {"input": ["r2", "r10"], 
                                "intermediate": ["r7", "r3", "r1"], 
                                "additional": ["HHRz5p", "HHRz3p", "sgRNA_handle", "HDVRz"]},
                       "XOR": {"input": ["r3", "r6"], 
                               "intermediate": ["r1", "r9", "r7"], 
                               "additional": ["HHRz5p", "HHRz3p", "sgRNA_handle", "HDVRz"]}}


    # Parallel
    # Count parts of each gate in respective fastq file
    counts_dict_list = []
    n_jobs = 10
    total_batches = int(round(len(full_file_list)/n_jobs, 0) + 1)
    sub_lists = grouper(n_jobs, full_file_list)
    for batchi, sub_file_list in enumerate(sub_lists):
        total_counts_df_list = Parallel(n_jobs=n_jobs, backend="multiprocessing")(delayed(countFile)(onefile, filedir, meta_file_stad, gate_seqs_dict, gate_dict, i, total_batches, batchi, len(sub_file_list), 80) for i, onefile in enumerate(sub_file_list))
        total_counts_df_list = [x for x in total_counts_df_list if x != {}]
        
        counts_dict_list.append(pd.DataFrame(total_counts_df_list))

        merged_meta_counts_df = pd.merge(meta_file_stad, pd.concat(counts_dict_list), how="outer", 
                                         left_on="fastq_R1_filename", right_on="fastq_R1_filename")
        
        merged_meta_counts_df.to_csv(os.path.join(maindir, output_file_name + ".csv"))
    
    print("FINSIHED!!")
    