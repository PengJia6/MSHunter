#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : benchmark.py
# Author : Peng Jia
# Date   : 2020.07.13
# Email  : pengjia@stu.xjtu.edu.cn
# Description: Build benchmark for microsatellite mutation calling
=============================================================================="""
import os
import re
import collections
import pysam
import multiprocessing
from src.global_dict import *
from src.units import *
from src.Window import Window

def benchmark_init(args):
    """
    argument procress
    """
    paras = {
        "command": "benchmark",
        "input": args.input[0],
        "output": args.output[0],
        "microsatellite": args.microsatellite[0],
        "reference": args.reference[0],
        "separator": args.separator[0],
        "prefix_len": args.prefix_len[0],
        "suffix_len": args.suffix_len[0],
        "debug": args.debug[0],
        "only_homopolymer": args.only_homopolymers[0],
        "minimum_support_reads": args.minimum_support_reads[0],
        "threads": args.threads[0],
        "batch": args.batch[0],
        "only_microsatellites": args.only_microsatellites[0],
        "ranges_of_repeat_times": {}
    }

    for i in args.minimum_repeat_times[0].split(";"):
        unit_range, repeat_range = i.split(":")
        if "-" in unit_range:
            unit_start, unit_end = tuple(map(int, unit_range.split("-")))
        else:
            unit_start = int(unit_range)
            unit_end = unit_start
        repeat_start = int(repeat_range)
        # print(unit_start,unit_end,"  ",repeat_start, repeatEnd)
        for ur in range(unit_start, unit_end + 1):
            if ur not in paras["ranges_of_repeat_times"]:
                paras["ranges_of_repeat_times"][ur] = {}
            paras["ranges_of_repeat_times"][ur]["min"] = repeat_start
        for i in args.maximum_repeat_times[0].split(";"):
            # print(i)
            unit_range, repeat_range = i.split(":")
            if "-" in unit_range:
                unit_start, unit_end = tuple(map(int, unit_range.split("-")))
            else:
                unit_start = int(unit_range)
                unit_end = unit_start
            repeat_start = int(repeat_range)
            # print(unit_start,unit_end,"  ",repeat_start, repeatEnd)
            for ur in range(unit_start, unit_end + 1):
                if ur not in paras["ranges_of_repeat_times"]:
                    paras["ranges_of_repeat_times"][ur] = {}
                paras["ranges_of_repeat_times"][ur]["max"] = repeat_start
    error_stat = False
    if os.path.exists(paras["input"]):
        logger.info("The input file is : '" + paras["input"] + "'.")
    else:
        logger.error('The input file '
                     + paras["input"] + ' is not exist, please check again')
        error_stat = True

    if os.path.isfile(paras["microsatellite"]):
        logger.info("The microsatellites file  is : " + paras["microsatellite"])
    else:
        logger.error('The microsatellites file '
                     + paras["microsatellite"] + ' is not exist, please check again')
        error_stat = True
    if os.path.isfile(paras["reference"]):
        logger.info("The reference file is : '" + paras["reference"] + "'.")
    else:
        paras["reference"] = "" if paras["reference"] == "." else paras["reference"]
        logger.error('The reference file ' + paras["reference"] + ' is not exist, please check again')
        error_stat = True
    if paras["input"][-4:] == "cram":
        paras["input_format"] = "cram"
        cramfile = pysam.AlignmentFile(paras["input"], mode="rb", reference_filename=paras["reference"])
        if not cramfile.has_index():
            logger.info("Build index for the input cram ...")
            pysam.index(paras["input"])
        cramfile.close()
    else:
        paras["input_format"] = "hap1"
        bam_file = pysam.AlignmentFile(paras["input"], mode="rb")
        if not bam_file.has_index():
            logger.info("Build index for the input bam ...")
            pysam.index(paras["input"])
        bam_file.close()
    paras["output_vcf"] = paras["output"] + ".vcf.gz"
    if not os.path.exists(paras["output_vcf"]):
        logger.info("The output is : " + paras["output_vcf"] + ".")
    else:
        logger.error(
            '[ERROR] The output ' + paras["output_vcf"] +
            ' is still exist! in case of overwrite files in this workspace, '
            'please check your script!')
        if not paras["debug"]:
            error_stat = True
    if error_stat: return False
    set_value("paras", paras)
    return True


def bm_process_one_ms_site(msDetail):
    msDetail.get_dis()
    msDetail.get_pileup_info()
    return msDetail


def bm_write_vcf_init(outputpath):
    outputfile = pysam.VariantFile(outputpath, "w")
    bam_file = pysam.AlignmentFile(get_value("paras")["input"], "rb")
    contigs = bam_file.references
    contigsLen = bam_file.lengths
    chromList = get_value("chrom_list")
    contigs_len_dict = {}
    sortedContig = []
    for contig, length in zip(contigs, contigsLen):
        if contig in chromList:
            sortedContig.append(contig)
            contigs_len_dict[contig] = length
    for contig in sortedContig:
        outputfile.header.add_line(
            "##contig=<ID={chrom},length={length}>".format(chrom=contig, length=contigs_len_dict[contig]))
    set_value("contigs_info", contigs_len_dict)
    outputfile.header.add_line('##INFO=<ID=chrom,Number=1,Type=String,Description="Chromosome">')
    outputfile.header.add_line('##INFO=<ID=pos,Number=1,Type=Integer,Description="Position">')
    outputfile.header.add_line('##INFO=<ID=ms_start,Number=1,Type=Integer,Description='
                               '"Position start of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=ms_end,Number=1,Type=Integer,Description='
                               '"Position end of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=motif,Number=1,Type=String,Description="Repeat unit">')
    outputfile.header.add_line('##INFO=<ID=repeat_times,Number=1,Type=Integer,Description='
                               '"Repeat times of motif in reference">')
    outputfile.header.add_line('##INFO=<ID=motif_len,Number=1,Type=Integer,Description="Repeat unit length">')
    outputfile.header.add_line('##INFO=<ID=ref_repeat_length,Number=1,Type=Integer,Description='
                               '"length of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=start_pre,Number=1,Type=Integer,Description="Start position for analysis">')
    outputfile.header.add_line('##INFO=<ID=end_suf,Number=1,Type=Integer,Description="End position for analysis">')
    outputfile.header.add_line('##INFO=<ID=mut_start,Number=1,Type=Integer,Description="Start position of mutaiton">')
    outputfile.header.add_line('##INFO=<ID=mut_end,Number=1,Type=Integer,Description="End position of mution">')
    outputfile.header.add_line('##INFO=<ID=query_repeat_length,Number=1,Type=Integer,Description='
                               '"Evaluation repeat length of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=dis_stat,Number=1,Type=String,Description='
                               '"True,the distribution is available">')
    outputfile.header.add_line('##INFO=<ID=check,Number=1,Type=String,Description='
                               '"True,the site is available for benchmark">')
    outputfile.header.add_line('##INFO=<ID=check_stats,Number=1,Type=String,Description='
                               '"Why this site is not available for benchmark">')
    outputfile.header.add_line('##INFO=<ID=allele,Number=1,Type=Integer,Description="Allele number in this site">')
    outputfile.header.add_line('##INFO=<ID=dis,Number=1,Type=String,Description='
                               'Distribution of repeat length>')
    outputfile.header.add_line('##INFO=<ID=var_type,Number=1,Type=String,Description='
                               '"Variant typeComplex, SNP, DEL, INS">')
    outputfile.header.add_line('##INFO=<ID=var_type_list,Number=1,Type=String,Description='
                               '"Variant type, Complex, SNP, DEL, INS">')
    outputfile.header.add_line('##INFO=<ID=var_detail,Number=1,Type=String,Description='
                               '"Variant Details, rec1|rec2 rec1: pos:type:content">')
    return outputfile


def bm_write_vcf_close(outputfile):
    outputfile.close()
    pysam.tabix_index(get_value("paras")["output_vcf"], force=True, preset="vcf")


def benchmark(parase):
    if not benchmark_init(parase):
        logger.error("Benchmark init ERROR!")
        return -1
        # return if the process the arguments errors
    args = get_value("paras")
    out_vcf = args["output_vcf"]
    output_file = bm_write_vcf_init(out_vcf)
    contigs_info = get_value("contigs_info")
    df_microsatellites = load_microsatellites(args)
    for contig, contig_len in contigs_info.items():
        logger.info("\t--------------------------------------------------------------------------------")
        logger.info("\tProcessing " + contig + "...")
        this_contig_microsatellite = df_microsatellites[df_microsatellites["chr"] == contig].sort_values("pos")
        window_ms = []
        ms_num = 0
        for ms_id, info in this_contig_microsatellite.iterrows():
            ms_num += 1
            info["prefix_len"] = args["prefix_len"]
            info["suffix_len"] = args["suffix_len"]
            info["reference"] = args["reference"]
            window_ms.append(info)
            if ms_num % (args["batch"] * args["threads"]) == 0:
                window = Window(contig, window_ms)
                window.run_window(output_file)
                window_ms = []
        if len(window_ms) > 0:
            window = Window(contig, window_ms)
            del window_ms
            window.run_window(output_file)
    bm_write_vcf_close(output_file)
