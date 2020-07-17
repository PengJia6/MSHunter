#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : benchmark_merge.py
# Author : Peng Jia
# Date   : 2020.07.16
# Email  : pengjia@stu.xjtu.edu.cn
# Description: Merge two haplotype calling result
=============================================================================="""
from src.global_dict import *
import os
import pysam


def benchmark_merge_init(args):
    """
        argument procress
        """
    paras = {}
    paras["hap1"] = args.hap1[0]
    paras["hap2"] = args.hap2[0]
    paras["output"] = args.output[0]
    paras["reference"] = args.reference[0]
    error_stat = False
    if os.path.exists(paras["hap1"]):
        print("[INFO] The haplotype 1 microsatellite calling result is : '" + paras["hap1"] + "'.")
    else:
        print('[ERROR] The aplotype 1 microsatellite calling result '
              + paras["input"] + ' is not exist, please check again')
        error_stat = True
    if os.path.exists(paras["hap2"]):
        print("[INFO] The haplotype 2 microsatellite calling result is : '" + paras["hap2"] + "'.")
    else:
        print('[ERROR] The aplotype 2 microsatellite calling result '
              + paras["input"] + ' is not exist, please check again')
        error_stat = True

    set_value("paras", paras)
    if error_stat: return False
    return True


def merge():
    paras = get_value("paras")
    hap1 = pysam.VariantFile(paras["hap1"])
    hap2 = pysam.VariantFile(paras["hap2"])
    mergerd_file = pysam.VariantFile(paras["output"], "w", header=hap1.header)
    # outputfile.header.add_line('##INFO=<ID=MutationType,Number=1,Type=String,Description="Mutation type">')
    mergerd_file.header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    # outputfile.header.add_line('##FORMAT=<ID=DP,Number=1,Type=String,Description="Allele Depth">')

    for hap1_rec, hap2_rec in zip(hap1.fetch(), hap2.fetch()):
        merged_rec = mergerd_file.new_record()
        merged_rec.contig = hap1_rec.contig
        merged_rec.id = hap1_rec.id
        mut_start = min(hap1_rec.pos, hap2_rec.pos)
        mut_end=max(hap1_rec.info["mut_end"],hap2_rec.info["mut_end"])
        if hap2_rec.info["mut_start"]!=hap2_rec.info["mut_start"]:
            print("kdkdk")

        ref_prefix=""
        hap1_prefix=""
        hap2_prefix=""
        ref_suffix = ""
        hap1_suffix = ""
        hap2_suffix = ""
        if hap1_rec.info["mut_start"]>mut_start:
            ref_prefix=hap2_rec.ref[0:hap1_rec.info["mut_start"]-mut_start]
            hap1_prefix=hap2_rec.ref[0:hap1_rec.info["mut_start"]-mut_start]
        if hap2_rec.info["mut_start"] > mut_start:
            # ref_prefix = hap1_rec.ref[0:hap2_rec.info["mut_start"] - mut_start]
            hap2_prefix = hap1_rec.ref[0:hap2_rec.info["mut_start"] - mut_start]
        if hap1_rec.info["mut_end"]< mut_end:
            ref_suffix=hap2_rec.ref[-(hap1_rec.info["mut_end"]-mut_end):]
            hap1_suffix=hap2_rec.ref[-(hap1_rec.info["mut_end"]-mut_end):]
        if hap2_rec.info["mut_end"]< mut_end:
            # ref_suffix=hap1_rec.ref[-(hap2_rec.info["mut_end"]-mut_end):]
            hap2_suffix=hap1_rec.ref[-(hap2_rec.info["mut_end"]-mut_end):]
        merged_rec.ref=ref_prefix+hap1_rec.ref+ref_suffix
        hap1_alt=hap1_prefix+hap1_rec.alts[0]+hap1_suffix
        hap2_alt=hap2_prefix+hap2_rec.alts[0]+hap2_suffix
        if hap1_alt!=hap2_alt:
            print("hap1",hap1_rec)
            print("hap2",hap2_rec)
        # merged_rec.alts=




        # merged_rec.pos=mut_start
        # # merged_rec=hap2_rec
        # # merged_rec.info["ms_start"]=0
        # print(merged_rec)
        #
        # print(hap2_rec)
        # print("i", i)
        # print("j", j)


def benchmark_merge(parase):
    if not benchmark_merge_init(parase):
        return -1
    merge()
