#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : genotype.py
# Author : Peng Jia
# Date   : 2020.07.13
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
import os
from src.bam2dis import *
from src.call import *
from src.errEval import *

def genotype_init(args):
    """
    argument procress
    """
    paras = {}
    paras["input"] = args.input[0]
    paras["output"] = args.output[0]
    paras["microsatellite"] = args.microsatellite[0]
    paras["reference"] = args.reference[0]
    paras["separator"] = args.separator[0]
    paras["tech"] = args.technology[0]
    paras["hap"] = args.haplotype_bam[0]
    paras["prefix_len"] = args.prefix_len[0]
    paras["suffix_len"] = args.suffix_len[0]
    paras["debug"] = args.debug[0]
    paras["only_homopolymer"] = args.only_homopolymers[0]
    paras["minimum_support_reads"] = args.minimum_support_reads[0]
    paras["minimum_mapping_quality"] = args.minimum_mapping_quality[0]
    paras["allow_mismatch"] = args.allow_mismatch[0]
    paras["threads"] = args.threads[0]
    paras["batch"] = args.batch[0]
    paras["minimum_phasing_reads"] = args.minimum_phasing_reads[0]
    paras["ranges_of_repeat_times"] = {}

    for i in args.minimum_repeat_times[0].split(";"):
        unitRange, repeatRange = i.split(":")
        if "-" in unitRange:
            unitStart, unitEnd = tuple(map(int, unitRange.split("-")))
        else:
            unitStart = int(unitRange)
            unitEnd = unitStart
        repeatStart = int(repeatRange)
        # print(unitStart,unitEnd,"  ",repeatStart, repeatEnd)
        for ur in range(unitStart, unitEnd + 1):
            if ur not in paras["ranges_of_repeat_times"]:
                paras["ranges_of_repeat_times"][ur] = {}
            paras["ranges_of_repeat_times"][ur]["min"] = repeatStart
        for i in args.maximum_repeat_times[0].split(";"):
            # print(i)
            unitRange, repeatRange = i.split(":")
            if "-" in unitRange:
                unitStart, unitEnd = tuple(map(int, unitRange.split("-")))
            else:
                unitStart = int(unitRange)
                unitEnd = unitStart
            repeatStart = int(repeatRange)
            # print(unitStart,unitEnd,"  ",repeatStart, repeatEnd)
            for ur in range(unitStart, unitEnd + 1):
                if ur not in paras["ranges_of_repeat_times"]:
                    paras["ranges_of_repeat_times"][ur] = {}
                paras["ranges_of_repeat_times"][ur]["max"] = repeatStart
    error_stat = False
    if os.path.exists(paras["input"]):
        print("[INFO] The input file is : '" + paras["input"] + "'.")
    else:
        print('[ERROR] The input file '
              + paras["input"] + ' is not exist, please check again')
        error_stat = True
    if os.path.isfile(paras["microsatellite"]):
        print("[INFO] The microsatellites file  is : " + paras["microsatellite"])
    else:
        print('[ERROR] The microsatellites file '
              + paras["microsatellite"] + ' is not exist, please check again')
        error_stat = True
    if os.path.isfile(paras["reference"]):
        print("[INFO] The reference file is : '" + paras["reference"] + "'.")
    else:
        paras["reference"] = "" if paras["reference"] == "." else paras["reference"]
        print('[ERROR] The reference file ' + paras["reference"] + ' is not exist, please check again')
        error_stat = True
    if paras["input"][-4:] == "cram":
        paras["input_format"] = "cram"
        cramfile = pysam.AlignmentFile(paras["input"], mode="rb", reference_filename=paras["reference"])
        if not cramfile.has_index():
            print("[INFO] Build index for the input cram ...")
            pysam.index(paras["input"])
        cramfile.close()
    else:
        paras["input_format"] = "bam"
        bamfile = pysam.AlignmentFile(paras["input"], mode="rb")
        if not bamfile.has_index():
            print("[INFO] Build index for the input bam ...")
            pysam.index(paras["input"])
        bamfile.close()

    if not os.path.exists(paras["output"]):
        print("[INFO] The output is : " + paras["output"] + ".")
    else:
        print(
            '[ERROR] The output ' + paras["output"] +
            ' is still exist! in case of overwrite files in this workspace, '
            'please check your script!')
        if not paras["debug"]:
            error_stat = True
    if error_stat: return False
    output_path = paras["output"]
    output_path = output_path if output_path[-1] == "/" else output_path + "/"
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    paras["output"] = output_path
    input_path = paras["input"]
    input_path = input_path[:-1] if input_path[-1] == "/" else input_path
    case = input_path.split("/")[-1].strip(".bam")
    case = case.strip(".cram")
    paras["output_dis"] = paras["output"] + case + ".dis.vcf.gz"
    paras["output_tmp"] = paras["output"] + case + "_tmp"
    if not os.path.exists(paras["output_tmp"]):
        os.makedirs(paras["output_tmp"])
    paras["output_model"] = paras["output"] + case + ".model"
    paras["output_call"] = paras["output"] + case + ".vcf.gz"
    set_value("case", case)
    set_value("paras", paras)
    return True


def genotype(parase):
    if genotype_init(parase):
        bam2dis()
        errEval()
        call()
