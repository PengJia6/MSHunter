#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : global_dict.py
# Author : Peng Jia
# Date   : 2020.07.13
# Email  : pengjia@stu.xjtu.edu.cn
# Description: Define global variables
=============================================================================="""

def global_init():
    global _global_dict
    _global_dict = {}
    _global_dict["ms_number"] = 0
    _global_dict["tools_version"] = "1.0"
    _global_dict["tools_name"] = "mstools"
    _global_dict["chrom_list"] = [str(i) for i in range(1, 23)] + \
                                 ["chr" + str(i) for i in range(1, 23)] + \
                                 ["X", "Y", "chrX", "chrY", "chrM", "MT"]
    _global_dict["default"] = {
        "bam2dis": {
            "threads": 4,
            "minimum_mapping_quality": 1,
            "minimum_support_reads": 2,
            "batch": 2000,
            "debug": False,
            "separator": "tab",  # comma,tab,space
            "only_homopolymers": False,
            "allow_mismatch": True,
            "minimum_repeat_times": "1:8;2-5:5",
            "maximum_repeat_times": "1-5:100",

        },
        "errEval": {
            "threads": 4,
            "batch": 2000,
        },

        "call": {
        },
        "genotype": {
            "reference": ".",
            "threads": 4,
            "minimum_mapping_quality": 1,
            "minimum_support_reads": 2,
            "batch": 2000,
            "debug": False,
            "separator": "tab",  # comma,tab,space
            "only_homopolymers": False,
            "allow_mismatch": True,
            "minimum_repeat_times": "1:8;2-5:5",
            "maximum_repeat_times": "1-5:100",
            "prefix_len": 500,
            "suffix_len": 500,
            "kmer_size": 5,
            "minimum_phasing_reads": 3,
            "tech": "ccs",
            "hap": True,
        },
        "benchmark": {
            "reference": ".",
            "threads": 4,
            "minimum_mapping_quality": 1,
            "minimum_support_reads": 2,
            "batch": 2000,
            "debug": False,
            "separator": "tab",  # comma,tab,space
            "only_homopolymers": False,
            "allow_mismatch": True,
            "minimum_repeat_times": "1:8;2-5:5",
            "maximum_repeat_times": "1-5:100",
            "prefix_len": 5,
            "suffix_len": 5,
            "kmer_size": 5,
            "minimum_phasing_reads": 3,
            "tech": "ccs",
            "hap": True,
            "only_microsatellites":True
        },

    }


def set_value(name, value):
    _global_dict[name] = value


def get_value(name, defValue=None):
    try:
        return _global_dict[name]
    except KeyError:
        print("[ERROR] No variable", name, "in global_dict")
        return defValue
