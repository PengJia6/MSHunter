#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : benchmark_merge.py
# Author : Peng Jia
# Date   : 2020.07.16
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
from src.global_dict import *
import os


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
        print("[INFO] The haplotype 1 microsatellite calling result is : '" + paras["input"] + "'.")
    else:
        print('[ERROR] The aplotype 1 microsatellite calling result '
              + paras["input"] + ' is not exist, please check again')
        error_stat = True
    if os.path.exists(paras["hap2"]):
        print("[INFO] The haplotype 2 microsatellite calling result is : '" + paras["input"] + "'.")
    else:
        print('[ERROR] The aplotype 2 microsatellite calling result '
              + paras["input"] + ' is not exist, please check again')
        error_stat = True

    set_value("paras", paras)
    if error_stat: return False
    return True


def benchmark_merge(parase):
    if not benchmark_merge_init(parase):
        return -1
