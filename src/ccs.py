#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : ccs.py
# Author : Peng Jia
# Date   : 2020.08.14
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
from src.units import *
from src.pre_stat import pre_stat
from src.estimate_error import estimate_error
from src.call_variants import call_variants





def genotype_ccs(paras):
    df_microsatellites = load_microsatellites(paras)
    pre_stat(df_microsatellites)
    estimate_error()
    call_variants(df_microsatellites)

    # pre_stat(paras, df_microsatellites)

    # print(df_microsatellites_download_sample)
    # print("dklldfldl",paras)
    pass
