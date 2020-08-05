#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : Microsatellite.py
# Author : Peng Jia
# Date   : 2020.08.04
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
from src.global_dict import *


class Microsatellite:

    def __init__(self, ms_info,
                 # prefix_len=get_value("paras")["prefix_len"],
                 # suffix_len=get_value("paras")["suffix_len"]
                 ):
        self.chrom = ms_info["chr"]
        self.start = ms_info["pos"]
        self.ms_id = self.chrom + str(self.start)
        self.repeat_times = ms_info["repeatTimes"]
        self.repeat_unit = ms_info["motif"]
        self.repeat_unit_len = ms_info["motifLen"]
        self.repeat_len = self.repeat_times * self.repeat_unit_len
        self.end = self.start + self.repeat_len
        self.start_pre = self.start - ms_info["prefix_len"]
        self.end_suf = self.end + ms_info["suffix_len"]
        self.reads_info = {}
        # print(ms_info)

    # def update_ms_id(self,ms_id):
    #     if ms_id not in self.reads_info
    #     pass
