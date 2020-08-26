#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : pre_stat.py
# Author : Peng Jia
# Date   : 2020.08.20
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO TODO
=============================================================================="""
import pandas as pd
from src.global_dict import *
from src.Microsatellite import Microsatellite


def pre_stat(paras, df_microsatellites):
    # reference=paras["reference"]
    pre_stat = paras["output_tmp"] + get_value("case") + ".stat"
    df_microsatellites_download_sample = microsatellites_sampling(df_microsatellites, paras)
    for repeat_unit, info in df_microsatellites_download_sample.items():
        for repeat_times, ms_infos in info.items():
            num = 0
            for id, info in ms_infos.iterrows():
                ms_id = info["chr"] + "_" + str(int(info["pos"]))
                num += 1
                info["reference"] = paras["reference"]
                info["prefix_len"] = paras["prefix_len"]
                info["suffix_len"] = paras["suffix_len"]
                ms = Microsatellite(info)
                dis = ms.get_reads_info()
                print(ms_id, repeat_unit, repeat_times, dis)
                # print(ms, info)
            # print(ms_infos)

    # print(df_microsatellites_download_sample)

    return


def microsatellites_sampling(df_microsatellite, paras, sample_num=1000):
    df_microsatellite_downsample = {}
    repeat_ranges = paras["ranges_of_repeat_times"]
    for repeat_unit_length in sorted(repeat_ranges.keys()):
        df_microsatellite_repeat_unit = df_microsatellite[df_microsatellite["motifLen"] == repeat_unit_length]
        if len(df_microsatellite_repeat_unit) > 0 and repeat_unit_length not in df_microsatellite_downsample:
            df_microsatellite_downsample[repeat_unit_length] = {}
        for repeat_times in range(repeat_ranges[repeat_unit_length]["min"], repeat_ranges[repeat_unit_length]["max"]):
            df_microsatellite_repeat_times = df_microsatellite_repeat_unit[
                df_microsatellite_repeat_unit["repeatTimes"] == repeat_times]
            if len(df_microsatellite_repeat_times) > sample_num:
                df_microsatellite_repeat_times = df_microsatellite_repeat_times.sample(sample_num)
            if len(df_microsatellite_repeat_times) > 0:
                df_microsatellite_downsample[repeat_unit_length][repeat_times] = df_microsatellite_repeat_times
    return df_microsatellite_downsample
