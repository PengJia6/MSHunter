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
from src.units import *
from src.Microsatellite import Microsatellite
from multiprocess.pool import Pool


def process_one_ms(info):
    ms = Microsatellite(info)
    dis = ms.get_reads_info()
    if len(dis) > 0:
        dis_str = "|".join([":".join([str(a), str(b)]) for a, b in dis.items()])
        dis_mean, dis_std = dis_stats(dis, ["mean", 'std'])
        return ms.ms_id, ms.repeat_unit, ms.repeat_times, dis_mean, dis_std, dis_str
    else:
        return None, None, None, None


def pre_stat(paras, df_microsatellites):
    # reference=paras["reference"]
    path_pre_stat = paras["output_tmp"].rstrip("/") + "/" + get_value("case") + ".stat"
    file_all_stat = open(path_pre_stat, "w")
    df_microsatellites_download_sample = microsatellites_sampling(df_microsatellites, paras)
    for repeat_unit, info in df_microsatellites_download_sample.items():
        for repeat_times, ms_infos in info.items():
            logger.info("Processing   repeat unit: " + str(repeat_unit) + " repeat times: " + str(repeat_times))
            infos = []
            for id, info in ms_infos.iterrows():
                info["reference"] = paras["reference"]
                info["prefix_len"] = paras["prefix_len"]
                info["suffix_len"] = paras["suffix_len"]
                infos.append(info)
            pool = Pool(processes=paras["threads"])
            res_infos = pool.map(process_one_ms, infos)
            pool.close()
            pool.join()
            suffix_str = "." + str(repeat_unit) + "." + str(repeat_times)
            file = open(path_pre_stat + suffix_str, "w")
            this_repeat_means = []
            this_repeat_stds = []
            num = 0
            for res in res_infos:
                if None not in res:
                    file.write("\t".join(list(map(str, res))))
                    this_repeat_means.append(res[3])
                    this_repeat_stds.append(res[4])

            file.close()
            this_repeat_mean_mean = np.mean(this_repeat_means)
            this_repeat_mean_std = np.std(this_repeat_means)
            this_repeat_std_mean = np.mean(this_repeat_stds)
            this_repeat_std_std = np.std(this_repeat_stds)
            this_info_list = list(map(str, [repeat_unit, repeat_times,
                                            this_repeat_mean_mean, this_repeat_mean_std,
                                            this_repeat_std_mean, this_repeat_std_std,
                                            ]))
            file_all_stat.write("\t".join(this_info_list) + "\n")
    file_all_stat.close()
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
