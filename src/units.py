#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : units.py
# Author : Peng Jia
# Date   : 2020.07.13
# Email  : pengjia@stu.xjtu.edu.cn
# Description: Useful functions
=============================================================================="""
import numpy as np
import time
import pandas as pd
from src.global_dict import *


# TODO check and normalize
def removeZeroDict(dict):
    newdict = {}
    for key, value in dict.items():
        if value > 0.000001:
            newdict[key] = value
    return newdict


# the fastest
def getDisdistance(dict1, dict2):
    dictkey = list(dict1.keys()) + list(dict2.keys())
    for key in dictkey:
        if key not in dict1:
            dict1[key] = 0
        if key not in dict2:
            dict2[key] = 0
    sum = 0
    for key in dictkey:
        err = dict1[key] - dict2[key]
        sum += (err * err)
    return round(np.sqrt(sum), 6)


# TODO check and normalize
def getDisdistance2(dict1, dict2):
    dictkey = list(dict1.keys()) + list(dict2.keys())
    list1 = []
    list2 = []
    for key in dictkey:
        if key not in dict1:
            list1.append(0)
        else:
            list1.append(dict1[key])
        if key not in dict2:
            list2.append(0)
        else:
            list2.append(dict2[key])
    # print(np.linalg.norm(np.array(list1)-np.array(list2)))
    # print(np.sqrt(np.sum(np.square(np.array(list1)-np.array(list2)))))
    return np.linalg.norm(np.array(list1) - np.array(list2))


# TODO check and normalize

def getDisdistance3(dict1, dict2):
    dictkey = list(dict1.keys()) + list(dict2.keys())
    list1 = []
    list2 = []
    for key in dictkey:
        if key not in dict1:
            list1.append(0)
        else:
            list1.append(dict1[key])
        if key not in dict2:
            list2.append(0)
        else:
            list2.append(dict2[key])
    # print(np.linalg.norm(np.array(list1)-np.array(list2)))
    return np.sqrt(np.sum(np.square(np.array(list1) - np.array(list2))))


def load_microsatellites(args):
    """
    Description: load Microsatellite information
    Stat: PASS
    """
    logger.info("Loading microsatellite file...")
    ms = args["microsatellite"]
    separator = args["separator"]
    # ID,chr,pos,motif,motifLen,repeatTimes,prefix,suffix
    if separator == "comma":
        df_microSatellites = pd.read_csv(ms, index_col=0)
        # print(df_microSatellites.columns)
    elif separator == "tab":
        df_microSatellites = pd.read_table(ms, header=0)
        columns = df_microSatellites.columns
        if "chromosome" in columns:
            df_microSatellites.rename(columns={"chromosome": "chr"}, inplace=True)
            # df_microSatellites["chr"] = df_microSatellites["chromosome"]
            # del df_microSatellites["chromosome"]
        if "location" in columns:
            df_microSatellites.rename(columns={"location": "pos"}, inplace=True)
            # df_microSatellites["pos"] = df_microSatellites["location"]
            # del df_microSatellites["location"]
        if "repeat_unit_bases" in columns:
            df_microSatellites.rename(columns={"repeat_unit_bases": "motif"}, inplace=True)
            # df_microSatellites["motif"] = df_microSatellites["repeat_unit_bases"]
            # del df_microSatellites["repeat_unit_bases"]
        if "repeat_unit_length" in columns:
            df_microSatellites.rename(columns={"repeat_unit_length": "motifLen"}, inplace=True)
            # df_microSatellites["motifLen"] = df_microSatellites["repeat_unit_length"]
            # del df_microSatellites["repeat_unit_length"]
        if "repeat_times" in columns:
            df_microSatellites.rename(columns={"repeat_times": "repeatTimes"}, inplace=True)
            # df_microSatellites["repeatTimes"] = df_microSatellites["repeat_times"]
            # del df_microSatellites["repeat_times"]
        if "left_flank_bases" in columns:
            df_microSatellites.rename(columns={"left_flank_bases": "prefix"}, inplace=True)
            # df_microSatellites["prefix"] = df_microSatellites["left_flank_bases"]
            # del df_microSatellites["left_flank_bases"]
        if "right_flank_bases" in columns:
            df_microSatellites.rename(columns={"right_flank_bases": "suffix"}, inplace=True)
            # df_microSatellites["suffix"] = df_microSatellites["right_flank_bases"]
            # del df_microSatellites["right_flank_bases"]
        # df_microSatellites.index = df_microSatellites["chr"] + "_" + df_microSatellites["pos"].astype(str)
    elif separator == "space":
        df_microSatellites = pd.read_table(ms, header=0, sep=" ")
    chromList = get_value("chrom_list")
    df_microSatellites = df_microSatellites[df_microSatellites['chr'].isin(chromList)]
    if args["only_homopolymer"]:
        df_microSatellites = df_microSatellites[df_microSatellites['motifLen'] == 1]

    repeatRange = args["ranges_of_repeat_times"]
    repeatUnitList = sorted(repeatRange.keys())

    df_microsatellite_pass = pd.DataFrame()
    for ul in repeatUnitList:
        minr = repeatRange[ul]["min"]
        maxr = repeatRange[ul]["max"]
        df_microsatellite_pass = pd.concat(
            [df_microsatellite_pass, df_microSatellites[(df_microSatellites["motifLen"] == ul) &
                                                        (df_microSatellites["repeatTimes"] >= minr) &
                                                        (df_microSatellites["repeatTimes"] <= maxr)
                                                        ]])
    # if args["debug"]:
    #     locis_num = 10000
    #     # df_microsatellite_pass=df_microsatellite_pass[df_microSatellites["pos"]>143200000]
    #     df_microsatellite_pass = df_microsatellite_pass.iloc[100000:locis_num + 100000, :]
    #     # df_microSatellites = df_microSatellites[df_microSatellites["motifLen"] == 1]
    #     # if len(df_microsatellite_pass) > locis_num:
    #     #     df_microsatellite_pass = df_microsatellite_pass.sample(locis_num)
    # df_microsatellite_pass.sort_index(inplace=True)

    logger.info("There are total " + str(len(df_microsatellite_pass)) + " microsatellites.")
    set_value("ms_number", len(df_microsatellite_pass))
    set_value("motifList", set(df_microsatellite_pass["motif"]))
    return df_microsatellite_pass


def get_max_support_index(input_dict):
    """
    @param input_dict: dict value with float or int type
    @type input_dict: dict
    @return: key of the max value
    @rtype: float / int
    """
    m = max(input_dict.keys(), key=(lambda x: input_dict[x]))
    return m


def dis_stats(dis, values=["mean", "std"]):
    repeat_length_list = []
    for repeat_len, times in dis.items():
        repeat_length_list.extend([repeat_len] * times)
    res_list = []
    for value in values:
        if value == "mean":
            this_value = np.mean(repeat_length_list)
        elif value == "std":
            this_value = np.std(repeat_length_list)
        elif value == "max":
            this_value = np.max(repeat_length_list)
        elif value == "min":
            this_value = np.min(repeat_length_list)
        else:
            this_value = None
        res_list.append(this_value)
    return res_list


def str2int(item):
    if len(item) > 0:
        return np.nanmean([ord(i) - 33 for i in item])
    else:
        return np.nan


def int2str(qual_int):
    if np.isnan(qual_int):
        return chr(33)
    else:
        return chr(int(qual_int) + 33)


# TODO chcck and normalize
def dis_sum(dict_list):
    keylist = []
    for item in dict_list:
        for key in item:
            if key not in keylist:
                keylist.append(key)
    res_dict = {}
    for key in keylist:
        res_dict[key] = 0
    for item in dict_list:
        for key in item:
            res_dict[key] += item[key]
    return res_dict


def change_dim_for_pool_map(input, threads):
    item_num = 0
    for win in input:
        item_num += len(win)
    item_per_win = item_num // threads + 1
    item_index = 0
    output = []
    win_sub = []
    for win in input:
        for item in win:
            item_index += 1
            win_sub.append(win_sub)
            if item_index % item_per_win == 0:
                output.append(win_sub)
                win_sub = []
    if len(win_sub) > 0:
        output.append(win_sub)
    return output


if __name__ == "__main__":
    a = get_max_support_index({1: 5, 6: 40, 3: 2})
    # a=get_max_support_index({})
    print(get_max_support_index({13: 8, 14: 11, 12: 1}))
#     print(a)
#     
# # TODO check and normalize
# if __name__ == "__main__":
#     ""
#     dict1 = {1: 1, 2: 2, 3: 3}
#     dict2 = {1: 2, 2: 4, 3: 8}
#     start = time.time()
#     step = 90000
#     for i in range(step):
#         getDisdistance(dict2, dict1)
#     end = time.time()
#     print(end - start)
#     start = time.time()
#     for i in range(step):
#         getDisdistance2(dict2, dict1)
#     end = time.time()
#     print(end - start)
#     start = time.time()
#     for i in range(step):
#         getDisdistance3(dict2, dict1)
#     end = time.time()
#     print(end - start)
#     print(getDisdistance(dict1, dict2))
#     print(getDisdistance2(dict1, dict2))
#     print(getDisdistance3(dict1, dict2))
