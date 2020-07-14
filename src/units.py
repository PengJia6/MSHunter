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


if __name__ == "__main__":
    ""
    dict1 = {1: 1, 2: 2, 3: 3}
    dict2 = {1: 2, 2: 4, 3: 8}
    start = time.time()
    step = 90000
    for i in range(step):
        getDisdistance(dict2, dict1)
    end = time.time()
    print(end - start)
    start = time.time()
    for i in range(step):
        getDisdistance2(dict2, dict1)
    end = time.time()
    print(end - start)
    start = time.time()
    for i in range(step):
        getDisdistance3(dict2, dict1)
    end = time.time()
    print(end - start)
    print(getDisdistance(dict1, dict2))
    print(getDisdistance2(dict1, dict2))
    print(getDisdistance3(dict1, dict2))
    # getDisdistance(dict2,dict1)


def load_microsatellites(args):
    """
    Args:
        args ():

    Returns:
    """
    print("[INFO] Loading microsatellite file...")
    ms = args["microsatellite"]
    separator = args["separator"]
    # ID,chr,pos,motif,motifLen,repeatTimes,prefix,suffix
    if separator == "comma":
        dfMicroSatellites = pd.read_csv(ms, index_col=0)
        # print(dfMicroSatellites.columns)
    elif separator == "tab":
        dfMicroSatellites = pd.read_table(ms, header=0)
        columns = dfMicroSatellites.columns
        if "chromosome" in columns:
            dfMicroSatellites.rename(columns={"chromosome": "chr"}, inplace=True)
            # dfMicroSatellites["chr"] = dfMicroSatellites["chromosome"]
            # del dfMicroSatellites["chromosome"]
        if "location" in columns:
            dfMicroSatellites.rename(columns={"location": "pos"}, inplace=True)
            # dfMicroSatellites["pos"] = dfMicroSatellites["location"]
            # del dfMicroSatellites["location"]
        if "repeat_unit_bases" in columns:
            dfMicroSatellites.rename(columns={"repeat_unit_bases": "motif"}, inplace=True)
            # dfMicroSatellites["motif"] = dfMicroSatellites["repeat_unit_bases"]
            # del dfMicroSatellites["repeat_unit_bases"]
        if "repeat_unit_length" in columns:
            dfMicroSatellites.rename(columns={"repeat_unit_length": "motifLen"}, inplace=True)
            # dfMicroSatellites["motifLen"] = dfMicroSatellites["repeat_unit_length"]
            # del dfMicroSatellites["repeat_unit_length"]
        if "repeat_times" in columns:
            dfMicroSatellites.rename(columns={"repeat_times": "repeatTimes"}, inplace=True)
            # dfMicroSatellites["repeatTimes"] = dfMicroSatellites["repeat_times"]
            # del dfMicroSatellites["repeat_times"]
        if "left_flank_bases" in columns:
            dfMicroSatellites.rename(columns={"left_flank_bases": "prefix"}, inplace=True)
            # dfMicroSatellites["prefix"] = dfMicroSatellites["left_flank_bases"]
            # del dfMicroSatellites["left_flank_bases"]
        if "right_flank_bases" in columns:
            dfMicroSatellites.rename(columns={"right_flank_bases": "suffix"}, inplace=True)
            # dfMicroSatellites["suffix"] = dfMicroSatellites["right_flank_bases"]
            # del dfMicroSatellites["right_flank_bases"]
        # dfMicroSatellites.index = dfMicroSatellites["chr"] + "_" + dfMicroSatellites["pos"].astype(str)
    elif separator == "space":
        dfMicroSatellites = pd.read_table(ms, header=0, sep=" ")
    chromList = get_value("chrom_list")
    dfMicroSatellites = dfMicroSatellites[dfMicroSatellites['chr'].isin(chromList)]
    if args["only_homopolymer"]:
        dfMicroSatellites = dfMicroSatellites[dfMicroSatellites['motifLen'] == 1]

    repeatRange = args["ranges_of_repeat_times"]
    repeatUnitList = sorted(repeatRange.keys())

    newDf = pd.DataFrame()
    for ul in repeatUnitList:
        minr = repeatRange[ul]["min"]
        maxr = repeatRange[ul]["max"]
        newDf = pd.concat([newDf, dfMicroSatellites[(dfMicroSatellites["motifLen"] == ul) &
                                                    (dfMicroSatellites["repeatTimes"] >= minr) &
                                                    (dfMicroSatellites["repeatTimes"] <= maxr)
                                                    ]])
    if args["debug"]:
        locis_num = 80000
        newDf = newDf.iloc[100:locis_num + 100, :]
        # dfMicroSatellites = dfMicroSatellites[dfMicroSatellites["motifLen"] == 1]
        # if len(newDf) > locis_num:
        #     newDf = newDf.sample(locis_num)
    newDf.sort_index(inplace=True)

    print("[INFO] There are total", len(newDf), "microsatellites.")
    set_value("ms_number", len(newDf))
    set_value("motifList", set(newDf["motif"]))
    # print(set(newDf["motif"]))
    return newDf


def get_max_support_index(input_dict):
    """
    @param input_dict: dict value with float or int type
    @type input_dict: dict
    @return: key of the max value
    @rtype: float / int
    """
    m = max(input_dict.keys(), key=(lambda x: input_dict[x]))
    return m
if __name__ =="__main__":
    a=get_max_support_index({1:5,6:40,3:2})
    a=get_max_support_index({})
    print(a)