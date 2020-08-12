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
        self.ms_id = self.chrom + "_" + str(self.start)
        self.repeat_times = ms_info["repeatTimes"]
        self.repeat_unit = ms_info["motif"]
        self.repeat_unit_len = ms_info["motifLen"]
        self.repeat_len = self.repeat_times * self.repeat_unit_len
        self.end = self.start + self.repeat_len
        self.start_pre = self.start - ms_info["prefix_len"]
        self.end_suf = self.end + ms_info["suffix_len"]
        self.reads_info = {}
        self.length_dis_reads = {}
        self.depth = 0
        self.check = True
        # print(ms_info)

    def set_reads_info(self, reads_info):
        self.reads_info = reads_info
        self.depth = len(self.reads_info)

    def deletion_merge(self):
        #  = []
        all_deletions = {}

        for read_id, read_info in self.reads_info.items():
            deletions_len = len(read_info.deletions)
            deletions = []
            if deletions_len == 0:
                pass
            elif deletions_len == 1:
                deletions = [[read_info.deletions[0][0], 1]]
            else:
                current_id = read_info.deletions[0][0]
                deletion_index = {current_id: 1}
                for pos_id in range(1, deletions_len):
                    # print(read_info.deletions[pos_id][0], read_info.deletions[pos_id - 1][0])
                    if read_info.deletions[pos_id][0] == read_info.deletions[pos_id - 1][0] + 1:
                        deletion_index[current_id] += 1
                        # current_id += 1
                    else:
                        current_id = read_info.deletions[pos_id][0]
                        deletion_index[current_id] = 1
                for pos in sorted(deletion_index.keys()):
                    deletions.append([pos, deletion_index[pos]])
                # print("hhh", deletions,len(read_info.deletions))
            # all_deletions[read_id] = deletions

            self.reads_info[read_id].deletions = deletions

    # def get_dis(self):
    #     # print(self.depth)
    #     # num = 0
    #     ms_dis = {}
    #     for read_id, read_info in self.reads_info.items():
    #         if read_info.repeat_length not in ms_dis:
    #             ms_dis[read_info.repeat_length] = 1
    #         else:
    #             ms_dis[read_info.repeat_length] += 1
    #     self.ms_dis = ms_dis
    # def check_noise(self):

    def one_hap_genotype(self):
        self.deletion_merge()
        mismatches = {}
        deletions = {}
        insertions = {}
        ms_dis = {}
        for read_id, read_info in self.reads_info.items():
            for mut in read_info.mismatches:
                if mut[0] not in mismatches:
                    mismatches[mut[0]] = {}
                if mut[1] not in mismatches[mut[0]]:
                    mismatches[mut[0]][mut[1]] = 1
                else:
                    mismatches[mut[0]][mut[1]] += 1
            for mut in read_info.insertions:
                if mut[0] not in insertions:
                    insertions[mut[0]] = {}
                if mut[1] not in insertions[mut[0]]:
                    insertions[mut[0]][mut[1]] = 1
                else:
                    insertions[mut[0]][mut[1]] += 1
            for mut in read_info.deletions:
                if mut[0] not in deletions:
                    deletions[mut[0]] = {}
                if mut[1] not in deletions[mut[0]]:
                    deletions[mut[0]][mut[1]] = 1
                else:
                    deletions[mut[0]][mut[1]] += 1
        for mut in deletions.values():
            if len(mut) > 1:
                self.check = False
        for mut in insertions.values():
            if len(mut) > 1:
                self.check = False
        for mut in mismatches.values():
            if len(mut) > 1:
                self.check = False
        if len(mismatches)>0:
            print(mismatches)
            print(self.start)
            print(self.start_pre)



            # if read_info.mismatches
            # print(read_info.deletions)
            # print(read_info.insertions)
            # print(read_info.mismatches)
            if read_info.repeat_length not in ms_dis:
                ms_dis[read_info.repeat_length] = 1
            else:
                ms_dis[read_info.repeat_length] += 1
        self.ms_dis = ms_dis
