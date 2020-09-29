#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : Read.py
# Author : Peng Jia
# Date   : 2020.08.05
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
from src.units import *
import pysam
import numpy as np


class Read:
    """
    Description: Read
    """

    def __init__(self, read_id, alignment, reference, chrom, tech=""):
        self.chrom = chrom
        self.read_name = alignment.query_name
        self.align_start = alignment.reference_start
        self.align_end = alignment.reference_end
        self.this_read_str = alignment.query_sequence.upper()
        self.tech = tech
        # print(alignment)
        if tech == "contig":
            self.this_read_quals = []
        else:
            self.this_read_quals = "".join([chr(i + 33) for i in alignment.query_qualities])
        self.strand = False if alignment.is_reverse else True  # True: forward False: reverse
        self.this_read_list = []
        self.this_quals_list = []
        self.this_ref_str = ""
        self.this_ref_list = []
        self.read_id = read_id
        self.reference = reference
        self.support_microsatellites = []
        # self.alignment = alignment
        self.microsatellites = {}
        self.direction = False if alignment.is_reverse else True
        self.hap = int(alignment.get_tag("HP")) if alignment.has_tag("HP") else 0
        # print(self.hap)
        self.cigartuples = alignment.cigartuples
        self.mut_info = {}
        self.repeat_lengths = {}

    # def get_microsatellite_detail(self, ms_info):
    #     self. = ms_info

    def get_read_str(self):
        pass
        # TODO:  To be optimized

        self.this_ref_str = pysam.FastaFile(self.reference).fetch(self.chrom, start=self.align_start,
                                                                  end=self.align_end).upper()
        self.this_ref_list = list(self.this_ref_str)
        sub_read_str = []
        sub_read_quals = []

        # read_mut_info = ReadInfo()
        # read_mut_info.direction = self.direction
        # read_mut_info.hap = self.hap
        read_pos = 0
        for cigartuple in self.cigartuples:
            if cigartuple[0] in [0, 7, 8]:  # 0 : M : match or mishmatch ; 7: :=:match; 8:X:mismatch
                match_read = list(self.this_read_str[read_pos:read_pos + cigartuple[1]])
                match_quals = list(self.this_read_quals[read_pos:read_pos + cigartuple[1]])
                sub_read_str.extend(match_read)
                sub_read_quals.extend(match_quals)
                read_pos += cigartuple[1]
            elif cigartuple[0] in [1, 4, 5]:  # 1:I:inserion ;4:S:soft clip 5:H:hardclip
                if cigartuple[0] == 1:
                    if len(sub_read_str) == 0:
                        continue
                    else:
                        sub_read_str[-1] += self.this_read_str[read_pos:read_pos + cigartuple[1]]
                        if self.tech == "contig": continue
                        sub_read_quals[-1] += self.this_read_quals[read_pos:read_pos + cigartuple[1]]
                elif cigartuple[0] == 5:
                    continue
                read_pos += cigartuple[1]
            elif cigartuple[0] in [2, 3]:  # 2:D; 3:N: skip region of reference
                sub_read_str.extend([""] * cigartuple[1])
                if self.tech == "contig": continue
                sub_read_quals.extend([""] * cigartuple[1])
            else:
                return -1
        self.this_read_list = sub_read_str
        self.this_read_str = ""
        self.this_quals_list = sub_read_quals
        self.this_read_quals = ""

    def get_repeat_length(self, ms_start, ms_end):  # give a start and end
        query_repeat_length = len(
            "".join(self.this_read_list[ms_start - 1 - self.align_start:ms_end - self.align_start - 1]))
        return query_repeat_length

    def get_repeat_length_all_ms(self):  # return all microsatellite covered
        self.microsatellites_num = len(self.microsatellites)
        repeat_lengths = {}
        for ms_id, ms in self.microsatellites.items():
            ms_start = ms.start
            ms_end = ms.end
            query_repeat_length = len(
                "".join(self.this_read_list[ms_start - 1 - self.align_start:ms_end - self.align_start - 1]))
            repeat_lengths[ms_id] = query_repeat_length
        self.repeat_lengths = repeat_lengths

    def get_quals(self, q_start, q_end):
        quals_list = self.this_quals_list[q_start - self.align_start - 1:q_end - self.align_start - 1]
        return quals_list, self.strand

    def get_ms_info_one_read(self):
        self.microsatellites_num = len(self.microsatellites)
        # print(self.read_id, self.microsatellites_num,len(self.support_microsatellites))
        read_muts = {}
        repeat_lengths = {}
        for ms_id, ms in self.microsatellites.items():
            ms_start = ms.start
            ms_end = ms.end
            ms_start_pre = ms.start_pre
            ms_end_suf = ms.end_suf
            query_repeat_length = len(
                "".join(self.this_read_list[ms_start - 1 - self.align_start:ms_end - self.align_start - 1]))
            mismatches = []
            deletions = []
            insertions = []
            ref_pos = ms_start_pre - 2
            # if len(self.this_ref_list) == 13746:
            #     print("13746", ms_id)
            # if len(self.this_ref_list) <= ms_end_suf - self.align_start:
            #     print('index out', ms_id, self.read_id, )
            #     print('align', self.align_start, self.align_end)
            #     print('ms_pos', ms_start_pre, ms_end_suf)
            # print(ms_start_pre - self.align_start - 1, ms_end_suf - self.align_start)
            # print(len(self.this_ref_list))
            for pot in range(ms_start_pre - self.align_start - 1, ms_end_suf - self.align_start):
                ref_pos += 1
                this_read_base = self.this_read_list[pot]
                this_ref_base = self.this_ref_list[pot]
                if this_read_base == this_ref_base:
                    continue
                else:
                    this_read_base_len = len(this_read_base)
                    this_ref_base_len = len(this_ref_base)
                    if this_read_base_len == this_ref_base_len:
                        mismatches.append([ref_pos, this_read_base])
                    else:
                        if this_read_base_len < this_ref_base_len:
                            deletions.append([ref_pos, this_read_base])
                        else:
                            insertions.append([ref_pos, this_read_base])
            read_muts[ms_id] = Read_Mutation(repeat_length=query_repeat_length, strand=self.strand, hap=self.hap,
                                             mismatches=mismatches, deletions=deletions, insertions=insertions,
                                             )
            repeat_lengths[ms_id] = query_repeat_length
        self.repeat_lengths = repeat_lengths
        self.mut_info = read_muts
