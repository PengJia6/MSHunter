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
import numpy as np
import pysam
from src.units import *
from src.Read import Read


class Microsatellite:
    """
    Description: For microsatellite
    """

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
        self.reference = ms_info["reference"]
        self.repeat_len = self.repeat_times * self.repeat_unit_len
        self.end = self.start + self.repeat_len
        self.start_pre = self.start - ms_info["prefix_len"]
        self.end_suf = self.end + ms_info["suffix_len"]
        # self.ref_list = []
        self.reads_info = {}
        self.muts = {}
        self.length_dis_reads = {}
        self.depth = 0
        self.check = True
        self.check_status = []
        self.mut = Mutation()
        self.ref_str = ""
        self.alt_str = ""
        self.dis_stat = True
        self.mut_start = self.start
        self.mut_end = self.end
        self.ms_dis_info = {}
        self.ms_dis = {}
        self.ms_dis_forward = {}
        self.ms_dis_reversed = {}
        self.ms_dis_hap1 = {}
        self.ms_dis_hap2 = {}
        self.ms_dis_hap0 = {}
        self.query_repeat_length = 0
        self.deletions = {}
        self.insertions = {}
        self.mismatches = {}
        self.reads_phased = True  # True if reads in this regions is phased
        self.model_stat = True  # True if model is built in estimate process

        # print(ms_info)

    def set_reads_info(self, reads_info):
        self.reads_info = reads_info
        self.depth = len(self.reads_info)

    def set_muts_info(self, muts):
        self.muts = muts
        self.depth = len(self.muts)

    def set_read_dis_info(self, reads_info):
        # self.reads_info = reads_info
        dis = {}
        dis_stand = {True: {}, False: {}}
        dis_hap = {0: {}, 1: {}, 2: {}}
        # print(reads_info)

        for read_id, ms_info in reads_info.items():
            repeat_length, stand, hap = ms_info
            if repeat_length not in dis:
                dis[repeat_length] = 0
            dis[repeat_length] += 1
            if repeat_length not in dis_hap[hap]:
                dis_hap[hap][repeat_length] = 0
            dis_hap[hap][repeat_length] += 1
            if repeat_length not in dis_stand[stand]:
                dis_stand[stand][repeat_length] = 0
            dis_stand[stand][repeat_length] += 1
        self.depth = len(dis)
        self.dis_stat = True if self.depth > 0 else False
        self.ms_dis = dis
        self.ms_dis_hap0 = dis_hap[0]
        self.ms_dis_hap1 = dis_hap[1]
        self.ms_dis_hap2 = dis_hap[2]
        self.ms_dis_forward = dis_stand[True]
        self.ms_dis_reversed = dis_stand[False]
        self.query_repeat_length = get_max_support_index(dis)
        dis_hap0_num = len(dis_hap[0])
        dis_hap1_num = len(dis_hap[1])
        dis_hap2_num = len(dis_hap[2])
        if abs(dis_hap1_num - dis_hap2_num) > self.depth * 0.4:  # TODO add in input arguments
            self.reads_phased = False
        elif dis_hap0_num < self.depth * 0.5:  # TODO add in input arguments
            self.reads_phased = False
        else:
            self.reads_phased = True

    def get_dis(self):
        samfile = pysam.Samfile(get_value("paras")["input"])
        reads = {}
        for alignment in samfile.fetch(self.chrom, self.start_pre - 1, self.end_suf + 1):
            # print(read)
            if alignment.is_unmapped or alignment.is_duplicate or alignment.is_secondary or alignment.is_supplementary:
                continue
            if alignment.reference_start > self.start_pre - 1 or alignment.reference_end < self.end_suf + 1:
                continue
            if len(alignment.query_sequence) < 2:
                continue
            read_ms = Read(read_id="", alignment=alignment, reference=self.reference, chrom=self.chrom)
            read_ms.get_read_str()
            q_repeat_length = read_ms.get_repeat_length(self.start, self.end)
            if q_repeat_length not in reads:
                reads[q_repeat_length] = 1
            else:
                reads[q_repeat_length] += 1
        return reads

    def get_dis_qual(self):
        samfile = pysam.Samfile(get_value("paras")["input"])
        reads = {}
        quals = {}
        prefix_forward = []
        suffix_forward = []
        ms_forward = []
        prefix_reversed = []
        suffix_reversed = []
        ms_reversed = []

        num_forward = 0
        num_reversed = 0
        for alignment in samfile.fetch(self.chrom, self.start_pre - 1, self.end_suf + 1):
            # print(read)
            if alignment.is_unmapped or alignment.is_duplicate or alignment.is_secondary or alignment.is_supplementary:
                continue
            if alignment.reference_start > self.start_pre - 1 or alignment.reference_end < self.end_suf + 1:
                continue
            if len(alignment.query_sequence) < 2:
                continue
            read_ms = Read(read_id="", alignment=alignment, reference=self.reference, chrom=self.chrom)
            read_ms.get_read_str()
            q_repeat_length = read_ms.get_repeat_length(self.start, self.end)
            if q_repeat_length not in reads:
                reads[q_repeat_length] = 1
            else:
                reads[q_repeat_length] += 1
            prefix_qual = read_ms.get_quals(self.start_pre, self.start)
            suffix_qual = read_ms.get_quals(self.end, self.end_suf)
            ms_qual = read_ms.get_quals(self.start, self.end)
            if prefix_qual[1]:
                num_forward += 1
                prefix_forward.append(np.array(list(map(str2int, prefix_qual[0]))))
                suffix_forward.append(np.array(list(map(str2int, suffix_qual[0]))))
                ms_forward.append(np.array(list(map(str2int, ms_qual[0]))))
            else:
                num_reversed += 1
                prefix_reversed.append(np.array(list(map(str2int, prefix_qual[0]))))
                suffix_reversed.append(np.array(list(map(str2int, suffix_qual[0]))))
                ms_reversed.append(np.array(list(map(str2int, ms_qual[0]))))
        # print(set([len(i) for i in suffix]))
        quals["num_forward"] = num_forward
        quals["prefix_forward"] = np.array(prefix_forward)
        quals["suffix_forward"] = np.array(suffix_forward)
        quals["ms_forward"] = np.array(ms_forward)

        quals["num_reversed"] = num_reversed
        quals["prefix_reversed"] = np.array(prefix_reversed)
        quals["suffix_reversed"] = np.array(suffix_reversed)
        quals["ms_reversed"] = np.array(ms_reversed)
        return reads, quals

    def deletion_merge(self):
        for read_id, read_info in self.muts.items():
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
            self.muts[read_id].deletions = deletions

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

    def get_alt(self, alt_dict, offset):

        alt_list = list(self.ref_str)
        # print(self.ref_str)
        # print(len(alt_list))
        # print(alt_dict)
        # print(offset)
        # print("==============================")
        for pos, info in alt_dict:
            alt_list[pos - offset] = info
        return "".join(alt_list)

    def call_init(self):
        if self.depth == 0:
            self.check = False
            self.check_status.append("No_read_covered")
            self.dis_stat = False
            self.ref_str = pysam.FastaFile(self.reference).fetch(self.chrom, self.mut_start, self.mut_end + 1)
            return False
        model_all = get_value("model")
        if self.repeat_unit not in model_all:
            self.check = False
            self.check_status.append("No_error_model")
            self.model_stat = False
            return False

        model = model_all[self.repeat_unit]
        del model_all
        # self.model_stat = True
        repeat_length_info = {}
        for read_id, read_info in self.reads_info.items():
            print(self.ms_id, read_info.repeat_lengths)
            if self.ms_id not in read_info.repeat_lengths: continue
            print("hhhhh")
            repeat_length_info[read_id] = [read_info.repeat_lengths[self.ms_id], read_info.stand, read_info.hap]
        self.set_read_dis_info(repeat_length_info)

    # def call_init(self):
    #     if self.depth == 0:
    #         self.check = False
    #         self.check_status.append("No_read_covered")
    #         self.dis_stat = False
    #         self.ref_str = pysam.FastaFile(self.reference).fetch(self.chrom, self.mut_start, self.mut_end + 1)
    #         return
    #     self.deletion_merge()
    #     ms_dis = {}
    #     mut_dict_by_pos = {}
    #     for read_id, read_info in self.reads_info.items():
    #         hap = read_info.hap
    #         strand = read_info.strand
    #         for mut in read_info.mismatches:
    #             if mut[0] not in mut_dict_by_pos:
    #                 mut_dict_by_pos[mut[0]] = []
    #             mut_dict_by_pos[mut[0]].append([mut[0], "SNV", read_id, hap, strand, mut])
    #
    #         for mut in read_info.insertions:
    #             if mut[0] not in mut_dict_by_pos:
    #                 mut_dict_by_pos[mut[0]] = []
    #             mut_dict_by_pos[mut[0]].append([mut[0], "INS", read_id, hap, strand, mut])
    #
    #         for mut in read_info.deletions:
    #             if mut[0] not in mut_dict_by_pos:
    #                 mut_dict_by_pos[mut[0]] = []
    #             mut_dict_by_pos[mut[0]].append([mut[0], "DEL", read_id, hap, strand, mut])
    #
    #         if read_info.repeat_length not in ms_dis:
    #             ms_dis[read_info.repeat_length] = 1
    #             if
    #         else:
    #             ms_dis[read_info.repeat_length] += 1
    #     self.ms_dis = ms_dis
    #     self.query_repeat_length = get_max_support_index(ms_dis)
    #
    #
    #
    #
    #     reads_info = {}
    #     for pos, infos in mut_dict_by_pos.items():
    #         support = len(infos)
    #         support_forward = 0
    #         support_reversed = 0
    #         for info in infos:
    #             if info[4]:
    #                 support_forward += 1
    #             else:
    #                 support_reversed += 1
    #         if support < self.depth * 0.2: continue  # remove low frequency mutation
    #         if abs(support_forward - support_reversed) > support * 0.4: continue  # remove strand bias
    #         for info in infos:
    #             if info[2] not in reads_info:
    #                 reads_info[info[2]] = []
    #             reads_info[info[2]].append(info)
    #     pattern_id = 0
    #     patterns = {}
    #     for read_id, read_info in reads_info.items():
    #         pattern = {}
    #         hap = read_info[0][3]
    #         strand = read_info[0][4]
    #         for item in read_info:
    #             pattern[item[0]] = "_".join([str(item[0]), item[1], str(item[5][1])])
    #         new = True
    #         for i_pattern in patterns:
    #             if patterns[i_pattern].contain(pattern):
    #                 patterns[i_pattern].add_reads(read_id, hap, strand)
    #                 new = False
    #                 break
    #         if new:
    #             pattern_id += 1
    #             patterns[pattern_id] = PatternCluster(pattern)
    #             patterns[pattern_id].init_patttern(pattern, read_id, hap, strand)

    # print("===================================================")
    # print(read_id, reads_info)

    # TODO remove noise in reads and processing according reads
    # if len(patterns) == 0:
    #     pattern_id = 0
    #     patterns[pattern_id] = PatternCluster(pattern_id)
    # else:
    #     for pattern_id, pattern in patterns.items():
    #         print(pattern_id)
    #         if pattern
    def remove_noise(self):
        # TODO
        pass

    def ccs_genotype(self):
        # TODO
        pass

    def one_hap_genotype(self):
        if self.depth == 0:
            self.check = False
            self.check_status.append("No_read_covered")
            self.dis_stat = False
            self.ref_str = pysam.FastaFile(self.reference).fetch(self.chrom, self.mut_start, self.mut_end + 1)
            return
        self.deletion_merge()
        mismatches = {}
        deletions = {}
        insertions = {}
        ms_dis = {}
        mut_start = self.start
        mut_end = self.end - 1
        for read_id, read_info in self.muts.items():
            for mut in read_info.mismatches:
                mut_start = min(mut_start, mut[0])
                mut_end = max(mut_end, mut[0])
                if mut[0] not in mismatches:
                    mismatches[mut[0]] = {}
                if mut[1] not in mismatches[mut[0]]:
                    mismatches[mut[0]][mut[1]] = 1
                else:
                    mismatches[mut[0]][mut[1]] += 1
            for mut in read_info.insertions:
                mut_start = min(mut_start, mut[0])
                mut_end = max(mut_end, mut[0])
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
                # print(mut)
                mut_start = min(mut_start, mut[0])
                mut_end = max(mut_end, mut[0] + mut[1])

            if read_info.repeat_length not in ms_dis:
                ms_dis[read_info.repeat_length] = 1
            else:
                ms_dis[read_info.repeat_length] += 1
        self.mut_start = mut_start
        self.mut_end = mut_end
        self.ms_dis = ms_dis
        self.query_repeat_length = get_max_support_index(ms_dis)
        for mut in deletions.values():
            if len(mut) > 1:
                self.check = False
                self.check_status.append("More_alleles_in_MS")
        for mut in insertions.values():
            if len(mut) > 1:
                self.check = False
                self.check_status.append("More_alleles_in_MS")
        for mut in mismatches.values():
            if len(mut) > 1:
                self.check = False
                self.check_status.append("More_alleles_in_MS")

        # mutation = Mutation()
        # if self.check:
        alt_list = []
        for mut_pos, info in mismatches.items():
            gt_list = list(info.keys())
            alt_list.append([mut_pos, gt_list[0]])
            if mut_pos < self.start:
                self.mut.var_pre[mut_pos] = ["SNV", gt_list]
            elif mut_pos < self.end:
                self.mut.var_ms[mut_pos] = ["SNV", gt_list]
            else:
                self.mut.var_suf[mut_pos] = ["SNV", gt_list]
        for mut_pos, info in insertions.items():
            gt_list = list(info.keys())
            alt_list.append([mut_pos, gt_list[0]])
            if mut_pos < self.start - 1:
                self.mut.var_pre[mut_pos] = ["INS", gt_list]
            elif mut_pos < self.end:
                self.mut.var_ms[mut_pos] = ["INS", gt_list]
            else:
                self.mut.var_suf[mut_pos] = ["INS", gt_list]
        for mut_pos, info in deletions.items():
            gt_list = list(info.keys())
            for del_pos in range(mut_pos, mut_pos + gt_list[0]):
                alt_list.append([del_pos, ""])
            del_end = mut_pos + gt_list[0] - 1
            if del_end < self.start:
                self.mut.var_pre[mut_pos] = ["DEL", gt_list]
            elif mut_pos >= self.end:
                self.mut.var_suf[mut_pos] = ["DEL", gt_list]
            else:
                if del_end < self.end or mut_pos >= self.start:
                    self.mut.var_suf[mut_pos] = ["DEL", gt_list]
                else:
                    self.check = False
                    self.check_status.append("DEL_Span")
        self.mut.compute()
        self.ref_str = pysam.FastaFile(self.reference).fetch(self.chrom, self.mut_start - 1, self.mut_end + 1)
        self.alt_str = "." if self.mut.var_type == "None" else self.get_alt(alt_list, offset=self.mut_start - 1)
