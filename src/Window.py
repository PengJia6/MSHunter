#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : Window.py
# Author : Peng Jia
# Date   : 2020.08.04
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
import multiprocessing
from src.global_dict import *
from src.Microsatellite import Microsatellite
from src.Read import Read
import pysam


# import Read

class Window:

    def __init__(self, contig, ms_info_list):
        self.contig = contig

        self.paras = get_value("paras")
        self.ms_list = ms_info_list

        self.bam_path = self.paras["input"]
        self.threads = self.paras["threads"]
        self.win_start = ms_info_list[0]["pos"] - self.paras["prefix_len"]
        self.win_end = ms_info_list[-1]["pos"] + ms_info_list[-1]["repeatTimes"] * ms_info_list[-1]["motifLen"] + \
                       self.paras["suffix_len"]
        self.reads = {}
        self.reads_num = 0
        self.microsatellites = {}
        self.microsatellites_id = {}
        logger.info("Processing " + contig + " " + str(self.win_start) + "-" +
                    str(self.win_end) + "\t Microsatellites: " + str(len(ms_info_list)))

        # self.ms_list_pos=[ () for info in ms_info_list]

    def init_one_microsatellite(self, ms_info):
        ms = Microsatellite(ms_info)
        return ms

    def init_microsatellites(self):
        pool = multiprocessing.Pool(processes=self.threads)
        microsatellites = pool.map(self.init_one_microsatellite, self.ms_list)
        self.microsatellites = {ms_info.ms_id: ms_info for ms_info in microsatellites}
        pool.close()
        pool.join()

    def init_reads(self):
        reads = {}
        sam_file = pysam.AlignmentFile(self.paras["input"], mode="rb", reference_filename=self.paras["reference"])
        # pysam.AlignmentFile(self.bam_path, mode="rb", reference_filename=self.reference_path)
        for ms_id, ms_info in self.microsatellites.items():
            for alignment in sam_file.fetch(ms_info.chrom, ms_info.start_pre - 1, ms_info.end_suf + 1):
                if alignment.is_unmapped or alignment.is_duplicate or alignment.is_secondary:
                    continue
                if alignment.reference_start > ms_info.start_pre or alignment.reference_end < ms_info.end_suf:
                    continue
                if len(alignment.query_sequence) < 2:
                    continue
                read_id = alignment.query_name + "_" + str(alignment.reference_start)
                if read_id not in reads:
                    reads[read_id] = Read(read_id=read_id,
                                          chrom=self.contig,
                                          alignment=alignment,
                                          reference=self.paras["reference"])
                if ms_info.ms_id not in reads[read_id].support_microsatellites:
                    reads[read_id].support_microsatellites.append(ms_info.ms_id)
            # print(reads)
            self.reads = reads
        self.reads_num = len(self.reads)
        logger.info("Processing " + self.contig + " " + str(self.win_start) + "-" +
                    str(self.win_end) + "\tReads: " + str(self.reads_num))

    def get_one_read_info(self, read):

        read.microsatellites = {ms_id: self.microsatellites[ms_id] for ms_id in read.support_microsatellites}
        read.get_read_str()
        read.get_ms_length_one_read()
        return read
        pass  # self.ms_list = [ms.ms_id for ms in result_list]

    def get_reads_info(self):
        pool = multiprocessing.Pool(processes=self.threads)
        # print(self.reads)
        result_list = pool.map(self.get_one_read_info, self.reads.values())
        pool.close()
        pool.join()
        # result_list = []
        # for item in self.reads.values():
        #     result_list.append(self.get_one_read_info(item))

        self.reads = {read.read_id: read for read in result_list}
        pass

    def merge_reads_info(self):
        microsatellites_dict = {ms_id: {} for ms_id in self.microsatellites}
        for read_id, read in self.reads.items():
            for ms_id, ms_read_mut in read.microsatellites.items():
                microsatellites_dict[ms_id][read_id] = ms_read_mut
        # print(microsatellites_dict)
        for ms_id, reads_info in microsatellites_dict.items():
            # print("kkk",reads_info)
            self.microsatellites[ms_id].reads_info = reads_info

    def genotype_one_microsatellite(self, microsatellite):
        microsatellite.get_dis()
        return microsatellite

    def genotype_microsatellite(self):
        pool = multiprocessing.Pool(processes=self.threads)
        microsatellites = pool.map(self.genotype_one_microsatellite, self.microsatellites.values())
        print(len(microsatellites))
        pool.close()
        pool.join()
        self.microsatellites = {ms.ms_id: ms for ms in microsatellites}
        # for i in microsatellites:
        #     if i.depth>0:
        #         print(i.depth)
        # pass

    def write_to_vcf(self):
        pass

    def run_window(self):
        self.init_microsatellites()  # 并行
        self.init_reads()  # 扫描read 确实其对应的 MS
        self.get_reads_info()  # 处理read 并行
        self.merge_reads_info()  # 合并read信息为MS信息
        self.genotype_microsatellite()  # 变异检测 并行
        self.write_to_vcf()  # 一条一条写入

        # print()
        pass
