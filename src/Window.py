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


class Window:

    def __init__(self, ms_info_list):
        self.contig = ms_info_list[0]["chr"]

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
        self.microsatellites_id = [it["chr"] + "_" + str(it["pos"]) for it in ms_info_list]
        self.vcf_recs = []
        # logger.info("\t--------------------------------------------------------------------------------")
        # logger.info("\tProcessing " + contig + ":" + str(self.win_start) + "-" + str(self.win_end))
        # logger.info("\tNo. of Microsatellites: " + str(len(ms_info_list)))
        # logger.info("Processing " + contig + " " + str(self.win_start) + "-" +
        #             str(self.win_end) + "\t Microsatellites: " + str(len(ms_info_list)))

        # self.ms_list_pos=[ () for info in ms_info_list]

    def init_one_microsatellite(self, ms_info):
        # ref_path=get_value("paras")["reference"]
        ms = Microsatellite(ms_info)
        return ms

    def init_microsatellites(self):
        # pool = multiprocessing.Pool(processes=self.threads)
        # microsatellites = pool.map(self.init_one_microsatellite, self.ms_list)
        # pool.close()
        # pool.join()
        microsatellites = []
        for ms in self.ms_list:
            microsatellites.append(self.init_one_microsatellite(ms))
        self.microsatellites = {ms_info.ms_id: ms_info for ms_info in microsatellites}

    def init_reads(self):
        reads = {}
        sam_file = pysam.AlignmentFile(self.paras["input"], mode="rb", reference_filename=self.paras["reference"])
        # pysam.AlignmentFile(self.bam_path, mode="rb", reference_filename=self.reference_path)
        # TODO optimize
        for ms_id, ms_info in self.microsatellites.items():
            for alignment in sam_file.fetch(ms_info.chrom, ms_info.start_pre - 1, ms_info.end_suf + 1):
                if alignment.is_unmapped or alignment.is_duplicate or alignment.is_secondary or alignment.is_supplementary:
                    continue
                if alignment.reference_start > ms_info.start_pre - 1 or alignment.reference_end < ms_info.end_suf + 1:
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
        self.reads = reads
        self.reads_num = len(self.reads)

    def get_one_read_info(self, read):
        read.microsatellites = {ms_id: self.microsatellites[ms_id] for ms_id in read.support_microsatellites}
        read.get_read_str()
        read.get_ms_info_one_read()
        # read.get_mut_info()
        return read
        pass  # self.ms_list = [ms.ms_id for ms in result_list]

    def get_reads_info(self):
        # logger.info("\tScan reads covered microsatellites ... ")
        # pool = multiprocessing.Pool(processes=self.threads)
        # # print(self.reads)
        # result_list = pool.map(self.get_one_read_info, self.reads.values())
        # pool.close()
        # pool.join()
        result_list = []
        for item in self.reads.values():
            result_list.append(self.get_one_read_info(item))

        self.reads = {read.read_id: read for read in result_list}
        pass

    def merge_reads_info(self):
        # logger.info("\tMerge microsatellites infomation from different reads... ")
        microsatellites_dict = {ms_id: {} for ms_id in self.microsatellites}
        for read_id, read in self.reads.items():
            for ms_id, ms_read_mut in read.mut_info.items():
                microsatellites_dict[ms_id][read_id] = ms_read_mut
        self.reads = {}
        for ms_id, reads_info in microsatellites_dict.items():
            self.microsatellites[ms_id].set_reads_info(reads_info)

    def genotype_one_microsatellite(self, microsatellite):
        # microsatellite.get_dis()
        microsatellite.one_hap_genotype()
        return microsatellite

    def genotype_microsatellite_ccs_contig(self):
        # logger.info("\tMicrosatellites genotyping ... ")
        #
        # pool = multiprocessing.Pool(processes=self.threads)
        # microsatellites = pool.map(self.genotype_one_microsatellite, self.microsatellites.values())
        # pool.close()
        # pool.join()
        microsatellites = []
        for microsatellite in self.microsatellites.values():
            microsatellites.append(self.genotype_one_microsatellite(microsatellite))
        self.microsatellites = {ms.ms_id: ms for ms in microsatellites}

    def write_to_vcf_ccs_contig(self, file_output):
        # logger.info("\tWrite to vcf ... ")
        recs = []
        for ms_id in self.microsatellites_id:
            # print(ms_id)
            # print(self.microsatellites)
            ms = self.microsatellites[ms_id]
            vcfrec = file_output.new_record()
            # print("infoKey",vcfrec.info.keys())
            vcfrec.contig = ms.chrom
            # vcfrec.stop = pos + ms.repeat_times * len(ms.motif)
            vcfrec.pos = ms.mut_start
            vcfrec.ref = ms.ref_str
            vcfrec.alts = (ms.alt_str,) if ms.alt_str != "" else (".",)
            vcfrec.id = ms.ms_id
            vcfrec.stop = ms.mut_end
            vcfrec.info["ms_start"] = ms.start
            vcfrec.info["ms_end"] = ms.end
            vcfrec.info["motif"] = ms.repeat_unit
            vcfrec.info["repeat_times"] = ms.repeat_times
            vcfrec.info["motif_len"] = ms.repeat_unit_len
            vcfrec.info["ref_repeat_length"] = ms.repeat_len
            vcfrec.info["start_pre"] = ms.start_pre
            vcfrec.info["end_suf"] = ms.end_suf
            vcfrec.info["mut_start"] = ms.mut_start
            vcfrec.info["mut_end"] = ms.mut_end

            vcfrec.info["query_repeat_length"] = ms.query_repeat_length
            vcfrec.info["dis_stat"] = str(ms.dis_stat)
            vcfrec.info["check"] = str(ms.check)
            vcfrec.info["check_stats"] = "|".join(ms.check_status)
            vcfrec.info["dis"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis.items()])
            # vcfrec.info["allele"] = ms.allele
            # if ms.check:
            vcfrec.info["var_type"] = ms.mut.var_type
            vcfrec.info["var_type_list"] = ms.mut.var_type_detail
            # print(ms.mut_type.var_prefix,ms.mut_type.var_ms,ms.mut_type.var_suffix)
            vcfrec.info["var_detail"] = ms.mut.var_detail_str
            # file_output.write(vcfrec)
            recs.append(vcfrec)
        return recs

    def run_window_benchmark(self):
        self.init_microsatellites()  # 并行
        self.init_reads()  # 扫描read 确实其对应的 MS
        self.get_reads_info()  # 处理read 并行
        self.merge_reads_info()  # 合并read信息为MS信息
        if self.paras["command"] == "benchmark":
            self.genotype_microsatellite_ccs_contig()  # 变异检测 并行

            # self.write_to_vcf_ccs_contig(file_output)  # 一条一条写入
