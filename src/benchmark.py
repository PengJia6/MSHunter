#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : benchmark.py
# Author : Peng Jia
# Date   : 2020.07.13
# Email  : pengjia@stu.xjtu.edu.cn
# Description: Build benchmark for microsatellite mutation calling
=============================================================================="""
import os
import re
import collections
import pysam
import multiprocessing
from src.global_dict import *
from src.units import load_microsatellites


# from src.call import *
# from src.errEval import *


class MSHAP:
    # TODO : normalize the argument
    # TODO : finished all member variables
    """
      chrom : chromsome
      pos_start: start position of microsatellite
      motif: repeat unit of microsatellite 
      motif_len = repeat unit length of microsatellite
      repeat_times =  the repeat times of microsatellite
      pos_end :end position of microsatellite
      bam_path: bam file path
    """
    chrom = ""
    pos_start = 0
    motif = 0
    motif_len = 0
    repeat_times = 0
    pos_end = 0
    bam_path = ""
    reference_path = ""
    start_pre = 0
    end_suf = 0
    ref_repeat_length = 0
    prefix_len = 10  # TODO add in command
    suffix_len = 10  # TODO add in command
    repeat_length_dis = {}
    query_repeat_length = 0
    dis_stat = False
    more_than_one_alleles = False
    more_than_one_alleles_ms = False
    check = True
    check_stats = []
    ms_mismatch = []
    ms_var_type = []
    mutation_id = ""
    ref_str = ""
    alt_str = ""
    low_support = False
    depth = 0

    def __init__(self, chrom, pos_start, pos_end, motif, motif_len, repeat_times,
                 bam_path, reference_path):
        """
        @param chrom: chromsome
        @param pos_start: start position
        @param pos_end: end position
        @param motif: repeat unit
        @param motif_len: repeat unit length
        @param repeat_times: repeat times
        @param bam_path: bam path
        @param reference_path: reference path
        """
        self.chrom = chrom
        self.pos_start = pos_start
        self.motif = motif
        self.motif_len = motif_len
        self.repeat_times = repeat_times
        self.pos_end = pos_end
        self.bam_path = bam_path
        self.reference_path = reference_path
        self.start_pre = pos_start - self.prefix_len
        self.end_suf = pos_end + self.suffix_len
        self.ref_repeat_length = repeat_times * motif_len

    def get_reads_alignment(self, bam_file):

        alignment_list = [alignment for alignment in bam_file.fetch(self.chrom, self.pos_start - 1, self.pos_end + 1)]
        # bam_file.close()
        reads_com = []
        for alignment in alignment_list:
            if alignment.is_unmapped or alignment.is_duplicate or alignment.is_secondary:
                continue
            if alignment.reference_start > self.pos_start or alignment.reference_end < self.pos_end:
                continue
            reads_com.append(alignment)
        return reads_com

    def get_dis2(self):

        bam_file = pysam.AlignmentFile(self.bam_path, mode="rb",
                                       reference_filename=self.reference_path)
        # print(bam_file.count(self.chrom, start=self.pos_start, stop=self.pos_end))
        # return 2
        # print(type(bam_file.pileup(self.chrom,self.start_pre,self.end_suf,truncate=True)))

        variants = {}
        pos = self.start_pre
        for pot in bam_file.pileup(self.chrom,
                                   self.start_pre,
                                   self.end_suf,
                                   truncate=True,
                                   fastafile=pysam.FastaFile(self.reference_path)):
            # pots.append(pot)
            pos += 1
            pot_alleles = list(set(map(lambda x: x.upper(),
                                       pot.get_query_sequences(mark_matches=True,
                                                               mark_ends=False,
                                                               add_indels=True))))
            pot_alleles = collections.Counter(pot_alleles).most_common()
            if len(pot_alleles) > 1:
                self.more_than_one_alleles = True
            if pot_alleles[0][0] in [",", ".", "*"]:
                continue
            if pot_alleles[0][0] in ["A", "G", "C", "T"]:
                variants[pos - 1] = {"Mismatch": pot_alleles[0][0]}
                # print("Mismatch",pot_alleles)
                continue
            p = re.compile("[\\+\\-][0-9]+")
            # print(pot_alleles[0])
            indelf = p.findall(pot_alleles[0][0])[0]
            indel_type = "I" if pot_alleles[0][0][1] == "-" else "D"
            # print(indel_type,indelf)
            indel_len = int(indelf[1:])
            indel_str = pot_alleles[0][0][-indel_len:]
            # print(indel_len,indel_type,indel_str)
            if pot_alleles[0][0][0] in ["A", "G", "C", "T"]:
                variants[pos - 1] = {"Mismatch": pot_alleles[0][0][0], indel_type: [indel_len, indel_str]}
            else:
                variants[pos - 1] = {indel_type: [indel_len, indel_str]}

    def get_dis(self):
        bam_file = pysam.AlignmentFile(self.bam_path, mode="rb", reference_filename=self.reference_path)
        repeat_length_dict = {}
        reads = self.get_reads_alignment(bam_file)
        for alignment in reads:
            repeat_length = self.get_repeat_length(alignment)
            if repeat_length < 0: continue
            if repeat_length not in repeat_length_dict:
                repeat_length_dict[repeat_length] = 0
            repeat_length_dict[repeat_length] += 1
        self.repeat_length_dis = repeat_length_dict
        self.depth = len(repeat_length_dict)
        if self.depth > 1:
            self.more_than_one_alleles = True
            self.more_than_one_alleles_ms = True
            self.check = False
            self.check_stats.append("More_alleles_in_MS")
        if len(repeat_length_dict) > 0:
            self.query_repeat_length = list(repeat_length_dict.values())[0]
            self.dis_stat = True
        else:
            self.check = False
            self.check_stats.append("No_read_covered")
        if not self.check:
            self.ms_var_type.append("Fuzzy")
            return -1

        fa_file = pysam.FastaFile(self.reference_path)
        pos = self.pos_start
        segment_pos = 0
        self.ref_str = fa_file.fetch(self.chrom, self.pos_start, self.pos_end)
        alt_str = []
        for pot in bam_file.pileup(self.chrom,
                                   self.pos_start,
                                   self.pos_end,
                                   truncate=True,
                                   fastafile=fa_file):
            pot_alleles = list(set(map(lambda x: x.upper(),
                                       pot.get_query_sequences(mark_matches=True,
                                                               mark_ends=False,
                                                               add_indels=True))))
            pot_alleles = collections.Counter(pot_alleles).most_common()
            if pot_alleles[0][0][0] in [",", "."]:
                alt_str.append(self.ref_str[segment_pos])
            elif pot_alleles[0][0][0] in ["A", "G", "C", "T"]:
                alt_str.append(pot_alleles[0][0][0])
                self.ms_mismatch.append([pos, pot_alleles[0][0][0], fa.fetch(self.chrom, pos, pos + 1)])
            if len(pot_alleles[0][0]) > 1 and pot_alleles[0][0][1] == "+":
                p = re.compile("[\\+\\-][0-9]+")
                indel_f = p.findall(pot_alleles[0][0])[0]
                # print(indel_type,indelf)
                indel_len = int(indel_f[1:])
                indel_str = pot_alleles[0][0][-indel_len:]
                alt_str.append(indel_str)
                # print("indser",pot_alleles)
            pos += 1
            segment_pos += 1
        bam_file.close()
        self.alt_str = "".join(alt_str)
        if self.ref_repeat_length != self.query_repeat_length:
            indel_type = "DEL" if self.ref_repeat_length > self.query_repeat_length else "INS"
            if len(self.ms_mismatch) > 0:
                self.ms_var_type.append("Complex")
                self.ms_var_type.append(indel_type)
                self.ms_var_type.append("SNP")
            else:
                self.ms_var_type.append(indel_type)
        else:
            if len(self.ms_mismatch) > 0:
                self.ms_var_type.append("SNP")
            else:
                self.ms_var_type.append("None")
        print(self.ref_str,self.alt_str)

    def compile(self, pots):
        seq = []
        seq_error = {}
        for pot in pots:
            thisseq = []
            for pp in pot:
                thispp = pp.upper()
                if thispp not in thisseq:
                    thisseq = thispp

    def get_alignment_detail(self):
        # print(self.dis_stat)
        if not self.dis_stat:
            return -1

        bam_file = pysam.AlignmentFile(self.bam_path, mode="rb", reference_filename=self.reference_path)

        # print(type(bam_file.pileup(self.chrom,self.start_pre,self.end_suf,truncate=True)))

        variants = {}
        pos = self.start_pre
        for pot in bam_file.pileup(self.chrom,
                                   self.start_pre,
                                   self.end_suf,
                                   truncate=True,
                                   fastafile=pysam.FastaFile(self.reference_path)):
            # pots.append(pot)
            pos += 1
            pot_alleles = list(set(map(lambda x: x.upper(),
                                       pot.get_query_sequences(mark_matches=True,
                                                               mark_ends=False,
                                                               add_indels=True))))
            pot_alleles = collections.Counter(pot_alleles).most_common()
            if len(pot_alleles) > 1:
                self.more_than_one_alleles = True
            if pot_alleles[0][0] in [",", ".", "*"]:
                continue
            if pot_alleles[0][0] in ["A", "G", "C", "T"]:
                variants[pos - 1] = {"Mismatch": pot_alleles[0][0]}
                # print("Mismatch",pot_alleles)
                continue
            p = re.compile("[\\+\\-][0-9]+")
            # print(pot_alleles[0])
            indelf = p.findall(pot_alleles[0][0])[0]
            indel_type = "I" if pot_alleles[0][0][1] == "-" else "D"
            # print(indel_type,indelf)
            indel_len = int(indelf[1:])
            indel_str = pot_alleles[0][0][-indel_len:]
            # print(indel_len,indel_type,indel_str)
            if pot_alleles[0][0][0] in ["A", "G", "C", "T"]:
                variants[pos - 1] = {"Mismatch": pot_alleles[0][0][0], indel_type: [indel_len, indel_str]}
            else:
                variants[pos - 1] = {indel_type: [indel_len, indel_str]}
        if len(variants) > 0:
            print(variants)

            # print(pot_alleles)

            # print(pot_alleles)

            # pots.append()
            # print(type(pot))
            # print(pot)
            # print(pot.indel())
            # if pot.is_del:
            #     print("del",pot.get_query_sequences())
            # if pot.is_refskip:
            #     print("refskip", pot.get_query_sequences())
            # print(pot.get_query_sequences())
            pass
        # for i in pots:
        #     if i not in ["A","C","G","T"]:
        #         print(pots)
        #         break
        # print(len(pots),self.start_pre-self.end_suf)
        bam_file.close()

        # print()

    def pos_convert_ref2read(self, ref_block: list, read_block: list, pos: int, direction="start") -> tuple:
        """
        @param direction:  start of end
        @param ref_block:
        @param read_block:
        @param start:
        @param end:
        @return:
        """
        # print("=================================")

        if direction == "start":
            block_index = 0
            for sub_ref_block in ref_block:
                if sub_ref_block[0] <= pos <= sub_ref_block[1]:
                    # print(sub_ref_block,pos)
                    break
                block_index += 1

            if pos == ref_block[block_index][1]:  # M|D M|I D|M
                read_pos = read_block[block_index][1]
                if ref_block[block_index][2] == 2:
                    self.start_pre = min([self.start_pre, ref_block[block_index][0]])
                return pos, read_pos

            if ref_block[block_index][2] == 0:  # match
                read_pos = read_block[block_index][0] + (pos - ref_block[block_index][0])
                return pos, read_pos
            elif ref_block[block_index][2] == 2:  # deletion
                pos = pos
                read_pos = read_block[block_index][1]
                self.start_pre = min([read_block[block_index][0], self.start_pre])

                # pos = ref_block[block_index - 1][1] - 1
                # read_pos = read_block[block_index - 1][1] - 1
                return pos, read_pos
            else:
                print("gjjgjdgffkfkfkfpooooo")
                pos = ref_block[block_index - 1][1] - 1
                read_pos = read_block[block_index - 1][1] - 1
                return pos, read_pos

        if direction == "end":
            block_index = len(ref_block) - 1
            for sub_ref_block in ref_block[::-1]:
                if sub_ref_block[0] <= pos <= sub_ref_block[1]:
                    # print(sub_ref_block,pos)
                    break
                block_index = block_index - 1
            # if pos==ref_block[block_index][1]:
            #     print(block_index,len(ref_block))
            #     print( pos,ref_block[block_index-1],ref_block[block_index],)
            #     print("end")
            if pos == ref_block[block_index][0]:  # D|M I|M M|D
                read_pos = read_block[block_index][0]
                if ref_block[block_index][2] == 2:
                    self.end_suf = max([self.end_suf, ref_block[block_index][1]])
                return pos, read_pos

            if ref_block[block_index][2] == 0:
                read_pos = read_block[block_index][0] + (pos - ref_block[block_index][0])
                return pos, read_pos
            elif ref_block[block_index][2] == 2:
                pos = pos
                read_pos = read_block[block_index][0]
                self.end_suf = max([self.end_suf, read_block[block_index][1]])
                # pos = ref_block[block_index + 1][0] + 1
                # read_pos = read_block[block_index + 1][0] + 1
                return pos, read_pos
            else:
                print("lfllflflfllflflffkfjfhvjfnjgffjfjjjjjjjj")
                # ref_block[block_index][0] == pos:
                # print("llff")
                pos = pos + 1
                read_pos = read_block[block_index - 1][0] - 1
                return pos, read_pos

    def get_repeat_length(self, alignment):
        align_start = alignment.reference_start
        ref_block = []
        read_block = []
        read_pos = 0
        ref_pos = align_start
        for cigartuple in alignment.cigartuples:
            if cigartuple[0] in [0, 7, 8]:  # 0 : M : match or mishmatch ; 7: :=:match; 8:X:mismatch
                ref_block.append((ref_pos, ref_pos + cigartuple[1], 0))
                read_block.append((read_pos, read_pos + cigartuple[1], 0))
                read_pos += cigartuple[1]
                ref_pos += cigartuple[1]
            elif cigartuple[0] in [1, 4, 5]:  # 1:I:inserion ;4:S:soft clip 5:H:hardclip
                ref_block.append((ref_pos, ref_pos + 0, 1))
                read_block.append((read_pos, read_pos + cigartuple[1], 1))
                read_pos += cigartuple[1]
            elif cigartuple[0] in [2, ]:  # 2:D; 3:N: skip region of reference
                ref_block.append((ref_pos, ref_pos + cigartuple[1], 2))
                read_block.append((read_pos, read_pos, 2))
                ref_pos += cigartuple[1]
            else:
                return -1

        if self.pos_start >= ref_block[-1][1] or self.pos_start <= ref_block[0][0]:
            # print("pos2")
            # print(self.pos_start)
            # print(ref_block)
            return -1
        if self.pos_end >= ref_block[-1][1] or self.pos_end <= ref_block[0][0]:
            # print("pos3")
            return -1

        ref_start, read_start = self.pos_convert_ref2read(ref_block, read_block, self.pos_start, direction="start")
        ref_end, read_end = self.pos_convert_ref2read(ref_block, read_block, self.pos_end, direction="end")
        rpt = self.repeat_times * self.motif_len + ((read_end - read_start) - (ref_end - ref_start))

        return rpt

    def dis_sum(self, dict_list):
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

    def get_snp_info(self):
        pre_content = 5
        suf_content = 5
        start = self.pos_start
        end = self.pos_end
        fafile = pysam.FastaFile(self.reference_path)
        start_pos = start - pre_content
        end_pos = end + suf_content
        ref_seq = fafile.fetch(self.chrom, start_pos, end_pos)
        bam_file = pysam.AlignmentFile(self.bam_path)
        print('++++++++++++++++++')
        outfile = pysam.AlignmentFile("-", "w", template=bam_file, index_filename=self.bam_path + ".bai")
        for pot in bam_file.fetch(self.chrom, start_pos, end_pos):
            outfile.write(pot)
        for pot in outfile.pileup(self.chrom, start_pos, end_pos, truncate=True, index_filename=self.bam_path + ".bai"):
            print(pot)

        return

    def getrepeat_times2(self, alignment):

        """
        :param alignment:
        :param motif:
        :param motif_len:
        :param prefix:
        :param suffix:
        :return:
        """

        # self.getrepeat_times2(alignment)
        if alignment.mapping_quality < self.min_mapping_qual:
            return -1
        readString = alignment.query
        prefixState = readString.find(self.prefix)
        if prefixState < 0: return -1
        suffixState = readString.rfind(self.suffix)
        if suffixState < 0: return -3
        if prefixState + 5 >= suffixState: return -2
        while prefixState >= 0:
            count = 0
            start = prefixState + 5
            while start == readString.find(self.motif, start):
                count += 1
                start = readString.find(self.motif, start) + self.motif_len
            if (self.motif_len == 1 and count >= 1) or (self.motif_len > 1 and count >= 1):
                if start == readString.find(self.suffix, start):
                    return count
            prefixState = readString.find(self.prefix, prefixState + 1)
        return -4


def benchmark_init(args):
    """
    argument procress
    """
    paras = {}
    paras["input"] = args.input[0]
    paras["output"] = args.output[0]
    paras["microsatellite"] = args.microsatellite[0]
    paras["reference"] = args.reference[0]
    paras["separator"] = args.separator[0]
    paras["prefix_len"] = args.prefix_len[0]
    paras["suffix_len"] = args.suffix_len[0]
    paras["debug"] = args.debug[0]
    paras["only_homopolymer"] = args.only_homopolymers[0]
    paras["minimum_support_reads"] = args.minimum_support_reads[0]
    paras["threads"] = args.threads[0]
    paras["batch"] = args.batch[0]
    paras["only_microsatellites"] = args.only_microsatellites[0]

    paras["ranges_of_repeat_times"] = {}

    for i in args.minimum_repeat_times[0].split(";"):
        unitRange, repeatRange = i.split(":")
        if "-" in unitRange:
            unitStart, unitEnd = tuple(map(int, unitRange.split("-")))
        else:
            unitStart = int(unitRange)
            unitEnd = unitStart
        repeatStart = int(repeatRange)
        # print(unitStart,unitEnd,"  ",repeatStart, repeatEnd)
        for ur in range(unitStart, unitEnd + 1):
            if ur not in paras["ranges_of_repeat_times"]:
                paras["ranges_of_repeat_times"][ur] = {}
            paras["ranges_of_repeat_times"][ur]["min"] = repeatStart
        for i in args.maximum_repeat_times[0].split(";"):
            # print(i)
            unitRange, repeatRange = i.split(":")
            if "-" in unitRange:
                unitStart, unitEnd = tuple(map(int, unitRange.split("-")))
            else:
                unitStart = int(unitRange)
                unitEnd = unitStart
            repeatStart = int(repeatRange)
            # print(unitStart,unitEnd,"  ",repeatStart, repeatEnd)
            for ur in range(unitStart, unitEnd + 1):
                if ur not in paras["ranges_of_repeat_times"]:
                    paras["ranges_of_repeat_times"][ur] = {}
                paras["ranges_of_repeat_times"][ur]["max"] = repeatStart
    error_stat = False
    if os.path.exists(paras["input"]):
        print("[INFO] The input file is : '" + paras["input"] + "'.")
    else:
        print('[ERROR] The input file '
              + paras["input"] + ' is not exist, please check again')
        error_stat = True

    if os.path.isfile(paras["microsatellite"]):
        print("[INFO] The microsatellites file  is : " + paras["microsatellite"])
    else:
        print('[ERROR] The microsatellites file '
              + paras["microsatellite"] + ' is not exist, please check again')
        error_stat = True
    if os.path.isfile(paras["reference"]):
        print("[INFO] The reference file is : '" + paras["reference"] + "'.")
    else:
        paras["reference"] = "" if paras["reference"] == "." else paras["reference"]
        print('[ERROR] The reference file ' + paras["reference"] + ' is not exist, please check again')
        error_stat = True
    if paras["input"][-4:] == "cram":
        paras["input_format"] = "cram"
        cramfile = pysam.AlignmentFile(paras["input"], mode="rb", reference_filename=paras["reference"])
        if not cramfile.has_index():
            print("[INFO] Build index for the input cram ...")
            pysam.index(paras["input"])
        cramfile.close()
    else:
        paras["input_format"] = "hap1"
        bam_file = pysam.AlignmentFile(paras["input"], mode="rb")
        if not bam_file.has_index():
            print("[INFO] Build index for the input bam ...")
            pysam.index(paras["input"])
        bam_file.close()
    if not os.path.exists(paras["output"]):
        print("[INFO] The output is : " + paras["output"] + ".")
    else:
        print(
            '[ERROR] The output ' + paras["output"] +
            ' is still exist! in case of overwrite files in this workspace, '
            'please check your script!')
        if not paras["debug"]:
            error_stat = True
    if error_stat: return False
    output_path = paras["output"]
    output_path = output_path if output_path[-1] == "/" else output_path + "/"
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    paras["output"] = output_path
    input_path = paras["input"]
    input_path = input_path[:-1] if input_path[-1] == "/" else input_path
    case = input_path.split("/")[-1].strip(".bam")
    case = case.strip(".cram")
    paras["output_dis"] = paras["output"] + case + ".dis.vcf.gz"
    paras["output_tmp"] = paras["output"] + case + "_tmp"
    if not os.path.exists(paras["output_tmp"]):
        os.makedirs(paras["output_tmp"])
    paras["output_model"] = paras["output"] + case + ".model"
    paras["output_call"] = paras["output"] + case + ".vcf.gz"
    set_value("case", case)
    set_value("paras", paras)
    return True


def bm_processOneMs(msDetail):
    msDetail.get_dis()
    msDetail.get_dis2()
    if not get_value("default")["benchmark"]["only_microsatellites"]:
        msDetail.get_alignment_detail()
    # msDetail.calcuShiftProbability()
    # msDetail.get_snp_info()
    return msDetail


def bm_write_vcf_init(outputpath):
    outputfile = pysam.VariantFile(outputpath, "w")
    bam_file = pysam.AlignmentFile(get_value("paras")["input"], "rb")
    contigs = bam_file.references
    contigsLen = bam_file.lengths
    chromList = get_value("chrom_list")
    contigs_len_dict = {}
    sortedContig = []
    for contig, length in zip(contigs, contigsLen):
        if contig in chromList:
            sortedContig.append(contig)
            contigs_len_dict[contig] = length
    for contig in sortedContig:
        outputfile.header.add_line(
            "##contig=<ID={chrom},length={length}>".format(chrom=contig, length=contigs_len_dict[contig]))
    set_value("contigsInfo", contigs_len_dict)
    outputfile.header.add_line('##INFO=<ID=chrom,Number=1,Type=String,Description="Chromosome">')
    outputfile.header.add_line('##INFO=<ID=pos,Number=1,Type=Integer,Description="Position">')
    outputfile.header.add_line('##INFO=<ID=Start,Number=1,Type=Integer,Description="Position start">')
    outputfile.header.add_line('##INFO=<ID=End,Number=1,Type=Integer,Description="Position End">')
    outputfile.header.add_line('##INFO=<ID=motif,Number=1,Type=String,Description="Repeat unit">')
    outputfile.header.add_line('##INFO=<ID=repeat_times,Number=1,Type=Integer,Description="Repeat imes">')
    outputfile.header.add_line('##INFO=<ID=prefix,Number=1,Type=String,Description="Prefix of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=suffix,Number=1,Type=String,Description="Suffix of microsatellite">')
    # outputfile.header.add_line('##INFO=<ID=varType,Number=1,Type=String,Description="Variants Type">')
    outputfile.header.add_line(
        '##INFO=<ID=MSVarType,Number=1,Type=String,Description="Variants type in microsatellite region">')
    outputfile.header.add_line(
        '##INFO=<ID=upstreamVarType,Number=1,Type=String,Description="Variants type in upstream regions of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=dis,Number=1,Type=String,Description='
                               'Distribution of repeat length>')
    outputfile.header.add_line('##INFO=<ID=disStat,Number=1,Type=String,Description="Distribution Stat">')
    return outputfile


def bm_write_vcf_close(outputfile):
    outputfile.close()
    pysam.tabix_index(get_value("paras")["output_dis"], force=True, preset="vcf")


def bm_write_vcf(outputfile, dataList):
    # print(header.contigs)
    # print("write", len(dataList))
    for msDetail in dataList:
        chrom = msDetail.chrom
        pos = int(msDetail.pos_start)
        vcfrec = outputfile.new_record()
        # print("infoKey",vcfrec.info.keys())
        vcfrec.contig = chrom
        # vcfrec.stop = pos + msDetail.repeat_times * len(msDetail.motif)
        vcfrec.pos = pos
        vcfrec.ref = msDetail.ref_str
        vcfrec.alt = msDetail.alt_str
        vcfrec.info["chrom"] = chrom
        vcfrec.info["pos"] = pos
        vcfrec.info["Start"] = pos
        vcfrec.stop = pos + msDetail.repeat_times * len(msDetail.motif)
        vcfrec.info["End"] = pos + msDetail.repeat_times * len(msDetail.motif)
        vcfrec.info["motif"] = msDetail.motif
        vcfrec.info["repeat_times"] = msDetail.repeat_times
        vcfrec.info["prefix"] = msDetail.prefix
        vcfrec.info["suffix"] = msDetail.suffix
        vcfrec.info["dis"] = ":".join(
            [str(key) + "-" + str(value) for key, value in msDetail.repeat_length_dis.items()])
        vcfrec.info["disStat"] = str(msDetail.disStat)

        outputfile.write(vcfrec)


def bm_multiRun(thread, datalist):
    pool = multiprocessing.Pool(processes=thread)
    result_list = pool.map(bm_processOneMs, datalist)
    pool.close()
    pool.join()
    # result_list = []
    # for ms in datalist:
    #     result_list.append(bm_processOneMs(ms))

    # print("input",len(datalist))
    # print("output",len(result_list))

    return result_list


def getDis(args={}, upstreamLen=5, downstreamLen=5):
    dis = args["output_dis"]
    input_format = args["input_format"]
    thread = args["threads"]
    batch = args["batch"]
    outputfile = bm_write_vcf_init(dis)
    contigs_info = get_value("contigsInfo")
    dfMicroSatellites = load_microsatellites(args)
    curentMSNum = 0
    tmpWindow = []
    ms_number = get_value("ms_number")
    fafile = pysam.FastaFile(args["reference"])
    prefix_len = args["prefix_len"]
    suffix_len = args["suffix_len"]
    # print(len(dfMicroSatellites))
    for ms_id, info in dfMicroSatellites.iterrows():
        curentMSNum += 1
        # if curentMSNum < 40000 and args["debug"]:
        #     continue
        chrom = info["chr"]
        if chrom not in contigs_info:
            continue
        pos_start = int(info["pos"])
        repeat_times = int(info["repeatTimes"])
        motif = info["motif"]
        motif_len = len(motif)
        pos_end = pos_start + motif_len * repeat_times
        # queryStart = pos_start - upstreamLen
        # queryEnd = pos_end + downstreamLen
        # prefix_str = fafile.fetch(chrom, pos_start - prefix_len, pos_start)
        # suffix_str = fafile.fetch(chrom, pos_end, pos_end + suffix_len)
        # print(prefix_str,info["prefix"],motif)
        # print(suffix_str,info["suffix"],motif)
        thisMSDeail = MSHAP(chrom=chrom,
                            pos_start=pos_start,
                            pos_end=pos_end,
                            motif=info["motif"],
                            motif_len=motif_len,
                            repeat_times=repeat_times,
                            bam_path=args["input"],
                            # input_format=input_format,

                            # prefix_str=prefix_str,
                            # suffix_str=suffix_str,
                            reference_path=args["reference"],
                            )
        tmpWindow.append(thisMSDeail)

    if curentMSNum % (batch * thread) == 0:
        print("[INFO] Bam2dis: Total", ms_number, "microsatelite, processing:", curentMSNum - batch * thread + 1,
              "-", curentMSNum, "(" + str(round(curentMSNum / ms_number * 100, 2)) + "%)")
        result_list = bm_multiRun(thread=thread, datalist=tmpWindow)
        bm_write_vcf(outputfile, result_list)
        tmpWindow = []
    result_list = bm_multiRun(thread=thread, datalist=tmpWindow)
    fafile.close()
    bm_write_vcf(outputfile, result_list)
    bm_write_vcf_close(outputfile)
    print("[INFO] Bam2dis: Total", ms_number, "microsatelite, finish all!")


def benchmark(parase):
    if benchmark_init(parase):
        args = get_value("paras")
        getDis(args)
        print("hhh")
