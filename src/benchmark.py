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
from src.units import *
from src.Window import Window
#
#
# class MutationType:
#     """
#     Description: Mutation
#     pos_del_prefix: deletion position in prefix
#     pos_del_ms:  deletion position in ms region
#     pos_del_suffix: deletion position in suffix
#     pos_ins_prefix: insertion position in prefix
#     pos_ins_suffix: insertion position in suffix
#     pos_ins_ms: insertion position in ms region
#     pos_snp_prefix: mismatch position in prefix
#     pos_snp_ms: mismatch position in ms regions
#     pos_snp_suffix: mismatch position in suffix
#     var_type_prefix: variant type in prefix
#     var_type_suffix: variant type in suffix
#     var_type_ms: variant type in ms region
#     var_type: variant type in all region
#     var_list: variant type list
#     """
#
#     def __init__(self):
#         self.var_prefix = []
#         self.var_suffix = []
#         self.var_ms = []
#         self.var_type_prefix = []
#         self.var_type_suffix = []
#         self.var_type_ms = []
#         self.var_type = ""
#         self.var_type_list = []
#
#     def prefix_var_type(self):
#         for var in self.var_prefix:
#             if len(var[1]) == len(var[2]):
#                 self.var_type_prefix.append("SNP")
#             elif len(var[1]) > len(var[2]):
#                 self.var_type_prefix.append("INS")
#             elif len(var[1]) < len(var[2]):
#                 self.var_type_prefix.append("DEL")
#
#     def suffix_var_type(self):
#         for var in self.var_suffix:
#             if len(var[1]) == len(var[2]):
#                 self.var_type_suffix.append("SNP")
#             elif len(var[1]) > len(var[2]):
#                 self.var_type_suffix.append("INS")
#             elif len(var[1]) < len(var[2]):
#                 self.var_type_suffix.append("DEL")
#
#     def ms_var_type(self, offset):
#
#         if offset < 0:
#             self.var_type_ms.append("DEL")
#         if offset > 0:
#             self.var_type_ms.append("INS")
#
#         for var in self.var_ms:
#             if len(var[1]) == len(var[2]):
#                 self.var_type_ms.append("SNP")
#
#     def comput(self, offset):
#         self.ms_var_type(offset)
#         self.prefix_var_type()
#         self.suffix_var_type()
#         self.var_type_list = self.var_type_prefix + self.var_type_ms + self.var_type_suffix
#         var_num = len(self.var_type_list)
#         if var_num > 1:
#             self.var_type = "Complex"
#         elif var_num == 1:
#             self.var_type = self.var_type_list[0]
#         else:
#             self.var_type = "None"
#
#
# class MSHAP:
#     """
#     @chrom : chromsome
#     @pos_start: start position of microsatellite
#     @motif: repeat unit of microsatellite
#     @motif_len = repeat unit length of microsatellite
#     @repeat_times =  the repeat times of microsatellite
#     @pos_end :end position of microsatellite
#     @bam_path: bam file path
#     @reference_path: reference path
#     @start_pre: start pos of this call, default 10bp upstream the microsatellite
#     @end_suf:  start pos of this call, default 10bp downstream the microsatellite
#     @ref_repeat_length: the reference repeat length
#     @prefix_len: how many bps to detect upstream
#     @suffix_len = how many bps to detect downstream
#     @repeat_length_dis: repeat length distribution
#     @query_repeat_length： the repeat length supported by the most reads
#     dis_stat: True, at least one read covered this microsatellite
#     more_than_one_alleles: more than one reads covered this microsatellite, and have more alleles
#     more_than_one_alleles_ms: more than one reads covered this microsatellite, and have more alleles in microsatellite region
#     check: True, could be used as benchmark
#     check_stats: Why do not be applied as benchmark
#     mismatch: position of mishmatch
#     ms_var_type: variant type in ms region (DEL,INS,SNP)
#     up_var_type: variant type in upstream of ms region (DEL,INS,SNP)
#     down_var_type: variant type in downstream of ms region (DEL,INS,SNP)
#     microsatellite_id: microsatellite_id (chrom_pos)
#     ref_str: reference string
#     alt_str: alternative string
#     allele: allele in this haplotype
#     mut_start: mutation start
#     mut_end: mutation end
#     mut_type: mutation type and details
#     """
#
#     def __init__(self, chrom, pos_start, pos_end, motif, motif_len, repeat_times,
#                  bam_path, reference_path,
#                  prefix_len=0,
#                  suffix_len=0):
#         """
#         @param chrom: chromsome
#         @param pos_start: start position
#         @param pos_end: end position
#         @param motif: repeat unit
#         @param motif_len: repeat unit length
#         @param repeat_times: repeat times
#         @param bam_path: bam path
#         @param reference_path: reference path
#         @prefix_len: {prefix_len} bps upstream of microsatellite to analysis
#         @suffix_len: {suffix_len} bps downstream of microsatellite to analysis
#         """
#         self.chrom = str(chrom)
#         self.pos_start = pos_start
#         self.motif = motif
#         self.motif_len = motif_len
#         self.repeat_times = repeat_times
#         self.pos_end = pos_end
#         self.bam_path = bam_path
#         self.reference_path = reference_path
#         self.ref_repeat_length = repeat_times * motif_len
#         self.mirosatellite_id = chrom + "_" + str(pos_start)
#         self.prefix_len = prefix_len
#         self.suffix_len = suffix_len
#         self.start_pre = pos_start - self.prefix_len
#         self.end_suf = pos_end + self.suffix_len
#         self.mut_start = pos_start - 1
#         self.mut_end = pos_end
#         self.repeat_length_dis = {}
#         self.query_repeat_length = self.ref_repeat_length
#         self.dis_stat = False
#         self.more_than_one_alleles = False
#         self.more_than_one_alleles_ms = False
#         self.check = True
#         self.check_stats = []
#         self.ref_str = "."
#         self.alt_str = "."
#         self.allele = 0
#         self.mut_type = None
#
#     def get_reads_alignment(self, bam_file):
#
#         alignment_list = [alignment for alignment in bam_file.fetch(self.chrom, self.pos_start - 1, self.pos_end + 1)]
#         reads_com = []
#         for alignment in alignment_list:
#             if alignment.is_unmapped or alignment.is_duplicate or alignment.is_secondary:
#                 continue
#             if alignment.reference_start > self.pos_start or alignment.reference_end < self.pos_end:
#                 continue
#             reads_com.append(alignment)
#         return reads_com
#
#     def get_dis(self):
#         """
#         Description: get the distribution of the microsateliite repeat length
#         """
#         bam_file = pysam.AlignmentFile(self.bam_path, mode="rb", reference_filename=self.reference_path)
#         repeat_length_dict = {}
#         reads = self.get_reads_alignment(bam_file)
#         bam_file.close()
#         for alignment in reads:
#             repeat_length = self.get_repeat_length(alignment)
#             if repeat_length < 0: continue
#             if repeat_length not in repeat_length_dict:
#                 repeat_length_dict[repeat_length] = 0
#             repeat_length_dict[repeat_length] += 1
#         self.repeat_length_dis = repeat_length_dict
#         self.allele = len(repeat_length_dict)
#         if self.allele < 1:
#
#             self.check = False
#             self.check_stats.append("No_read_covered")
#         else:
#             self.dis_stat = True
#             # print(repeat_length_dict)
#             self.query_repeat_length = get_max_support_index(repeat_length_dict)
#             if self.allele > 1:
#                 self.more_than_one_alleles = True
#                 self.more_than_one_alleles_ms = True
#                 self.check = False
#                 self.check_stats.append("More_alleles_in_MS")
#         if not self.check:
#             return -1
#         else:
#             return 1
#
#     def get_pileup_info(self):
#         """
#         Description: get the detail mutational information of upstream and downstream
#         """
#         mut = MutationType()
#         left_pos = self.end_suf
#         right_pos = self.start_pre
#         alt_str = []
#         fa_file = pysam.FastaFile(self.reference_path)
#
#         self.ref_str = fa_file.fetch(self.chrom, self.start_pre, self.end_suf)
#
#         if self.check:
#             bam_file = pysam.AlignmentFile(self.bam_path, mode="rb", reference_filename=self.reference_path)
#             self.ref_str = fa_file.fetch(self.chrom, self.start_pre, self.end_suf)
#             pos = self.start_pre
#             segment_pos = 0
#
#             for pot in bam_file.pileup(self.chrom,
#                                        self.start_pre,
#                                        self.end_suf,
#                                        truncate=True,
#                                        fastafile=fa_file):
#                 pot_alleles = list(set(map(lambda x: x.upper(),
#                                            pot.get_query_sequences(mark_matches=True,
#                                                                    mark_ends=False,
#                                                                    add_indels=True))))
#                 pot_alleles = collections.Counter(pot_alleles).most_common()
#                 if pot_alleles[0][0][0] in [",", "."]:
#                     alt_str.append(self.ref_str[segment_pos])
#                 elif pot_alleles[0][0][0] in ["*"]:
#                     alt_str.append("")
#                 elif pot_alleles[0][0][0] in ["A", "G", "C", "T"]:
#                     left_pos = min(left_pos, pos)
#                     right_pos = max(right_pos, pos)
#                     alt_str.append(pot_alleles[0][0][0])
#                     if pos < self.pos_start:
#                         mut.var_prefix.append([pos, pot_alleles[0][0][0], self.ref_str[segment_pos]])
#                     elif pos >= self.pos_end:
#                         mut.var_suffix.append([pos, pot_alleles[0][0][0], self.ref_str[segment_pos]])
#                     else:
#                         mut.var_ms.append([pos, pot_alleles[0][0][0], self.ref_str[segment_pos]])
#                 if len(pot_alleles[0][0]) > 1:
#                     if pot_alleles[0][0][1] == "+":
#                         left_pos = min(left_pos, pos)
#                         right_pos = max(right_pos, pos)
#                         p = re.compile("[\\+\\-][0-9]+")
#                         indel_f = p.findall(pot_alleles[0][0])[0]
#                         indel_len = int(indel_f[1:])
#                         indel_str = pot_alleles[0][0][-indel_len:]
#                         alt_str[-1] = alt_str[-1] + indel_str
#                         if pos < self.pos_start - 1:
#                             mut.var_prefix.append([pos, indel_str, ""])
#                         elif pos >= self.pos_end:
#                             mut.var_suffix.append([pos, indel_str, ""])
#                         else:
#                             mut.var_ms.append([pos, indel_str, ""])
#                     else:
#                         p = re.compile("[\\+\\-][0-9]+")
#                         indel_f = p.findall(pot_alleles[0][0])[0]
#                         indel_len = int(indel_f[1:])
#                         indel_str = pot_alleles[0][0][-indel_len:]
#                         del_end = pos + 1 + indel_len
#                         left_pos = min(left_pos, pos + 1)
#                         right_pos = max(right_pos, del_end)
#                         if del_end < self.pos_start:
#                             mut.var_prefix.append([pos + 1, '', indel_str])
#                         elif pos + 1 >= self.pos_end:
#                             mut.var_suffix.append([pos + 1, '', indel_str])
#                         else:
#
#                             if del_end < self.pos_end or pos + 1 >= self.pos_start:
#                                 mut.var_ms.append([pos + 1, '', indel_str])
#                             else:
#                                 self.check = False
#                                 self.check_stats.append("DEL_Span")
#                             # else:
#                             #     mut.var_ms.append([pos + 1, '', indel_str[0:self.pos_end - (pos + 1)]])  # TODO
#                             #     mut.var_suffix.append([self.pos_end, '', indel_str[self.pos_end - (pos + 1):]])  # TODO
#
#                 pos += 1
#                 segment_pos += 1
#             bam_file.close()
#
#             self.mut_start = min(left_pos, self.pos_start - 1)
#             self.mut_end = max(right_pos, self.pos_end)
#             self.alt_str = "".join(alt_str[self.mut_start - self.start_pre:self.mut_end - self.start_pre])
#         else:
#             pass
#         mut.comput(self.query_repeat_length - self.ref_repeat_length)
#         self.ref_str = self.ref_str[self.mut_start - self.start_pre:self.mut_end - self.start_pre]
#
#         self.mut_type = mut
#         # if len(self.alt_str) < 1:
#         #     print(self.ref_str)
#         #     print(self.alt_str)
#         #     print(self.pos_start)
#         #     print(self.mut_start - self.start_pre, self.mut_end - self.start_pre)
#         #     print(self.check_stats)
#         #     print()
#         fa_file.close()
#
#     def pos_convert_ref2read(self, ref_block: list, read_block: list, pos: int, direction="start") -> tuple:
#
#         """
#         Description: get the read position according the reference position and cigar staring
#         @param direction:  start of end
#         @param ref_block:
#         @param read_block:
#         @param start:
#         @param end:
#         @return:
#         """
#         if direction == "start":
#             block_index = 0
#             for sub_ref_block in ref_block:
#                 if sub_ref_block[0] <= pos <= sub_ref_block[1]:
#                     # print(sub_ref_block,pos)
#                     break
#                 block_index += 1
#
#             if pos == ref_block[block_index][1]:  # M|D M|I D|M
#                 read_pos = read_block[block_index][1]
#                 if ref_block[block_index][2] == 2:
#                     self.start_pre = min([self.start_pre, ref_block[block_index][0] - 1])
#                 return pos, read_pos
#
#             if ref_block[block_index][2] == 0:  # match
#                 read_pos = read_block[block_index][0] + (pos - ref_block[block_index][0])
#                 return pos, read_pos
#             elif ref_block[block_index][2] == 2:  # deletion
#                 pos = pos
#                 read_pos = read_block[block_index][1]
#                 self.start_pre = min([ref_block[block_index][0] - 1, self.start_pre])
#
#                 # pos = ref_block[block_index - 1][1] - 1
#                 # read_pos = read_block[block_index - 1][1] - 1
#                 return pos, read_pos
#             else:
#                 pos = ref_block[block_index - 1][1] - 1
#                 read_pos = read_block[block_index - 1][1] - 1
#                 return pos, read_pos
#
#         if direction == "end":
#             block_index = len(ref_block) - 1
#             for sub_ref_block in ref_block[::-1]:
#                 if sub_ref_block[0] <= pos <= sub_ref_block[1]:
#                     # print(sub_ref_block,pos)
#                     break
#                 block_index = block_index - 1
#             if pos == ref_block[block_index][0]:  # D|M I|M M|D
#                 read_pos = read_block[block_index][0]
#                 if ref_block[block_index][2] == 2:
#                     self.end_suf = max([self.end_suf, ref_block[block_index][1] + 1])
#                 return pos, read_pos
#
#             if ref_block[block_index][2] == 0:
#                 read_pos = read_block[block_index][0] + (pos - ref_block[block_index][0])
#                 return pos, read_pos
#             elif ref_block[block_index][2] == 2:
#                 pos = pos
#                 read_pos = read_block[block_index][0]
#                 self.end_suf = max([self.end_suf, ref_block[block_index][1] + 1])
#                 # pos = ref_block[block_index + 1][0] + 1
#                 # read_pos = read_block[block_index + 1][0] + 1
#                 return pos, read_pos
#             else:
#                 pos = pos + 1
#                 read_pos = read_block[block_index - 1][0] - 1
#                 return pos, read_pos
#
#     def get_repeat_length(self, alignment):
#         """
#         Description: get the repeat length according to the aligned read.
#         """
#         align_start = alignment.reference_start
#         ref_block = []
#         read_block = []
#         read_pos = 0
#         ref_pos = align_start
#         for cigartuple in alignment.cigartuples:
#             if cigartuple[0] in [0, 7, 8]:  # 0 : M : match or mishmatch ; 7: :=:match; 8:X:mismatch
#                 ref_block.append((ref_pos, ref_pos + cigartuple[1], 0))
#                 read_block.append((read_pos, read_pos + cigartuple[1], 0))
#                 read_pos += cigartuple[1]
#                 ref_pos += cigartuple[1]
#             elif cigartuple[0] in [1, 4, 5]:  # 1:I:inserion ;4:S:soft clip 5:H:hardclip
#                 ref_block.append((ref_pos, ref_pos + 0, 1))
#                 read_block.append((read_pos, read_pos + cigartuple[1], 1))
#                 read_pos += cigartuple[1]
#             elif cigartuple[0] in [2, ]:  # 2:D; 3:N: skip region of reference
#                 ref_block.append((ref_pos, ref_pos + cigartuple[1], 2))
#                 read_block.append((read_pos, read_pos, 2))
#                 ref_pos += cigartuple[1]
#             else:
#                 return -1
#
#         if self.pos_start >= ref_block[-1][1] or self.pos_start <= ref_block[0][0]:
#             return -1
#         if self.pos_end >= ref_block[-1][1] or self.pos_end <= ref_block[0][0]:
#             return -1
#
#         ref_start, read_start = self.pos_convert_ref2read(ref_block, read_block, self.pos_start, direction="start")
#         ref_end, read_end = self.pos_convert_ref2read(ref_block, read_block, self.pos_end, direction="end")
#         rpt = self.ref_repeat_length + ((read_end - read_start) - (ref_end - ref_start))
#         return rpt


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
        unit_range, repeat_range = i.split(":")
        if "-" in unit_range:
            unit_start, unit_end = tuple(map(int, unit_range.split("-")))
        else:
            unit_start = int(unit_range)
            unit_end = unit_start
        repeat_start = int(repeat_range)
        # print(unit_start,unit_end,"  ",repeat_start, repeatEnd)
        for ur in range(unit_start, unit_end + 1):
            if ur not in paras["ranges_of_repeat_times"]:
                paras["ranges_of_repeat_times"][ur] = {}
            paras["ranges_of_repeat_times"][ur]["min"] = repeat_start
        for i in args.maximum_repeat_times[0].split(";"):
            # print(i)
            unit_range, repeat_range = i.split(":")
            if "-" in unit_range:
                unit_start, unit_end = tuple(map(int, unit_range.split("-")))
            else:
                unit_start = int(unit_range)
                unit_end = unit_start
            repeat_start = int(repeat_range)
            # print(unit_start,unit_end,"  ",repeat_start, repeatEnd)
            for ur in range(unit_start, unit_end + 1):
                if ur not in paras["ranges_of_repeat_times"]:
                    paras["ranges_of_repeat_times"][ur] = {}
                paras["ranges_of_repeat_times"][ur]["max"] = repeat_start
    error_stat = False
    if os.path.exists(paras["input"]):
        logger.info("The input file is : '" + paras["input"] + "'.")
    else:
        logger.error('The input file '
                     + paras["input"] + ' is not exist, please check again')
        error_stat = True

    if os.path.isfile(paras["microsatellite"]):
        logger.info("The microsatellites file  is : " + paras["microsatellite"])
    else:
        logger.error('The microsatellites file '
                     + paras["microsatellite"] + ' is not exist, please check again')
        error_stat = True
    if os.path.isfile(paras["reference"]):
        logger.info("The reference file is : '" + paras["reference"] + "'.")
    else:
        paras["reference"] = "" if paras["reference"] == "." else paras["reference"]
        logger.error('The reference file ' + paras["reference"] + ' is not exist, please check again')
        error_stat = True
    if paras["input"][-4:] == "cram":
        paras["input_format"] = "cram"
        cramfile = pysam.AlignmentFile(paras["input"], mode="rb", reference_filename=paras["reference"])
        if not cramfile.has_index():
            logger.info("Build index for the input cram ...")
            pysam.index(paras["input"])
        cramfile.close()
    else:
        paras["input_format"] = "hap1"
        bam_file = pysam.AlignmentFile(paras["input"], mode="rb")
        if not bam_file.has_index():
            logger.info("Build index for the input bam ...")
            pysam.index(paras["input"])
        bam_file.close()
    paras["output_vcf"] = paras["output"] + ".vcf.gz"
    if not os.path.exists(paras["output_vcf"]):
        logger.info("The output is : " + paras["output_vcf"] + ".")
    else:
        logger.error(
            '[ERROR] The output ' + paras["output_vcf"] +
            ' is still exist! in case of overwrite files in this workspace, '
            'please check your script!')
        if not paras["debug"]:
            error_stat = True
    if error_stat: return False
    set_value("paras", paras)
    return True


def bm_process_one_ms_site(msDetail):
    msDetail.get_dis()
    msDetail.get_pileup_info()
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
    set_value("contigs_info", contigs_len_dict)
    outputfile.header.add_line('##INFO=<ID=chrom,Number=1,Type=String,Description="Chromosome">')
    outputfile.header.add_line('##INFO=<ID=pos,Number=1,Type=Integer,Description="Position">')
    outputfile.header.add_line('##INFO=<ID=ms_start,Number=1,Type=Integer,Description='
                               '"Position start of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=ms_end,Number=1,Type=Integer,Description='
                               '"Position end of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=motif,Number=1,Type=String,Description="Repeat unit">')
    outputfile.header.add_line('##INFO=<ID=repeat_times,Number=1,Type=Integer,Description='
                               '"Repeat times of motif in reference">')
    outputfile.header.add_line('##INFO=<ID=motif_len,Number=1,Type=Integer,Description="Repeat unit length">')
    outputfile.header.add_line('##INFO=<ID=ref_repeat_length,Number=1,Type=Integer,Description='
                               '"length of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=start_pre,Number=1,Type=Integer,Description="Start position for analysis">')
    outputfile.header.add_line('##INFO=<ID=end_suf,Number=1,Type=Integer,Description="End position for analysis">')
    outputfile.header.add_line('##INFO=<ID=mut_start,Number=1,Type=Integer,Description="Start position of mutaiton">')
    outputfile.header.add_line('##INFO=<ID=mut_end,Number=1,Type=Integer,Description="End position of mution">')
    outputfile.header.add_line('##INFO=<ID=query_repeat_length,Number=1,Type=Integer,Description='
                               '"Evaluation repeat length of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=dis_stat,Number=1,Type=String,Description='
                               '"True,the distribution is available">')
    outputfile.header.add_line('##INFO=<ID=check,Number=1,Type=String,Description='
                               '"True,the site is available for benchmark">')
    outputfile.header.add_line('##INFO=<ID=check_stats,Number=1,Type=String,Description='
                               '"Why this site is not available for benchmark">')
    outputfile.header.add_line('##INFO=<ID=allele,Number=1,Type=Integer,Description="Allele number in this site">')
    outputfile.header.add_line('##INFO=<ID=dis,Number=1,Type=String,Description='
                               'Distribution of repeat length>')
    outputfile.header.add_line('##INFO=<ID=var_type,Number=1,Type=String,Description='
                               '"Variant typeComplex, SNP, DEL, INS">')
    outputfile.header.add_line('##INFO=<ID=var_type_list,Number=1,Type=String,Description='
                               '"Variant type, Complex, SNP, DEL, INS">')
    outputfile.header.add_line('##INFO=<ID=var_detail,Number=1,Type=String,Description='
                               '"Variant Details, prefix|ms|suffix, prefix: record1:record2, record1: pos-alt.ref">')
    return outputfile


def bm_write_vcf_close(outputfile):
    outputfile.close()
    pysam.tabix_index(get_value("paras")["output_vcf"], force=True, preset="vcf")


def bm_write_vcf(outputfile, dataList):
    for msDetail in dataList:
        vcfrec = outputfile.new_record()
        # print("infoKey",vcfrec.info.keys())
        vcfrec.contig = msDetail.chrom
        # vcfrec.stop = pos + msDetail.repeat_times * len(msDetail.motif)
        vcfrec.pos = msDetail.mut_start
        vcfrec.ref = msDetail.ref_str
        vcfrec.alts = (msDetail.alt_str,) if msDetail.alt_str != "" else ("N",)
        vcfrec.id = msDetail.mirosatellite_id
        vcfrec.stop = msDetail.pos_end
        vcfrec.info["ms_start"] = msDetail.pos_start
        vcfrec.info["ms_end"] = msDetail.pos_end
        vcfrec.info["motif"] = msDetail.motif
        vcfrec.info["repeat_times"] = msDetail.repeat_times
        vcfrec.info["motif_len"] = msDetail.motif_len
        vcfrec.info["ref_repeat_length"] = msDetail.ref_repeat_length
        vcfrec.info["start_pre"] = msDetail.start_pre
        vcfrec.info["end_suf"] = msDetail.end_suf
        vcfrec.info["mut_start"] = msDetail.mut_start
        vcfrec.info["mut_end"] = msDetail.mut_end
        vcfrec.info["query_repeat_length"] = msDetail.query_repeat_length
        vcfrec.info["dis_stat"] = str(msDetail.dis_stat)
        vcfrec.info["check"] = str(msDetail.check)
        vcfrec.info["check_stats"] = "|".join(msDetail.check_stats)
        vcfrec.info["dis"] = "|".join(
            [str(key) + ":" + str(value) for key, value in msDetail.repeat_length_dis.items()])
        vcfrec.info["allele"] = msDetail.allele
        # if msDetail.check:
        vcfrec.info["var_type"] = msDetail.mut_type.var_type
        vcfrec.info["var_type_list"] = '|'.join([":".join(msDetail.mut_type.var_type_prefix),
                                                 ":".join(msDetail.mut_type.var_type_ms),
                                                 ":".join(msDetail.mut_type.var_type_suffix)])
        # print(msDetail.mut_type.var_prefix,msDetail.mut_type.var_ms,msDetail.mut_type.var_suffix)
        vcfrec.info["var_detail"] = \
            "!".join([
                "|".join([":".join([str(one[0]), one[1], one[2]]) for one in msDetail.mut_type.var_prefix]),
                "|".join([":".join([str(one[0]), one[1], one[2]]) for one in msDetail.mut_type.var_ms]),
                "|".join([":".join([str(one[0]), one[1], one[2]]) for one in msDetail.mut_type.var_suffix]),
            ])
        outputfile.write(vcfrec)


def bm_multi_run(thread, datalist):
    # pool = multiprocessing.Pool(processes=thread)
    # result_list = pool.map(bm_process_one_ms_site, datalist)
    # pool.close()
    # pool.join()
    result_list = []
    for ms in datalist:
        result_list.append(bm_process_one_ms_site(ms))

    # print("input",len(datalist))
    # print("output",len(result_list))

    return result_list


def benchmark(parase):
    if not benchmark_init(parase):
        logger.error("Benchmark init ERROR!")
        return -1
        # return if the process the arguments errors
    args = get_value("paras")
    dis_out = args["output_vcf"]
    # thread = args["threads"]
    # batch = args["batch"]
    outputfile = bm_write_vcf_init(dis_out)
    contigs_info = get_value("contigs_info")
    df_microsatellites = load_microsatellites(args)
    for contig, contig_len in contigs_info.items():
        logger.info("Processing " + contig + "...")
        this_contig_microsatellite = df_microsatellites[df_microsatellites["chr"] == contig].sort_values("pos")
        window_ms = []
        ms_num = 0
        for ms_id, info in this_contig_microsatellite.iterrows():
            ms_num += 1
            info["prefix_len"] = args["prefix_len"]
            info["suffix_len"] = args["suffix_len"]
            window_ms.append(info)
            if ms_num % (args["batch"] * args["threads"]) == 0:
                window = Window(contig, window_ms)
                window_ms = []
                window.run_window()

                # print(info["pos"])
        if len(window_ms) > 0:
            window = Window(contig, window_ms)
            del window_ms
            window.run_window()
