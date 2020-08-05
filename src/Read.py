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
import pysam


class ReadInfo:
    """
    Description: Mutation
    pos_del_prefix: deletion position in prefix
    pos_del_ms:  deletion position in ms region
    pos_del_suffix: deletion position in suffix
    pos_ins_prefix: insertion position in prefix
    pos_ins_suffix: insertion position in suffix
    pos_ins_ms: insertion position in ms region
    pos_snp_prefix: mismatch position in prefix
    pos_snp_ms: mismatch position in ms regions
    pos_snp_suffix: mismatch position in suffix
    var_type_prefix: variant type in prefix
    var_type_suffix: variant type in suffix
    var_type_ms: variant type in ms region
    var_type: variant type in all region
    var_list: variant type list
    """

    def __init__(self):
        self.var_prefix = []
        self.var_suffix = []
        self.var_ms = []
        self.var_type_prefix = []
        self.var_type_suffix = []
        self.var_type_ms = []
        self.var_type = ""
        self.var_type_list = []
        self.deletion_ms = []
        self.insertion_ms = []
        self.mismatch_ms = []
        self.deletion_prefix = []
        self.insertion_prefix = []
        self.mismatch_prefix = []
        self.deletion_suffix = []
        self.insertion_suffix = []
        self.mismatch_suffix = []
        self.direction = True  # True read is forward, False: read is reversed
        self.rpl = 0  # repeat length
        self.read_name = ""
        self.read_str = ""
        self.read_list = []
        self.del_span = "None"  # all ,left,right,
        self.hap = 0  # 0: unknow , 1, 2

    def prefix_var_type(self):
        for var in self.var_prefix:
            if len(var[1]) == len(var[2]):
                self.var_type_prefix.append("SNP")
            elif len(var[1]) > len(var[2]):
                self.var_type_prefix.append("INS")
            elif len(var[1]) < len(var[2]):
                self.var_type_prefix.append("DEL")

    def suffix_var_type(self):
        for var in self.var_suffix:
            if len(var[1]) == len(var[2]):
                self.var_type_suffix.append("SNP")
            elif len(var[1]) > len(var[2]):
                self.var_type_suffix.append("INS")
            elif len(var[1]) < len(var[2]):
                self.var_type_suffix.append("DEL")

    def ms_var_type(self, offset):

        if offset < 0:
            self.var_type_ms.append("DEL")
        if offset > 0:
            self.var_type_ms.append("INS")

        for var in self.var_ms:
            if len(var[1]) == len(var[2]):
                self.var_type_ms.append("SNP")

    def comput(self, offset):
        self.ms_var_type(offset)
        self.prefix_var_type()
        self.suffix_var_type()
        self.var_type_list = self.var_type_prefix + self.var_type_ms + self.var_type_suffix
        var_num = len(self.var_type_list)
        if var_num > 1:
            self.var_type = "Complex"
        elif var_num == 1:
            self.var_type = self.var_type_list[0]
        else:
            self.var_type = "None"

    def show_info(self):
        # if len()
        print("read_name", self.read_name)
        print("read_str", self.read_str)
        print("deletion", self.deletion)
        print("insertion", self.insertion)
        print("mismatch", self.mismatch)
        print("direction", self.direction)
        print("rpl", self.rpl)
        print()


class Read:
    def __init__(self, read_id, alignment, reference, chrom):
        self.chrom = chrom
        # self.read_name = alignment.query_name
        self.align_start = alignment.reference_start
        self.align_end = alignment.reference_end
        self.this_read_str = alignment.query_sequence
        self.read_id = read_id
        self.reference = reference
        self.support_microsatellites = []
        # self.alignment = alignment
        self.microsatellite = {}
        self.direction = False if alignment.is_reverse else True
        self.hap = int(alignment.get_tag("HP")) if alignment.has_tag("HP") else 0
        self.cigartuples = alignment.cigartuples


    # def get_microsatellite_detail(self, ms_info):
    #     self. = ms_info


    def get_read_str(self):
        pass
        # TODO:  To be optimized

        this_ref_str = pysam.FastaFile(self.reference).fetch(self.chrom, start=self.align_start, end=self.align_end)

        sub_read_str = []
        sub_ref_str = []
        ref_block = []
        read_block = []
        read_mut_info = ReadInfo()
        read_mut_info.direction = self.direction
        read_mut_info.hap = self.hap
        read_pos = 0
        ref_pos = self.align_start
        ref_pos_str = 0
        for cigartuple in self.cigartuples:
            if cigartuple[0] in [0, 7, 8]:  # 0 : M : match or mishmatch ; 7: :=:match; 8:X:mismatch
                ref_block.append((ref_pos, ref_pos + cigartuple[1], 0))
                read_block.append((read_pos, read_pos + cigartuple[1], 0))
                match_ref = list(this_ref_str[ref_pos_str:ref_pos_str + cigartuple[1]])
                match_read = list(self.this_read_str[read_pos:read_pos + cigartuple[1]])
                sub_ref_str.extend(match_ref)
                sub_read_str.extend(match_read)
                # if ref_pos + cigartuple[1] < self.start_pre or ref_pos > self.end_suf:
                #     pass
                # else:
                #     pos = ref_pos - 1
                #     for ref_pot, read_pot in zip(match_ref, match_read):
                #         pos += 1
                #         if pos < self.start_pre or pos > self.end_suf: continue
                #         if read_pot == "N" or ref_pot == "N": continue
                #         if read_pot == ref_pot: continue
                #         if pos < self.pos_start:
                #             read_mut_info.mismatch_prefix.append([pos, ref_pot, read_pot])
                #         elif pos >= self.pos_end:
                #             read_mut_info.mismatch_suffix.append([pos, ref_pot, read_pot])
                #         else:
                #             read_mut_info.mismatch_ms.append([pos, ref_pot, read_pot])
                read_pos += cigartuple[1]
                ref_pos += cigartuple[1]
                ref_pos_str += cigartuple[1]
            elif cigartuple[0] in [1, 4, 5]:  # 1:I:inserion ;4:S:soft clip 5:H:hardclip
                ref_block.append((ref_pos, ref_pos + 0, 1))
                read_block.append((read_pos, read_pos + cigartuple[1], 1))
                if cigartuple[0] == 1:
                    # print(read_pos)
                    # print(alignment.cigartuples)
                    # print(sub_read_str[-1])
                    # print(self.this_read_str[read_pos:read_pos + cigartuple[1]])
                    if len(sub_read_str)<1:
                        print(self.cigartuples)
                    sub_read_str[-1] += self.this_read_str[read_pos:read_pos + cigartuple[1]]
                    # pos = ref_pos - 1
                    # if ref_pos + cigartuple[1] < self.start_pre or ref_pos > self.end_suf:
                    #     pass
                    # else:
                    #     if pos < self.start_pre or pos > self.end_suf: continue
                    #     if pos < self.pos_start - 1:
                    #         read_mut_info.insertion_prefix.append([pos, sub_ref_str[-1], sub_read_str[-1]])
                    #     elif pos >= self.pos_end:
                    #         read_mut_info.insertion_suffix.append([pos, sub_ref_str[-1], sub_read_str[-1]])
                    #     else:
                    #         read_mut_info.insertion_ms.append([pos, sub_ref_str[-1], sub_read_str[-1]])
                read_pos += cigartuple[1]
            elif cigartuple[0] in [2, ]:  # 2:D; 3:N: skip region of reference
                ref_block.append((ref_pos, ref_pos + cigartuple[1], 2))
                read_block.append((read_pos, read_pos, 2))
                delete_str = this_ref_str[ref_pos_str:ref_pos_str + cigartuple[1]]
                sub_ref_str.extend(list(delete_str))
                sub_read_str.extend([""] * cigartuple[1])
                # pos = ref_pos
                # end_pos = ref_pos + cigartuple[1]
                # if ref_pos < self.start_pre:  # ref_pos + cigartuple[1]
                #     if end_pos < self.start_pre:
                #         pass
                #     elif end_pos < self.end_suf:
                #         read_mut_info.del_span = "Left"
                #     else:
                #         read_mut_info.del_span = "All"
                # elif ref_pos > self.end_suf:
                #     pass
                # else:
                #     if end_pos >= self.end_suf:
                #         read_mut_info.del_span = "Right"
                #     else:
                #         # if pos < self.start_pre or pos > self.end_suf: continue
                #         if pos < self.pos_start:
                #             read_mut_info.deletion_prefix.append([pos, delete_str, ""])
                #         elif pos >= self.pos_end:
                #             read_mut_info.deletion_suffix.append([pos, delete_str, ""])
                #         else:
                #             read_mut_info.deletion_ms.append([pos, delete_str, ""])
                    # if
                ref_pos += cigartuple[1]
                ref_pos_str += cigartuple[1]
            else:
                return -1
        # this_repeat_length = len(
        #     "".join(sub_read_str[self.pos_start - 1 - self.align_start:self.pos_end - self.align_start - 1]))
        # read_mut_info.rpl = this_repeat_length
        # read_mut_info.read_name = alignment.query_name + "_" + ("1" if alignment.is_read1 else "2")
        # read_mut_info.read_str = "".join(sub_read_str[self.start_pre - align_start:self.end_suf - align_start])
        # read_mut_info.read_list = sub_read_str[self.start_pre - align_start:self.end_suf - align_start]
        return read_mut_info

# def update_ms_id(self, ms_id):
# if ms_id
