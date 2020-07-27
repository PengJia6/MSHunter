#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : main.py
# Author : Peng Jia
# Date   : 2020.07.13
# Email  : pengjia@stu.xjtu.edu.cn
# Description: main function and arguments processing
=============================================================================="""
import argparse
import os
import sys

curpath = os.path.abspath(os.path.dirname(sys.argv[0]))
sys.path.append(os.path.dirname(curpath))

from src.benchmark import *
from src.genotype import *
from src.benchmark_merge import benchmark_merge

print(" ".join(sys.argv))


def args_process():
    """
    argument procress
    """
    defaultPara = get_value("default")

    commands = []
    commandsParser = {}
    parser = argparse.ArgumentParser(description='mstools: Microsatellite genotyping toolbox.'
                                     # + ".show help of subcommand with '"
                                     # + get_value("tools_name") + " <subcommand> -h'"
                                     )
    parser.usage = get_value("tools_name") + " <command> [options]"
    parser.add_argument('-V', '--version', action='version',
                        version=get_value("tools_name") + get_value("tools_version"))
    subparsers = parser.add_subparsers(title="command", metavar="", dest='command')

    ###################################################################################################################
    # add arguments for genotype module
    parser_gt = subparsers.add_parser('genotype', help='Microsatellite genotyping')
    parser_gt.description = 'Microsatellite genotype.'
    commands.append("genotype")
    defaultPara_gt = defaultPara["genotype"]
    ##################################################################################
    # group input and output
    input_and_output = parser_gt.add_argument_group(title="Input and output")
    input_and_output.add_argument('-i', '--input', required=True, type=str, nargs=1,
                                  help="The path of input bam/cram file [required]")
    input_and_output.add_argument('-m', '--microsatellite', required=True, type=str, nargs=1,
                                  help="The path of the microsatellite regions [required]")
    input_and_output.add_argument('-o', '--output', required=True, type=str, nargs=1,
                                  help="The path of output file prefix [required]")
    input_and_output.add_argument('-r', '--reference', required=True, type=str, nargs=1,
                                  help="The path of reference file [required]")
    input_and_output.add_argument('-tech', '--technology', type=str, nargs=1, choices=["ccs", "clr", "ont", "ilm"],
                                  required=True,
                                  help='Sequencing technology [required]')
    # input_and_output.add_argument('-tech', '--technology', type=str, nargs=1, choices=["ccs", "clr", "ont", "ilm"],
    #                               default=[defaultPara_gt["tech"]],
    #                               help='Sequencing technology [default:'
    #                                    + str(defaultPara_gt["tech"]) + ']')
    input_and_output.add_argument('-hap', '--haplotype_bam', type=bool, nargs=1, choices=[True, False],
                                  default=[defaultPara_gt["hap"]],
                                  help=" Input bam file with haplotype tags [default:"
                                       + str(defaultPara_gt["hap"]) + "]")
    input_and_output.add_argument("-sep", '--separator', type=str, nargs=1, choices=["comma", "space", "tab"],
                                  default=[defaultPara_gt["separator"]],
                                  help='Separator for microsatellites file [default:'
                                       + str(defaultPara_gt["separator"]) + ']')
    ##################################################################################
    # group Analysis regions
    general_realign = parser_gt.add_argument_group(title="Analysis regions")
    general_realign.add_argument('-pl', '--prefix_len', type=int, nargs=1,
                                 default=[defaultPara_gt["prefix_len"]],
                                 help="Debug mode for developers [default:" +
                                      str(defaultPara_gt["prefix_len"]) + "]")
    general_realign.add_argument('-sl', '--suffix_len', type=int, nargs=1,
                                 default=[defaultPara_gt["suffix_len"]],
                                 help="Debug mode for developers [default:" +
                                      str(defaultPara_gt["suffix_len"]) + "]")
    # general_realign.add_argument('-ks', '--kmer_size', type=int, nargs=1,
    #                              default=[defaultPara_gt["kmer_size"]],
    #                              help="Debug mode for developers [default:" +
    #                                   str(defaultPara_gt["kmer_size"]) + "]")

    ##################################################################################
    # group general option
    general_option = parser_gt.add_argument_group(title="General option")
    general_option.add_argument('-d', '--debug', type=bool, nargs=1, choices=[True, False],
                                default=[defaultPara_gt["debug"]],
                                help="Debug mode for developers [default:" +
                                     str(defaultPara_gt["debug"]) + "]")
    general_option.add_argument('-oh', '--only_homopolymers', type=int, nargs=1, choices=[True, False],
                                default=[defaultPara_gt["only_homopolymers"]],
                                help="Only analyze homopolymer regions [default:"
                                     + str(defaultPara_gt["only_homopolymers"]) + "]")
    general_option.add_argument("-minr", '--minimum_repeat_times',
                                default=[defaultPara_gt["minimum_repeat_times"]],
                                type=str, nargs=1,
                                help="Minimum repeat times of microsatellites [default:"
                                     + defaultPara_gt["minimum_repeat_times"] + "]")
    general_option.add_argument('-maxr', '--maximum_repeat_times',
                                default=[defaultPara_gt["maximum_repeat_times"]], type=str, nargs=1,
                                help="Maximum repeat times of microsatellites [default:"
                                     + defaultPara_gt["maximum_repeat_times"] + "]")
    general_option.add_argument('-minh', '--minimum_phasing_reads',
                                default=[defaultPara_gt["minimum_phasing_reads"]], type=str, nargs=1,
                                help="Minimum reads for each haplotype reporting [default:"
                                     + str(defaultPara_gt["minimum_phasing_reads"]) + "]")
    general_option.add_argument('-q', '--minimum_mapping_quality', type=int, nargs=1,
                                default=[defaultPara_gt["minimum_mapping_quality"]],
                                help="minimum mapping quality of read [default:" +
                                     str(defaultPara_gt["minimum_mapping_quality"]) + "]")
    general_option.add_argument('-s', '--minimum_support_reads', type=int, nargs=1,
                                default=[defaultPara_gt["minimum_support_reads"]],
                                help="minimum support reads of an available microsatellite [default:" +
                                     str(defaultPara_gt["minimum_support_reads"]) + "]")

    ##################################################################################
    # group for bam2dis
    # bam2dis_option = parser_gt.add_argument_group(title="Option for bam2dis")

    # bam2dis_option.add_argument('-am', '--allow_mismatch', type=bool, nargs=1, choices=[True, False],
    #                             default=[defaultPara_gt["allow_mismatch"]],
    #                             help="allow mismatch when capture microsatellite [default:"
    #                                  + str(defaultPara_gt["allow_mismatch"]) + "]")

    ##################################################################################
    # group for multiple_thread

    multiple_thread = parser_gt.add_argument_group(title="Multiple thread")
    multiple_thread.add_argument('-t', '--threads', type=int, nargs=1,
                                 default=[defaultPara_gt["threads"]],
                                 help="The number of  threads to use [default:" +
                                      str(defaultPara_gt["threads"]) + "]")
    multiple_thread.add_argument('-b', '--batch', type=int, nargs=1,
                                 default=[defaultPara_gt["batch"]],
                                 help="The number of microsatellite one thread process [default:" +
                                      str(defaultPara_gt["batch"]) + "]")
    commandsParser["genotype"] = parser_gt




    ###################################################################################################################
    # add arguments for benchmark module
    parser_bm = subparsers.add_parser('benchmark', help='Microsatellite genotyping benchmark using phased assemblies')
    parser_bm.description = 'Microsatellite associate mutaion benchmark.'
    commands.append('benchmark')
    defaultPara_bm = defaultPara["benchmark"]
    ##################################################################################
    # group input and output
    bm_input_and_output = parser_bm.add_argument_group(title="Input and output")
    bm_input_and_output.add_argument('-i', '--input', required=True, type=str, nargs=1,
                                     help="Aligned contig of phased assembly, bam/cram [required]")
    bm_input_and_output.add_argument('-m', '--microsatellite', required=True, type=str, nargs=1,
                                     help="The path of the microsatellite regions [required]")
    bm_input_and_output.add_argument('-o', '--output', required=True, type=str, nargs=1,
                                     help="The path of output file prefix [required]")
    bm_input_and_output.add_argument('-r', '--reference', required=True, type=str, nargs=1,
                                     help="The path of reference file fa/fasta[required]")
    bm_input_and_output.add_argument("-sep", '--separator', type=str, nargs=1, choices=["comma", "space", "tab"],
                                     default=[defaultPara_bm["separator"]],
                                     help='Separator for microsatellites file [default:'
                                          + str(defaultPara_bm["separator"]) + ']')
    ##################################################################################
    # group read realignment
    # bm_general_realign = parser_bm.add_argument_group(title="Read realignment")
    #
    # bm_general_realign.add_argument('-ks', '--kmer_size', type=int, nargs=1,
    #                                 default=[defaultPara_bm["kmer_size"]],
    #                                 help="Debug mode for developers [default:" +
    #                                      str(defaultPara_bm["kmer_size"]) + "]")

    ##################################################################################
    # group general option
    bm_general_option = parser_bm.add_argument_group(title="General option")
    bm_general_option.add_argument('-d', '--debug', type=bool, nargs=1, choices=[True, False],
                                   default=[defaultPara_bm["debug"]],
                                   help="Debug mode for developers [default:" +
                                        str(defaultPara_bm["debug"]) + "]")
    bm_general_option.add_argument('-oh', '--only_homopolymers', type=int, nargs=1, choices=[True, False],
                                   default=[defaultPara_bm["only_homopolymers"]],
                                   help="Only analyze homopolymer regions [default:"
                                        + str(defaultPara_bm["only_homopolymers"]) + "]")
    bm_general_option.add_argument('-pl', '--prefix_len', type=int, nargs=1,
                                   default=[defaultPara_bm["prefix_len"]],
                                   help=" {prefix_len} bps upstream of microsatellite to analysis [default:" +
                                        str(defaultPara_bm["prefix_len"]) + "]")
    bm_general_option.add_argument('-sl', '--suffix_len', type=int, nargs=1,
                                   default=[defaultPara_bm["suffix_len"]],
                                   help=" {suffix_len} bps downstream of microsatellite to analysis [default:" +
                                        str(defaultPara_bm["suffix_len"]) + "]")
    bm_general_option.add_argument('-om', '--only_microsatellites', type=int, nargs=1, choices=[True, False],
                                   default=[defaultPara_bm["only_microsatellites"]],
                                   help="True, only detect variants in microsatelite microsatellite region;"
                                        " False, also detect upstream and downstream variants [default:"
                                        + str(defaultPara_bm["only_microsatellites"]) + "]")
    bm_general_option.add_argument("-minr", '--minimum_repeat_times',
                                   default=[defaultPara_bm["minimum_repeat_times"]],
                                   type=str, nargs=1,
                                   help="Minimum repeat times of microsatellites [default:"
                                        + defaultPara_bm["minimum_repeat_times"] + "]")
    bm_general_option.add_argument('-maxr', '--maximum_repeat_times',
                                   default=[defaultPara_bm["maximum_repeat_times"]], type=str, nargs=1,
                                   help="Maximum repeat times of microsatellites [default:"
                                        + defaultPara_bm["maximum_repeat_times"] + "]")
    bm_general_option.add_argument('-minh', '--minimum_phasing_reads',
                                   default=[defaultPara_bm["minimum_phasing_reads"]], type=str, nargs=1,
                                   help="Minimum reads for each haplotype reporting [default:"
                                        + str(defaultPara_bm["minimum_phasing_reads"]) + "]")

    ##################################################################################
    # group for bam2dis
    bm_bam2dis_option = parser_bm.add_argument_group(title="Option for bam2dis")
    bm_bam2dis_option.add_argument('-q', '--minimum_mapping_quality', type=int, nargs=1,
                                   default=[defaultPara_bm["minimum_mapping_quality"]],
                                   help="minimum mapping quality of read [default:" +
                                        str(defaultPara_bm["minimum_mapping_quality"]) + "]")
    bm_bam2dis_option.add_argument('-s', '--minimum_support_reads', type=int, nargs=1,
                                   default=[defaultPara_bm["minimum_support_reads"]],
                                   help="minimum support reads of an available microsatellite [default:" +
                                        str(defaultPara_bm["minimum_support_reads"]) + "]")
    bm_bam2dis_option.add_argument('-am', '--allow_mismatch', type=bool, nargs=1, choices=[True, False],
                                   default=[defaultPara_bm["allow_mismatch"]],
                                   help="allow mismatch when capture microsatellite [default:"
                                        + str(defaultPara_bm["allow_mismatch"]) + "]")

    ##################################################################################
    # group for multiple_thread

    multiple_thread = parser_bm.add_argument_group(title="Multiple thread")
    multiple_thread.add_argument('-t', '--threads', type=int, nargs=1,
                                 default=[defaultPara_bm["threads"]],
                                 help="The number of  threads to use [default:" +
                                      str(defaultPara_bm["threads"]) + "]")
    multiple_thread.add_argument('-b', '--batch', type=int, nargs=1,
                                 default=[defaultPara_bm["batch"]],
                                 help="The number of microsatellite one thread process [default:" +
                                      str(defaultPara_bm["batch"]) + "]")

    commandsParser["benchmark"] = parser_bm

    # bm_input_and_output.add_argument('-r', '--input', required=True, type=str, nargs=1,
    #                                  help="The path of input bam/cram file [required]")
    # bm_input_and_output.add_argument('-m', '--input', required=True, type=str, nargs=1,
    #                                  help="The path of input bam/cram file [required]")

    ###################################################################################################################
    # add arguments for benchmark module
    parser_bmm = subparsers.add_parser('benchmark_merge', help='Merge two haplotype result.')
    parser_bmm.description = 'Merge two haplotype microsatellite calling result.'
    commands.append('benchmark_merge')
    default_para_bmm = defaultPara["benchmark_merge"]

    ##################################################################################
    # group input and output
    bmm_input_and_output = parser_bmm.add_argument_group(title="Input and output")
    bmm_input_and_output.add_argument('-1', '--hap1', required=True, type=str, nargs=1,
                                      help="microsatellite calling result of haplotype 1, *vcf.gz [required]")
    bmm_input_and_output.add_argument('-2', '--hap2', required=True, type=str, nargs=1,
                                      help="microsatellite calling result of haplotype 2, *vcf.gz [required]")
    bmm_input_and_output.add_argument('-o', '--output', required=True, type=str, nargs=1,
                                      help="The path of output file prefix [required]")
    bmm_input_and_output.add_argument('-s', '--sample', required=False, type=str, nargs=1,
                                      default=[default_para_bmm["sample"]],
                                      help="Sample name[default:" + default_para_bmm["sample"] + "]")
    bmm_general_option = parser_bmm.add_argument_group(title="General option")
    bmm_general_option.add_argument('-d', '--debug', type=bool, nargs=1, choices=[True, False],
                                   default=[default_para_bmm["debug"]],
                                   help="Debug mode for developers [default:" +
                                        str(default_para_bmm["debug"]) + "]")
    if len(os.sys.argv) < 2:
        parser.print_help()
        return False

    if os.sys.argv[1] in ["-h", "--help", "-?"]:
        parser.print_help()
        return False
    if os.sys.argv[1] in ["-V", "-v", "--version"]:
        # parser.print_help()
        parser.parse_args()
        return False
    if os.sys.argv[1] not in commands:
        print("[Error] Command Error! ", os.sys.argv[1], "is not the available command")
        print("[Tips] Please input correct command such as " + ", ".join(commands) + "!")
        parser.print_help()
        # parser.parse_args()
        return False
    if len(os.sys.argv) == 2 and (os.sys.argv[1] in commands):
        commandsParser[os.sys.argv[1]].print_help()
        return False
    return parser


def main():
    """
    Main function.
    :return:
    """
    global_init()
    arg = args_process()

    if arg:
        parase = arg.parse_args()
        if parase.command == "genotype":
            genotype(parase)
            # genotype_ngs(parase)
        if parase.command == "benchmark":
            benchmark(parase)
        if parase.command == "benchmark_merge":
            benchmark_merge(parase)
        # if parase.command == "ngs":
        #     # genotype(parase)
        #     genotype_ngs(parase)


if __name__ == "__main__":
    main()
