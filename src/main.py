# =============================================================================
# Project : MShunter0.0.1
# Py Name: ScanMicrosatellites
# Author : main.py
# Date : 20-05-25
# Email : pengjia@stu.xjtu.edu.cn
# Description : ''
# =============================================================================

import argparse
import os
import sys

curpath = os.path.abspath(os.path.dirname(sys.argv[0]))
sys.path.append(os.path.dirname(curpath))
from src.bam2dis import *
from src.call import *
from src.errEval import *
from src.global_dict import *

print(" ".join(sys.argv))


def args_process():
    """
    argument procress
    """
    defaultPara = get_value("default")
    defaultPara_gt = defaultPara["genotype"]
    commands = ["genotype"]
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
                                  help="Required if cram file input")
    input_and_output.add_argument('-tech', '--technology', type=str, nargs=1, choices=["ccs", "clr", "ont", "ilm"],
                                  default=[defaultPara_gt["tech"]],
                                  help='Sequencing technology [default:'
                                       + str(defaultPara_gt["tech"]) + ']')
    input_and_output.add_argument('-hap', '--haplotype_bam', type=bool, nargs=1, choices=[True, False],
                                  default=[defaultPara_gt["hap"]],
                                  help=" Input bam file with haplotype tags [default:"
                                       + str(defaultPara_gt["hap"]) + "]")
    input_and_output.add_argument("-sep", '--separator', type=str, nargs=1, choices=["comma", "space", "tab"],
                                  default=[defaultPara_gt["separator"]],
                                  help='Separator for microsatellites file [default:'
                                       + str(defaultPara_gt["separator"]) + ']')
    ##################################################################################
    # group read realignment
    general_realign = parser_gt.add_argument_group(title="Read realignment")
    general_realign.add_argument('-pl', '--prefix_len', type=int, nargs=1,
                                 default=[defaultPara_gt["prefix_len"]],
                                 help="Debug mode for developers [default:" +
                                      str(defaultPara_gt["prefix_len"]) + "]")
    general_realign.add_argument('-sl', '--suffix_len', type=int, nargs=1,
                                 default=[defaultPara_gt["suffix_len"]],
                                 help="Debug mode for developers [default:" +
                                      str(defaultPara_gt["suffix_len"]) + "]")
    general_realign.add_argument('-ks', '--kmer_size', type=int, nargs=1,
                                 default=[defaultPara_gt["kmer_size"]],
                                 help="Debug mode for developers [default:" +
                                      str(defaultPara_gt["kmer_size"]) + "]")

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
    ##################################################################################
    # group for bam2dis
    bam2dis_option = parser_gt.add_argument_group(title="Option for bam2dis")
    bam2dis_option.add_argument('-q', '--minimum_mapping_quality', type=int, nargs=1,
                                default=[defaultPara_gt["minimum_mapping_quality"]],
                                help="minimum mapping quality of read [default:" +
                                     str(defaultPara_gt["minimum_mapping_quality"]) + "]")
    bam2dis_option.add_argument('-s', '--minimum_support_reads', type=int, nargs=1,
                                default=[defaultPara_gt["minimum_support_reads"]],
                                help="minimum support reads of an available microsatellite [default:" +
                                     str(defaultPara_gt["minimum_support_reads"]) + "]")
    bam2dis_option.add_argument('-am', '--allow_mismatch', type=bool, nargs=1, choices=[True, False],
                                default=[defaultPara_gt["allow_mismatch"]],
                                help="allow mismatch when capture microsatellite [default:"
                                     + str(defaultPara_gt["allow_mismatch"]) + "]")

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


def genotype_init(args):
    """
    argument procress
    """
    paras = {}
    paras["input"] = args.input[0]
    paras["output"] = args.output[0]
    paras["microsatellite"] = args.microsatellite[0]
    paras["reference"] = args.reference[0]
    paras["separator"] = args.separator[0]
    paras["tech"] = args.technology[0]
    paras["hap"] = args.haplotype_bam[0]
    paras["prefix_len"] = args.prefix_len[0]
    paras["suffix_len"] = args.suffix_len[0]
    paras["debug"] = args.debug[0]
    paras["only_homopolymer"] = args.only_homopolymers[0]
    paras["minimum_support_reads"] = args.minimum_support_reads[0]
    paras["minimum_mapping_quality"] = args.minimum_mapping_quality[0]
    paras["allow_mismatch"] = args.allow_mismatch[0]
    paras["threads"] = args.threads[0]
    paras["batch"] = args.batch[0]
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
        paras["input_format"] = "bam"
        bamfile = pysam.AlignmentFile(paras["input"], mode="rb")
        if not bamfile.has_index():
            print("[INFO] Build index for the input bam ...")
            pysam.index(paras["input"])
        bamfile.close()

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
    paras["output_dis"] = paras["output"] + case + ".dis.bcf"
    paras["output_tmp"] = paras["output"] + case + "_tmp"
    if not os.path.exists(paras["output_tmp"]):
        os.makedirs(paras["output_tmp"])
    paras["output_model"] = paras["output"] + case + ".model"
    paras["output_call"] = paras["output"] + case + ".bcf"
    set_value("case", case)
    set_value("paras", paras)
    return True


def genotype(parase):
    if genotype_init(parase):
        bam2dis()
        errEval()
        call()


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


if __name__ == "__main__":
    main()
