# =============================================================================
# Project : MShunter0.0.1
# Py Name: ScanMicrosatellites
# Author :
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
    defaultPara_bam2dis = defaultPara["bam2dis"]
    defaultPara_errEval = defaultPara["errEval"]
    defaultPara_call = defaultPara["call"]
    defaultPara_gt = defaultPara["genotype"]
    commands = ["bam2dis", "errEval", "call", "genotype"]
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
    # add arguments for bam2dis module
    parser_bam2dis = subparsers.add_parser('bam2dis', help='Get microsatellite distribution from bam file.')
    parser_bam2dis.description = 'Get microsatellite distribution from bam file.'
    parser_bam2dis.add_argument('-i', '--input', required=True, action='append', type=str,
                                help="The path of input bam file [required]")
    parser_bam2dis.add_argument('-o', '--output', required=True, action='append', type=str,
                                help="prefix of the output [required] ")
    parser_bam2dis.add_argument('-m', '--microsatellite', required=True, type=str, nargs=1,
                                help="path of the microsatellite list files [required]")
    parser_bam2dis.add_argument('-t', '--threads', type=int, nargs=1,
                                default=[defaultPara_bam2dis["threads"]],
                                help="mumber of additional threads to use [default:" +
                                     str(defaultPara_bam2dis["threads"]) + "]")
    parser_bam2dis.add_argument('-q', '--minimum_mapping_quality', type=int, nargs=1,
                                default=[defaultPara_bam2dis["minimum_mapping_quality"]],
                                help="minimum mapping quality of read [default:" +
                                     str(defaultPara_bam2dis["minimum_mapping_quality"]) + "]")
    parser_bam2dis.add_argument('-s', '--minimum_support_reads', type=int, nargs=1,
                                default=[defaultPara_bam2dis["minimum_support_reads"]],
                                help="minimum support reads of an available microsatellite [default:" +
                                     str(defaultPara_bam2dis["minimum_support_reads"]) + "]")
    parser_bam2dis.add_argument('-b', '--batch', type=int, nargs=1,
                                default=[defaultPara_bam2dis["batch"]],
                                help="batch size for one thread [default:" +
                                     str(defaultPara_bam2dis["batch"]) + "]")
    parser_bam2dis.add_argument('-d', '--debug', type=bool, nargs=1, choices=[True, False],
                                default=[defaultPara_bam2dis["debug"]],
                                help=" debug mode for developers [default:" +
                                     str(defaultPara_bam2dis["debug"]) + "]")
    parser_bam2dis.add_argument("-sep", '--separator', type=str, nargs=1, choices=["comma", "space", "tab"],
                                default=[defaultPara_bam2dis["separator"]],
                                help=' separator for microsatellites file [default:'
                                     + str(defaultPara_bam2dis["separator"]) + ']')
    parser_bam2dis.add_argument('-oh', '--only_homopolymers', type=int, nargs=1, choices=[True, False],
                                default=[defaultPara_bam2dis["allow_mismatch"]],
                                help="only analyze homopolymer regions [default:"
                                     + str(defaultPara_bam2dis["only_homopolymers"]) + "]")
    parser_bam2dis.add_argument('-am', '--allow_mismatch', type=int, nargs=1, choices=[True, False],
                                default=[defaultPara_bam2dis["allow_mismatch"]],
                                help="allow mismatch when capture microsatellite [default:"
                                     + str(defaultPara_bam2dis["allow_mismatch"]) + "]")
    parser_bam2dis.add_argument("-minr", '--minimum_repeat_times',
                                default=[defaultPara_bam2dis["minimum_repeat_times"]],
                                type=str, nargs=1,
                                help="minimum repeat times of microsatellites [default:"
                                     + defaultPara_bam2dis["minimum_repeat_times"] + "]")
    parser_bam2dis.add_argument('-maxr', '--maximum_repeat_times',
                                default=[defaultPara_bam2dis["maximum_repeat_times"]], type=str, nargs=1,
                                help="maximum repeat times of microsatellites [default:"
                                     + defaultPara_bam2dis["maximum_repeat_times"] + "]")

    commandsParser["bam2dis"] = parser_bam2dis
    ###################################################################################################################
    # add arguments for  "errEval" module
    parser_errEval = subparsers.add_parser('errEval', help='Evaluate the sequencing bias in microsatellite regions.')
    parser_errEval.description = 'Evaluate the sequencing bias in microsatellite regions.'
    parser_errEval.add_argument('-i', '--input', required=True, type=str,
                                help="The path of input dis.vcf file [required]")
    parser_errEval.add_argument('-o', '--output', required=True, type=str,
                                help="The path of output file prefix [required]")

    parser_errEval.add_argument('-s', '--minimum_support_reads', type=int, nargs=1, default=[5],
                                help="minimum support reads of an available microsatellite [default:20]")
    parser_errEval.add_argument('-l', '--max_repeat_times', type=int, nargs=1, default=[200],
                                help="maximum repeat times analysis [default:200]")
    parser_errEval.add_argument('-oh', '--only_homopolymers', type=int, nargs=1, default=[200],
                                help="only analyze homopolymer regions [default:200]")
    parser_errEval.add_argument('-t', '--threads', type=int, nargs=1,
                                default=[defaultPara_errEval["threads"]],
                                help="mumber of additional threads to use [default:" +
                                     str(defaultPara_errEval["threads"]) + "]")
    parser_errEval.add_argument('-b', '--batch', type=int, nargs=1,
                                default=[defaultPara_errEval["batch"]],
                                help="batch size for one thread [default:" +
                                     str(defaultPara_errEval["batch"]) + "]")
    commandsParser["errEval"] = parser_errEval
    ###################################################################################################################
    # add arguments for call module
    parser_call = subparsers.add_parser('call', help='Microsatellite genotyping')
    parser_call.description = 'Microsatellite genotyping.'
    parser_call.add_argument('-i', '--input', required=True, type=str,
                             help="The path of input dis.vcf file [required]")
    parser_call.add_argument('-o', '--output', required=True, type=str,
                             help="The path of output file prefix [required]")

    parser_call.add_argument('-s', '--minimum_support_reads', type=int, nargs=1, default=[5],
                             help="minimum support reads of an available microsatellite [default:20]")
    parser_call.add_argument('-l', '--max_repeat_times', type=int, nargs=1, default=[200],
                             help="maximum repeat times analysis [default:200]")
    parser_call.add_argument('-oh', '--only_homopolymers', type=int, nargs=1, default=[200],
                             help="only analyze homopolymer regions [default:200]")
    parser_call.add_argument('-t', '--threads', type=int, nargs=1,
                             default=[defaultPara_errEval["threads"]],
                             help="mumber of additional threads to use [default:" +
                                  str(defaultPara_errEval["threads"]) + "]")
    parser_call.add_argument('-b', '--batch', type=int, nargs=1,
                             default=[defaultPara_errEval["batch"]],
                             help="batch size for one thread [default:" +
                                  str(defaultPara_errEval["batch"]) + "]")
    commandsParser["call"] = parser_errEval
    ###################################################################################################################
    # add arguments for call module
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
    input_and_output.add_argument('-r', '--reference', required=False, type=str,
                                  help="Required if cram file input")
    input_and_output.add_argument("-sep", '--separator', type=str, nargs=1, choices=["comma", "space", "tab"],
                                  default=[defaultPara_gt["separator"]],
                                  help='Separator for microsatellites file [default:'
                                       + str(defaultPara_gt["separator"]) + ']')

    ##################################################################################
    # group general option
    general_option = parser_gt.add_argument_group(title="General option")
    general_option.add_argument('-d', '--debug', type=bool, nargs=1, choices=[True, False],
                                default=[defaultPara_bam2dis["debug"]],
                                help="Debug mode for developers [default:" +
                                     str(defaultPara_bam2dis["debug"]) + "]")
    general_option.add_argument('-oh', '--only_homopolymers', type=int, nargs=1, choices=[True, False],
                                default=[defaultPara_gt["allow_mismatch"]],
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
    bam2dis_option.add_argument('-am', '--allow_mismatch', type=int, nargs=1, choices=[True, False],
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

    # print(os.sys.argv)
    # print(parser.parse_args())
    # print(parser_errEval)
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
    paras["only_homopolymer"] = args.only_homopolymers[0]
    paras["threads"] = args.threads[0]
    paras["minimum_support_reads"] = args.minimum_support_reads[0]
    paras["minimum_mapping_quality"] = args.minimum_mapping_quality[0]
    paras["batch"] = args.batch[0]
    paras["debug"] = args.debug[0]
    paras["separator"] = args.separator[0]
    paras["separator"] = args.separator[0]

    paras["ranges_of_repeat_times"] = {}
    for i in args.minimum_repeat_times[0].split(";"):
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
              + paras["input"]
              + ' is not exist, please check again')
        error_stat = True

    if os.path.isfile(paras["microsatellite"]):
        print("[INFO] The microsatellites file  is : " + paras["microsatellite"])
    else:
        print('[ERROR] The microsatellites file '
              + paras["microsatellite"]
              + ' is not exist, please check again')
        error_stat = True
    if not os.path.exists(paras["output"]):
        print("[INFO] The output is : " + paras["output"] + ".")
        if not error_stat:
            print()
    else:
        print(
            '[ERROR] The output ' + paras["output"] +
            ' is still exist! in case of overwrite files in this workspace, '
            'please check your script!')
        # error_stat = True
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
    genotype_init(parase)
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
        if parase.command == "bam2dis":
            bam2dis(parase)
        elif parase.command == "errEval":
            errEval(parase)
        elif parase.command == "call":
            call(parase)
        elif parase.command == "genotype":
            genotype(parase)


if __name__ == "__main__":
    main()
