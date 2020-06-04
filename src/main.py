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
    commands = ["bam2dis", "errEval", "call"]
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
    parser_bam2dis.add_argument("-minr",'--minimum_repeat_times', default=[defaultPara_bam2dis["minimum_repeat_times"]],
                                type=str, nargs=1,
                                help="minimum repeat times of microsatellites [default:"
                                     + defaultPara_bam2dis["minimum_repeat_times"] + "]")
    parser_bam2dis.add_argument('-maxr','--maximum_repeat_times', default=[defaultPara_bam2dis["maximum_repeat_times"]], type=str, nargs=1,
                                help="maximum repeat times of microsatellites [default:"
                                     +defaultPara_bam2dis["maximum_repeat_times"]+"]")

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
    parser_errEval.description = 'Microsatellite genotyping.'
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

    # print(os.sys.argv)
    # print(parser.parse_args())
    # print(parser_errEval)
    ###################################################################################################################
    # add arguments for  "list2bed" module
    # convert MSIsensor list file to bed format

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
        if parase.command == "bam2dis":
            bam2dis(parase)
        elif parase.command == "errEval":
            errEval(parase)
        elif parase.command == "call":
            call(parase)


if __name__ == "__main__":
    main()
