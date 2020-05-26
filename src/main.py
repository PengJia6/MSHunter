# =============================================================================
# Project : MShunter0.0.1
# Py Name: ScanMicrosatellites
# Author :
# Date : 20-05-25
# Email : pengjia@stu.xjtu.edu.cn
# Description : ''
# =============================================================================
import argparse
from src.bam2dis import *
from src.call import *
from src.errEval import *
from src.global_dict import *
import sys
print(sys.argv)


def args_process():
    """
    argument procress
    """
    commands = ["bam2dis", "errEval", "call"]
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
                                help="The path of input bam file [required] （allow to specify multiple times）")
    parser_bam2dis.add_argument('-o', '--output', required=True, action='append', type=str,
                                help="prefix of the output [required] （allow to specify multiple times）")
    parser_bam2dis.add_argument('-m', '--microsatellite', required=True, type=str, nargs=1, default=["NA"],
                                help="path of the microsatellite list files [required]")
    parser_bam2dis.add_argument('-t', '--threads', type=int, nargs=1, default=[4],
                                help="mumber of additional threads to use [default:2]")
    parser_bam2dis.add_argument('-q', '--minimum_mapping_quality', type=int, nargs=1, default=[20],
                                help="minimum mapping quality of read [default:20]")
    parser_bam2dis.add_argument('-s', '--minimum_support_reads', type=int, nargs=1, default=[5],
                                help="minimum support reads of an available microsatellite [default:20]")
    parser_bam2dis.add_argument('-b', '--batch', type=int, nargs=1, default=[200],
                                help="batch size for one thread [default:2000]")

    ###################################################################################################################
    # add arguments for  "errEval" module

    ###################################################################################################################
    # add arguments for call module
    if os.sys.argv[1] in ["-h", "--help"]:
        parser.print_help()
        return False
    if os.sys.argv[1] in ["-V", "-v", "--version"]:
        parser.parse_args()
        return False
    if len(os.sys.argv) == 1 or os.sys.argv[1] not in commands:
        print("[Error] Command Error!")
        print("[Tips] Please input correct command such as " + ", ".join(commands) + "!")
        parser.print_help()
        parser.parse_args()
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
