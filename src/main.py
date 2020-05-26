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

# pool = multiprocessing.Pool(processes=cores)
global chromList
chromList = [str(i) for i in range(1, 23)] + ["chr" + str(i) for i in range(1, 23)] + ["chrM", "chrX", "chrY", "MT",
                                                                                       "X", "Y"]
def argumentProcress():
    """
    argument procress
    """

    parser = argparse.ArgumentParser(description='MSHunter: a Microsatellite Instability(MSI)'
                                                 ' detection tools using only tumor sequencing data!\n'
                                                 'You can test multiple sample one time in this tools')
    parser.add_argument('-i', '--input', required=True, type=str, nargs=1,
                        help="The path of input bam file [required]")
    parser.add_argument('-o', '--workspace', required=True, type=str, nargs=1, default=["./workspace"],
                        help="prefix of the output [required]")
    parser.add_argument('-m', '--microsatellite', required=True, type=str, nargs=1, default=["NA"],
                        help="path of the microsatellite list files [required]")
    parser.add_argument('-t', '--threads', type=int, nargs=1, default=[4],
                        help="mumber of additional threads to use [default:2]")
    parser.add_argument('-q', '--minimum_mapping_quality', type=int, nargs=1, default=[20],
                        help="minimum mapping quality of read [default:20]")
    parser.add_argument('-s', '--minimum_support_reads', type=int, nargs=1, default=[5],
                        help="minimum support reads of an available microsatellite [default:20]")
    parser.add_argument('-b', '--batch', type=int, nargs=1, default=[200],
                        help="batch size for one thread [default:2000]")
    args = parser.parse_args()
    ARGS = {}
    ARGS["input"] = args.input[0]
    ARGS["workspace"] = args.workspace[0] if args.workspace[0][-1] == "/" else args.workspace[0] + "/"
    ARGS["Microsatellite"] = args.microsatellite[0]
    ARGS["threads"] = args.threads[0]
    ARGS["minimum_support_reads"] = args.minimum_support_reads[0]
    ARGS["minimum_mapping_quality"] = args.minimum_mapping_quality[0]
    ARGS["batch"] = args.batch[0]
    ErrorStat = False
    if os.path.isfile(ARGS["input"]):
        print("[INFO] The input is : " + ARGS["input"])
    else:
        print('[ERROR] The input "' + ARGS["input"] + '" is not exist, please check again')
        ErrorStat = True
    if os.path.isfile(ARGS["Microsatellite"]):
        print("[INFO] The Microsatellites file  is : " + ARGS["Microsatellite"])
    else:
        print('[ERROR] The Microsatellites file "'
              + ARGS["Microsatellite"]
              + '" is not exist, please check again')
        ErrorStat = True
    if os.path.exists(ARGS["workspace"]):
        print(
            '[ERROR] The workspace is still exist! in case of overwrite files in this workspace, '
            'please give a new work space')
        # ErrorStat = True
        if ErrorStat: return False
    else:
        if ErrorStat: return False
        os.mkdir(ARGS["workspace"])
        # os.mkdir(ARGS["workspace"] + "detailInfo/")
        print("[INFO] The workspace is : " + ARGS["workspace"])
    return True

def main():
    if not argumentProcress():
        return 1

if __name__ == "__main__":
    if main() > 1:
        print("")
        print("HHH")
        # os.system("rm -r " + ARGS["workspace"])
