import os
import numpy as np

from src.global_dict import *


def errEval_args_init(args):
    """
    argument procress
    """
    paras = {}
    paras["input"] = args.input[0]
    paras["output"] = args.output[0]
    paras["threads"] = args.threads[0]
    paras["minimum_support_reads"] = args.minimum_support_reads[0]
    paras["minimum_mapping_quality"] = args.minimum_mapping_quality[0]
    paras["batch"] = args.batch[0]
    paras["max_repeat_times"] = args.max_repeat_times[0]
    paras["only_homopolymers"] = args.only_homopolymers[0]
    ErrorStat = False

    if os.path.isfile(paras["input"]):
        print("[INFO] The", "input is : " + paras["input"])
    else:
        print('[ERROR] The', 'input "' + paras["input"] + '" is not exist, please check again')
        ErrorStat = True

    if not os.path.isfile(paras["output"]):
        print("[INFO] The ", "output is : '" + paras["output"] + "'.")
    else:
        print(
            '[ERROR] The ',
            'output "' + paras["output"] + '" is still exist! in case of overwrite files in this workspace, '
                                           'please check your script!')
        ErrorStat = True

    if ErrorStat: return False

    set_value("paras", paras)
    return True


def classDisbyMotif(paras):
    path_dis = paras["input"]
    path_dis_parameter = paras["output"]
    min_support_reads = paras["minimum_support_reads"]
    thread = paras["threads"]
    windowSize = paras["batch"]
    dislen = paras["max_repeat_times"]

    onlyHomo = paras["only_homopolymers"]
    all_dis_parameter_tmp = {}
    path_dis_parameter = path_dis_parameter if path_dis_parameter[-1] == "/" else path_dis_parameter + "/"
    motifList = []

    def _write_tmp(all_dis_parameter_tmp):
        for motif in all_dis_parameter_tmp:
            if motif not in motifList:
                motifList.append(motif)
            file = open(path_dis_parameter + "tmp_motif_" + motif, "a")
            file.write(all_dis_parameter_tmp[motif])
            file.close()

    if not os.path.exists(path_dis_parameter):
        os.mkdir(path_dis_parameter)
    else:
        if len(os.listdir(path_dis_parameter)) > 0:
            print("[Err]: Plese make sure that path", path_dis_parameter, "is not exists!")
            return
    print("[Info] Scanning the distribution file of microsatellite!")
    processLable = False
    motifList = []
    with open(path_dis) as dis:
        linenum = 0
        tmpInfo = ""
        for line in dis:
            linenum += 1
            if linenum % 2 == 1:
                tmpInfo = ""
                lineinfo = line[:-1].split(" ")
                chrom = lineinfo[0]
                pos = lineinfo[1]
                loci = "_".join([chrom, pos])
                motif = lineinfo[3].split("[")[1].split("]")[0]
                motifLen = len(motif)
                repeatLen = lineinfo[3].split("[")[0]
                if motifLen > 1 and onlyHomo:
                    processLable = False
                else:
                    processLable = True
                if processLable:
                    tmpInfo = " ".join([chrom, pos, motif, repeatLen])
            else:
                if processLable:
                    processLable = False
                    disStr = line[:-1].split(":")[1][1:-1]
                    disList = list(map(int, disStr.split(" ")))[0:dislen]
                    disStr = " ".join(list(map(str, disList[:dislen])))
                    thisDis = np.sum(disList)
                    if thisDis < min_support_reads:
                        continue
                    if motif not in all_dis_parameter_tmp:
                        all_dis_parameter_tmp[motif] = ""
                    all_dis_parameter_tmp[motif] = all_dis_parameter_tmp[motif] + str(tmpInfo) + ":" + disStr + "\n"

            if linenum % windowSize == 0:
                _write_tmp(all_dis_parameter_tmp)
                #                 print("[Info] Finished: ",linenum)
                all_dis_parameter_tmp = {}
        _write_tmp(all_dis_parameter_tmp)


def get_errEval():
    print()


def errEval(parase):
    if not errEval_args_init(parase):
        print("[Error] Parameters error!")
        return -1
    args = get_value("paras")
    get_errEval(args)


if __name__ == "__main__":
    ""
