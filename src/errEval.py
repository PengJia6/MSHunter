import os

import numpy as np
import pysam

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
    # paras["minimum_mapping_quality"] = args.minimum_mapping_quality[0]
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

    if not os.path.exists(path_dis_parameter):
        os.mkdir(path_dis_parameter)
    else:
        if len(os.listdir(path_dis_parameter)) > 0:
            print("[Err]: Plese make sure that path", path_dis_parameter, "is not exists!")
            return
    print("[Info] Scanning the distribution file of microsatellite!")
    vcffile = pysam.VariantFile(path_dis)
    File_motif = {}
    recordNum = 0
    processLable = False
    for rec in vcffile.fetch('chr1', 100000, 200000):
        recordNum += 1
        recordInfo = rec.info
        motif = recordInfo["motif"]
        motifLen = len(motif)
        support_reads = recordInfo["support_reads"]
        if motifLen > 1 and onlyHomo:
            processLable = False
        else:
            processLable = True

        if processLable and (support_reads > min_support_reads):
            if motif not in File_motif:
                File_motif[motif] = pysam.VariantFile(path_dis_parameter + "tmp_motif_" + motif, 'w',
                                                      header=vcffile.header)
            File_motif[motif].write(rec)
    motifList = []
    for motif in File_motif:
        File_motif[motif].close()
        motifList.append(motif)
    return motifList


def getOneMotifProsess(paras, motif):
    # path_dis = paras["input"]
    path_dis_parameter = paras["output"]
    min_support_reads = paras["minimum_support_reads"]
    thread = paras["threads"]
    windowSize = paras["batch"]
    dislen = paras["max_repeat_times"]
    onlyHomo = paras["only_homopolymers"]
    motifDis_tmp = {}
    with open(path_dis_parameter + "tmp_motif_" + motif) as motifDis:
        linenum = 0
        for line in motifDis:
            linenum += 1
            lineinfo = line[:-1].split(":")
            #                 repeatLen=int(lininfo[0].split(" ")[3])
            [chrom, pos, motif, repeatLen] = lineinfo[0].split(" ")
            repeatLen = int(repeatLen)

            disList = list(map(int, lineinfo[1].split(" ")))
            if repeatLen not in motifDis_tmp:
                motifDis_tmp[repeatLen] = []
            thistmp = motifDis_tmp[repeatLen]
            thistmp.append(disList)
            motifDis_tmp[repeatLen] = thistmp
    repeatList = sorted(list(motifDis_tmp.keys()))
    summary = repeatList
    homoList = {}
    for repeatLen in repeatList:
        meanDis = np.array(motifDis_tmp[repeatLen]).mean(axis=0)
        normMeanDis = (meanDis / meanDis.sum()).round(4)
        homoList[repeatLen] = normMeanDis
    disLen = len(normMeanDis)

    for homo in range(50):
        if homo + 1 not in homoList:
            tmpthis = np.zeros(disLen)
            tmpthis[homo] = 1
            homoList[homo + 1] = tmpthis

    maxture = {}
    for first in range(50):
        for second in range(50):
            if first <= second:
                maxture[(first + 1) * 100 + (second + 1)] = list(
                    [float(i) for i in (homoList[first + 1] + homoList[second + 1]) / 2])
    homoList_final = {}
    for repeat in homoList:
        homoList_final[repeat] = list([float(i) for i in homoList[repeat]])

    with open(path_dis_parameter + "tmp_motif_" + motif) as motifDis:
        linenum = 0
        for line in motifDis:
            linenum += 1
            lineinfo = line[:-1].split(":")
            #                 repeatLen=int(lininfo[0].split(" ")[3])
            [chrom, chrompos, motif, repeatLen] = lineinfo[0].split(" ")
            disList = list(map(int, lineinfo[1].split(" ")))
            depth = sum(disList)
            if depth < 10:
                continue
            disList_norm = [i / depth for i in disList]
            start, end = 0, 0
            pos = 0
            for i in disList:
                pos += 1
                if i > 0:
                    start = np.max([pos - 10, 1])
                    break
            pos = len(disList)
            for i in disList[::-1]:
                pos -= 1
                if i > 0:
                    end = np.min([pos + 10, len(disList)])
                    break

            if linenum > 10:
                break
            distance_Euclidean = []
            label_list = []
            for first in range(start, end + 1):
                for second in range(start, end + 1):
                    if first <= second:
                        thisDis = maxture[first * 100 + second]
                        distance_Euclidean.append(np.linalg.norm(np.array(disList_norm) - np.array(thisDis)))
                        label_list.append(first * 100 + second)
            gtid = label_list[np.argmin(distance_Euclidean)]
            GT1 = (gtid // 100)
            GT2 = (gtid % 100)
            print(chrom, chrompos, motif, repeatLen)
            print(GT1, GT2)
            print(GT1, homoList_final[GT1])
            print(GT2, homoList_final[GT2])

            print("raw", disList)
            print("nor", disList_norm)
            print("DB", maxture[gtid])
            print("======================================================================")

    return {"summary": summary, "homoList": homoList_final, "maxture": maxture}


def get_errEval(args):
    motifList = classDisbyMotif(args)
    for motif in motifList:
        getOneMotifProsess(args.motif)

        print(motif)


def errEval(parase):
    if not errEval_args_init(parase):
        print("[Error] Parameters error!")
        return -1
    args = get_value("paras")
    get_errEval(args)


if __name__ == "__main__":
    ""
