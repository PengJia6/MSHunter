import os

import pysam
import matplotlib.pyplot as plt
from src.global_dict import *


def errEval_args_init(args):
    """
    argument procress
    """
    paras = {}
    paras["input"] = args.input
    paras["output"] = args.output
    paras["threads"] = args.threads[0]
    paras["minimum_support_reads"] = args.minimum_support_reads[0]
    # paras["minimum_mapping_quality"] = args.minimum_mapping_quality[0]
    paras["batch"] = args.batch[0]
    paras["max_repeat_times"] = args.max_repeat_times[0]
    paras["only_homopolymers"] = args.only_homopolymers[0]
    paras["output"] = paras["output"] if paras["output"][-1] == "/" else paras["output"] + "/"
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
    # path_dis_parameter = path_dis_parameter if path_dis_parameter[-1] == "/" else path_dis_parameter + "/"

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
    # processLable = False
    for rec in vcffile.fetch():
        recordNum += 1
        recordInfo = rec.info
        motif = recordInfo["motif"]
        motifLen = len(motif)
        support_reads = int(recordInfo["support_reads"])
        if motifLen > 1 and onlyHomo:
            processLable = False
        else:
            processLable = True

        if processLable and (support_reads > min_support_reads):
            if motif not in File_motif:
                File_motif[motif] = pysam.VariantFile(path_dis_parameter + "tmp_motif_" + motif, 'wb',
                                                      header=vcffile.header)
            File_motif[motif].write(rec)
    motifList = []
    for motif in File_motif:
        File_motif[motif].close()
        motifList.append(motif)
    return motifList


def getHomoNormalDis(motifDis_tmp, maxRepeat):
    homList = {}

    for homo in motifDis_tmp:
        # print(motifDis_tmp[homo])
        tmp_dis = {}
        for pot in range(1, maxRepeat + 1):
            tmp_dis[pot] = 0
        readNum = len(motifDis_tmp[homo])
        for oner in motifDis_tmp[homo]:
            for pot in oner:
                tmp_dis[pot] += oner[pot]
                # tmp_dis[]
            # print(oner)
        homList[homo] = {}
        for pot in tmp_dis:
            if tmp_dis[pot] > 0:
                homList[homo][pot] = tmp_dis[pot] / readNum

    for homo in range(1, maxRepeat + 1):
        if homo not in homList:
            homList[homo] = {homo: 1}
        x, y = [], []
        for i in homList[homo]:
            x.append(i)
            y.append(homList[homo][i])
        # print(x,y)
        # # print(homo,sumdis)
        if homo>9:
            plt.plot(x,y)
            plt.scatter(x,y)
            plt.vlines(homo,0,1)
            plt.show()

            plt.close()

        #     print(homList[homo])
        # print(homo,homList[homo])
    return homList


def getOneMotifProsess(paras, motif):
    # path_dis = paras["input"]
    path_dis_parameter = paras["output"]
    min_support_reads = paras["minimum_support_reads"]
    thread = paras["threads"]
    windowSize = paras["batch"]
    dislen = paras["max_repeat_times"]
    onlyHomo = paras["only_homopolymers"]
    motifDis_tmp = {}

    vcffile = pysam.VariantFile(path_dis_parameter + "tmp_motif_" + motif, "rb")
    maxRepeat = 0
    repeatList = []
    for rec in vcffile.fetch():
        recordInfo = rec.info
        # motif = recordInfo["motif"]
        repeatTimes = int(recordInfo["repeatTimes"])
        disList = [list(map(int, i.split(":"))) for i in recordInfo["dis"].split("|")]
        thismaxRepeat = max([i[0] for i in disList])
        maxRepeat = maxRepeat if maxRepeat >= thismaxRepeat else thismaxRepeat
        support_reads = int(recordInfo["support_reads"])
        disListnormal = {}
        for i in disList:
            disListnormal[i[0]] = i[1] / support_reads
        # disListnormal=[ [i[0],i[1]/support_reads] for i in  disList]
        # print(disList)
        # print(disListnormal)
        repeatList.append(repeatTimes)

        if repeatTimes not in motifDis_tmp:
            motifDis_tmp[repeatTimes] = []
        thistmp = motifDis_tmp[repeatTimes]
        thistmp.append(disListnormal)
        motifDis_tmp[repeatTimes] = thistmp

    homList = getHomoNormalDis(motifDis_tmp, maxRepeat)
    return homList

    maxture = {}

    for first in range(1, maxRepeat):
        for second in range(1, maxRepeat):
            if first <= second:

                firstDis = homList[first]
                secondDis = homList[second]
                thismaxture = {}
                print("+++++++++++++++++++++++++++++")
                print(first, firstDis)
                print(second, secondDis)
                for rp in set([i for i in firstDis] + [j for j in secondDis]):
                    if rp not in firstDis:
                        firstDis[rp] = 0
                    if rp not in secondDis:
                        secondDis[rp] = 0
                for rp in firstDis:
                    thismaxture[rp] = round((firstDis[rp] + secondDis[rp]) / 2)
                maxture[first * 100 + second] = thismaxture
    # return {"homo": homList, "maxture": maxture}
    ########### call
    for rec in vcffile.fetch():
        recordInfo = rec.info
        # motif = recordInfo["motif"]
        repeatTimes = int(recordInfo["repeatTimes"])
        disList = [list(map(int, i.split(":"))) for i in recordInfo["dis"].split("|")]
        thismaxRepeat = max([i[0] for i in disList])
        maxRepeat = maxRepeat if maxRepeat >= thismaxRepeat else thismaxRepeat
        support_reads = int(recordInfo["support_reads"])
        disListnormal = {}
        for i in disList:
            disListnormal[i[0]] = i[1] / support_reads

    # print(thismaxture)

    #
    # with open(path_dis_parameter + "tmp_motif_" + motif) as motifDis:
    #     linenum = 0
    #     for line in motifDis:
    #         linenum += 1
    #         lineinfo = line[:-1].split(":")
    #         #                 repeatLen=int(lininfo[0].split(" ")[3])
    #         [chrom, chrompos, motif, repeatLen] = lineinfo[0].split(" ")
    #         disList = list(map(int, lineinfo[1].split(" ")))
    #         depth = sum(disList)
    #         if depth < 10:
    #             continue
    #         disList_norm = [i / depth for i in disList]
    #         start, end = 0, 0
    #         pos = 0
    #         for i in disList:
    #             pos += 1
    #             if i > 0:
    #                 start = np.max([pos - 10, 1])
    #                 break
    #         pos = len(disList)
    #         for i in disList[::-1]:
    #             pos -= 1
    #             if i > 0:
    #                 end = np.min([pos + 10, len(disList)])
    #                 break
    #
    #         # if linenum > 10:
    #         #     break
    #         distance_Euclidean = []
    #         label_list = []
    #         for first in range(start, end + 1):
    #             for second in range(start, end + 1):
    #                 if first <= second:
    #                     thisDis = maxture[first * 100 + second]
    #                     distance_Euclidean.append(np.linalg.norm(np.array(disList_norm) - np.array(thisDis)))
    #                     label_list.append(first * 100 + second)
    #         gtid = label_list[np.argmin(distance_Euclidean)]
    #         GT1 = (gtid // 100)
    #         GT2 = (gtid % 100)
    #         print(chrom, chrompos, motif, repeatLen)
    #         print(GT1, GT2)
    #         print(GT1, homoList_final[GT1])
    #         print(GT2, homoList_final[GT2])
    #
    #         print("raw", disList)
    #         print("nor", disList_norm)
    #         print("DB", maxture[gtid])
    #         print("======================================================================")
    #
    # return {"summary": summary, "homoList": homoList_final, "maxture": maxture}

    # def get_errEval(args):


def errEval(parase):
    if not errEval_args_init(parase):
        print("[Error] Parameters error!")
        return -1
    args = get_value("paras")
    # motifList = classDisbyMotif(args)
    motifList = ['T', 'A', 'G', 'C']
    # print(motifList)
    for motif in motifList:
        maxture = getOneMotifProsess(args, motif)
    #
    #     print(motif)


if __name__ == "__main__":
    ""
