import os

import pysam
import yaml
import matplotlib.pyplot as plt
from src.global_dict import *
from src.units import *




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
                File_motif[motif] = pysam.VariantFile(path_dis_parameter + "tmp_motif_" + motif + ".bcf", 'wb',
                                                      header=vcffile.header)
            File_motif[motif].write(rec)
    motifList = []
    for motif in File_motif:
        File_motif[motif].close()
        motifList.append(motif)
    return motifList


def write_vcf_init_call(outputpath, inputpath):
    inputfile = pysam.VariantFile(inputpath, "rb")
    outputfile = pysam.VariantFile(outputpath, "wb", header=inputfile.header)
    outputfile.header.add_line(
        '##INFO=<ID=FirstAlleles,Number=1,Type=String,Description="The first allele type of this point">')
    outputfile.header.add_line(
        '##INFO=<ID=SecondAlleles,Number=1,Type=String,Description="The second allele type of this point">')
    outputfile.header.add_line('##INFO=<ID=CallQuality,Number=1,Type=String,Description="Genotype quality">')
    outputfile.header.add_line('##INFO=<ID=Mutation,Number=1,Type=String,Description="Mutation">')
    outputfile.header.add_line('##INFO=<ID=MutationType,Number=1,Type=String,Description="Mutation type">')
    outputfile.header.add_line('##FORMAT=<ID=DP,Number=1,Type=String,Description="Allele Depth">')
    outputfile.header.add_line('##FORMAT=<ID=AF,Number=1,Type=String,Description="Allele Frequency">')
    # outputfile.header.add_line('##INFO=<ID=,Number=1,Type=String,Description="Mutation type">')
    return outputfile


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
                homList[homo][pot] = round(tmp_dis[pot] / readNum, 6)

    for homo in range(1, maxRepeat + 1):
        if homo not in homList:
            homList[homo] = {homo: 1}
        x, y = [], []
        for i in homList[homo]:
            x.append(i)
            y.append(homList[homo][i])
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

    vcffile = pysam.VariantFile(path_dis_parameter + "tmp_motif_" + motif + ".bcf", "rb")
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
    # return homList

    maxture = {}

    for first in range(1, maxRepeat):
        for second in range(1, maxRepeat):
            if first <= second:

                firstDis = homList[first]
                secondDis = homList[second]
                thismaxture = {}
                # print("+++++++++++++++++++++++++++++")
                # print(first, firstDis)
                # print(second, secondDis)
                for rp in set([i for i in firstDis] + [j for j in secondDis]):
                    if rp not in firstDis:
                        firstDis[rp] = 0
                    if rp not in secondDis:
                        secondDis[rp] = 0
                for rp in firstDis:
                    thismaxture[rp] = round((firstDis[rp] + secondDis[rp]) / 2, 6)

                maxture[first * 1000 + second] = removeZeroDict(thismaxture)
    with open(path_dis_parameter + "tmp_motif_" + motif + ".model", "w") as f:
        yaml.dump({"homo": homList, "maxture": maxture}, f)
    return {"homo": homList, "maxture": maxture}


def call(paras, motif, outputvcf):
    path_dis_parameter = paras["output"]
    modelfile = open(path_dis_parameter + "tmp_motif_" + motif + ".model", "r")
    model = yaml.load(modelfile)
    vcffile = pysam.VariantFile(path_dis_parameter + "tmp_motif_" + motif + ".bcf", "rb")
    for rec in vcffile.fetch():
        recordInfo = rec.info
        # motif = recordInfo["motif"]
        repeatTimes = int(recordInfo["repeatTimes"])
        disList = [list(map(int, i.split(":"))) for i in recordInfo["dis"].split("|")]
        thismaxRepeat = max([i[0] for i in disList])
        # maxRepeat = maxRepeat if maxRepeat >= thismaxRepeat else thismaxRepeat
        support_reads = int(recordInfo["support_reads"])
        disListnormal = {}
        for i in disList:
            disListnormal[i[0]] = i[1] / support_reads
        # print(disListnormal)
        minAllele = min(disListnormal.keys()) - 2
        maxAllele = max(disListnormal.keys()) + 2
        distance = {}
        for first in range(minAllele, maxAllele + 1):
            for second in range(minAllele, maxAllele + 1):
                if first <= second:
                    modelid = first * 1000 + second
                    if modelid in model["maxture"]:
                        distance[modelid] = getDisdistance(model["maxture"][modelid], disListnormal)
        distance_tuple = sorted(distance.items(), key=lambda kv: (kv[1], kv[0]))
        rec.info["FirstAlleles"] = distance_tuple[0][0]
        print(distance_tuple)


def errEval(parase):
    if not errEval_args_init(parase):
        print("[Error] Parameters error!")
        return -1
    args = get_value("paras")
    # motifList = classDisbyMotif(args)
    motifList = ["A", "T", "C", "G"]

    outputvcf = write_vcf_init_call(args["output"][:-1] + ".bcf", args["input"])
    model={}
    for motif in motifList:
        model[motif] = getOneMotifProsess(args, motif)
    with open(args["output"][:-1] +  ".model", "w") as f:
        yaml.dump(model, f)

        # call(args, motif, outputvcf)

    #
    #     print(motif)


if __name__ == "__main__":
    ""
