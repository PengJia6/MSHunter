# =============================================================================
# Project : MShunter0.0.1
# Py Name: ScanMicrosatellites
# Author :
# Date : 20-05-25
# Email : pengjia@stu.xjtu.edu.cn
# Description : ''
# =============================================================================

import multiprocessing
import os

import pandas as pd
import pysam

from src.global_dict import *


class MSDeail:
    repeatDict = {}
    depth = 0
    depthCall = 0

    def __init__(self, chr_id, posStart, posEnd, queryStart, queryEnd, motif, motifLen, repeat_times, prefix, suffix,
                 bamfile, min_support_reads, min_mapping_qual):
        self.chrId = chr_id
        self.posStart = posStart
        self.motif = motif
        self.motifLen = motifLen
        self.repeatTimes = repeat_times
        self.posEnd = posEnd
        self.prefix = prefix
        self.suffix = suffix
        self.queryStart = queryStart
        self.queryEnd = queryEnd
        self.bamfile = bamfile
        self.min_support_reads = min_support_reads
        self.min_mapping_qual = min_mapping_qual


# , dfMicroSatellites, caseInfo


def bam2dis_args_init(args):
    """
    argument procress
    """
    paras = {}
    paras["input"] = args.input
    paras["output"] = args.output
    paras["Microsatellite"] = args.microsatellite[0]
    paras["threads"] = args.threads[0]
    paras["minimum_support_reads"] = args.minimum_support_reads[0]
    paras["minimum_mapping_quality"] = args.minimum_mapping_quality[0]
    paras["batch"] = args.batch[0]
    print(args.input)
    ErrorStat = False
    if len(paras["input"]) != len(set(paras["input"])) \
            or len(paras["output"]) != len(set(paras["output"])) \
            or len(paras["input"]) != len(paras["output"]):
        print(
            "[Error] If you specified option -i/--input and -o/--output multiple times,"
            " please make sure that the times of -i and -o should be equal and not repetitive\n"
            " correct: \n"
            "\t'-i case1.bam -i case2.bam -o case1.dis.vcf -o case2.dis.vcf'\n"
            "\t'-i case1.bam -o case1.dis.vcf -i case2.bam -o case2.dis.vcf'\n"
            " wrong: \n"
            "\t'-i case1.bam -i case2.bam -o case2.dis.vcf -o case1.dis.vcf'\n"
            "\t'-i case1.bam -o case1.dis.vcf -i case2.bam '\n"
        )
        ErrorStat = True
    bamnum = 0
    for onebam in paras["input"]:
        bamnum += 1
        if os.path.isfile(onebam):
            print("[INFO] The ", bamnum, " input is : " + onebam)
        else:
            print('[ERROR] The ', bamnum, 'input "' + onebam + '" is not exist, please check again')
            ErrorStat = True
            break
    if os.path.isfile(paras["Microsatellite"]):
        print("[INFO] The Microsatellites file  is : " + paras["Microsatellite"])
    else:
        print('[ERROR] The Microsatellites file "'
              + paras["Microsatellite"]
              + '" is not exist, please check again')
        ErrorStat = True
    bamnum = 0
    for onebam in paras["output"]:
        bamnum += 1
        if not os.path.isfile(onebam):
            print("[INFO] The ", bamnum, "output is : '" + onebam + "'.")
        else:
            print(
                '[ERROR] The ', bamnum,
                'output "' + onebam + '" is still exist! in case of overwrite files in this workspace, '
                                      'please check your script!')
            ErrorStat = True
            break
    if ErrorStat: return False

    set_value("paras", paras)
    return True


def loadMicroSatellite(ms):
    """
    :return:
    """
    dfMicroSatellites = pd.read_csv(ms, index_col=0)
    return dfMicroSatellites


def getRepeatTimes(alignment, motif, motifLen, prefix, suffix, min_support_reads, min_mapping_qual=0):
    """
    :param alignment:
    :param motif:
    :param motifLen:
    :param prefix:
    :param suffix:
    :return:
    """

    if alignment.mapping_quality < min_support_reads:
        return -1
    readString = alignment.query
    prefixState = readString.find(prefix)
    if prefixState < 0: return -1
    suffixState = readString.rfind(suffix)
    if suffixState < 0: return -3
    if prefixState + 5 >= suffixState: return -2
    while prefixState >= 0:
        count = 0
        start = prefixState + 5
        while start == readString.find(motif, start):
            count += 1
            start = readString.find(motif, start) + motifLen
        if (motifLen == 1 and count >= 1) or (motifLen > 1 and count >= 1):
            if start == readString.find(suffix, start):
                # print(count, "    ", prefix,motif, suffix,repeat)
                return count
        prefixState = readString.find(prefix, prefixState + 1)
    return -4


def processOneMs(msDetail):
    bamfile = pysam.AlignmentFile(msDetail.bamfile, "rb")
    alignmentList = [alignment for alignment in bamfile.fetch(msDetail.chrId, msDetail.queryStart, msDetail.queryEnd)]
    depth = len(alignmentList)
    if depth < msDetail.min_support_reads:
        return False
    repeatTimesDict = {}
    for alignment in alignmentList:
        if alignment.is_unmapped: continue
        thisRepeatTimes = getRepeatTimes(alignment, msDetail.motif, msDetail.motifLen, msDetail.prefix, msDetail.suffix,
                                         msDetail.min_support_reads)
        if thisRepeatTimes < 0: continue
        if thisRepeatTimes not in repeatTimesDict: repeatTimesDict[thisRepeatTimes] = 0
        repeatTimesDict[thisRepeatTimes] += 1
    if sum(repeatTimesDict.values()) < msDetail.min_support_reads:
        return False
    else:
        msDetail.repeatDict = repeatTimesDict
        msDetail.depth = depth
        msDetail.support_reads = sum(list(repeatTimesDict.values()))
        return msDetail


def multiRun(thread, datalist):
    pool = multiprocessing.Pool(processes=thread)
    result_list = pool.map(processOneMs, datalist)
    pool.close()
    pool.join()
    return result_list


def write_init(outputpath, sampleNameList):
    outputfile = open(outputpath, "w")
    clomum = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + sampleNameList)
    outputfile.write(clomum + "\n")
    return outputfile


def writeInfo(outputfile, dataList):
    for msDetail in dataList:
        if msDetail != False:
            chrom = msDetail.chrId
            pos = str(msDetail.posStart)
            idnum = "."
            ref = "[" + str(msDetail.repeatTimes) + "]" + msDetail.motif
            alt = "."
            qual = "."
            filter = "."
            Info = ";".join(
                ["chrom=" + chrom, "pos=" + pos, "motif=" + msDetail.motif, "repeatTimes=" + str(msDetail.repeatTimes),
                 "prefix=" + msDetail.prefix, "suffix=" + msDetail.suffix,
                 "depth=" + str(msDetail.depth),
                 "support_reads=" + str(msDetail.support_reads),
                 "dis=" + "|".join([str(key) + ":" + str(value) for key, value in msDetail.repeatDict.items()])]
            )
            format = "GT"
            value = "./."
            outputfile.write("\t".join([chrom, pos, idnum, ref, alt, qual, filter, Info, format, value]) + "\n")


def getDis(args={}, upstreamLen=5, downstreamLen=5):
    chromList = get_value("chrom_list")
    bam = args["input"]
    ms = args["Microsatellite"]
    dis = args["output"]
    thread = args["threads"]
    batch = args["batch"]
    outputfile = write_init(dis, [bam])
    # outputfile.close()
    dfMicroSatellites = loadMicroSatellite(ms)
    # bamfile = pysam.AlignmentFile(bam, "rb")
    curentMSNum = 0
    tmpWindow = []
    for id, info in dfMicroSatellites.iterrows():
        chr_id = info["chr"]
        if chr_id not in chromList:
            continue
        posStart = int(info["pos"])
        repeat_times = int(info["repeatTimes"])
        motif = info["motif"]
        motifLen = len(motif)
        posEnd = posStart + motifLen * repeat_times
        queryStart = posStart - upstreamLen
        queryEnd = posEnd + downstreamLen
        thisMSDeail = MSDeail(chr_id=chr_id,
                              posStart=posStart,
                              posEnd=posEnd,
                              motif=info["motif"],
                              motifLen=motifLen,
                              repeat_times=int(info["repeatTimes"]),
                              queryStart=queryStart,
                              queryEnd=queryEnd,
                              prefix=info['prefix'],
                              suffix=info['suffix'],
                              bamfile=args["input"],
                              min_support_reads=args["minimum_support_reads"],
                              min_mapping_qual=args["minimum_mapping_quality"]
                              )
        tmpWindow.append(thisMSDeail)
        curentMSNum += 1
        if curentMSNum > 1000:
            break
        if curentMSNum % (batch * thread) == 0:
            print("[Info] Processing:", curentMSNum - batch * thread + 1, "-", curentMSNum)
            result_list = multiRun(thread=thread, datalist=tmpWindow)
            writeInfo(outputfile, result_list)
            tmpWindow = []
    result_list = multiRun(thread=thread, datalist=tmpWindow)
    writeInfo(outputfile, result_list)
    outputfile.close()


def bam2dis(parase):
    if not bam2dis_args_init(parase):
        print("[Error] Parameters error!")
        return -1
    args = get_value("paras")
    inputs, outputs = args["input"], args["output"]
    bamnum = 0
    for inputbampath, outputebampath in zip(inputs, outputs):
        bamnum += 1
        args["input"] = inputbampath
        args["output"] = outputebampath
        print("[Info] Precessing", bamnum, "file")
        getDis(args=args)


if __name__ == "__main__":
    ""
