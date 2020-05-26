# =============================================================================
# Project : MShunter0.0.1
# Py Name: ScanMicrosatellites
# Author :
# Date : 20-05-25
# Email : pengjia@stu.xjtu.edu.cn
# Description : ''
# =============================================================================

import argparse
import multiprocessing
import os

import pandas as pd
import pysam

# pool = multiprocessing.Pool(processes=cores)
chromList = [str(i) for i in range(1, 23)] + ["chr" + str(i) for i in range(1, 23)] + ["chrM", "chrX", "chrY", "MT",
                                                                                       "X", "Y"]




class MSDeail:
    repeatDict = {}
    depth=0
    depthCall=0

    # processStat=False

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


global ARGS


# , dfMicroSatellites, caseInfo


def argumentProcress():
    """
    argument procress
    """
    global ARGS

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
        os.mkdir(ARGS["workspace"] + "detailInfo/")
        print("[INFO] The workspace is : " + ARGS["workspace"])
    return True


def loadMicroSatellite(ms):
    """
    :return:
    """
    dfMicroSatellites = pd.read_csv(ms, index_col=0)
    return dfMicroSatellites


def getRepeatTimes(alignment, motif, motifLen, prefix, suffix):
    """
    :param alignment:
    :param motif:
    :param motifLen:
    :param prefix:
    :param suffix:
    :return:
    """
    global ARGS
    if alignment.mapping_quality < ARGS["minimum_mapping_quality"]:
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
    depth=len(alignmentList)
    if depth < msDetail.min_support_reads:
        return False
    repeatTimesDict = {}
    for alignment in alignmentList:
        if alignment.is_unmapped: continue
        thisRepeatTimes = getRepeatTimes(alignment, msDetail.motif, msDetail.motifLen, msDetail.prefix, msDetail.suffix)
        if thisRepeatTimes < 0: continue
        if thisRepeatTimes not in repeatTimesDict: repeatTimesDict[thisRepeatTimes] = 0
        repeatTimesDict[thisRepeatTimes] += 1
    if sum(repeatTimesDict.values()) < msDetail.min_support_reads:
        return False
    else:
        msDetail.repeatDict = repeatTimesDict
        msDetail.depth=depth
        msDetail.support_reads=sum(list(repeatTimesDict.values()))
        # print(msDetail.repeatDict)
        return msDetail


def multiRun(thread, datalist):
    pool = multiprocessing.Pool(processes=thread)
    result_list = pool.map(processOneMs, datalist)
    pool.close()
    pool.join()

    # result_list = []
    # # print(datalist[1])
    # for msDetail in datalist:
    #     result_list.append(processOneMs(msDetail))

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
                 "depth="+str(msDetail.depth),
                 "support_reads="+str(msDetail.support_reads),
                 "dis=" + "|".join([str(key)+":"+str(value) for key, value in msDetail.repeatDict.items()])]
            )
            format="GT"
            value="./."

            outputfile.write("\t".join([chrom,pos,idnum,ref,alt,qual,filter,Info,format,value]) + "\n")
        # outputfile.write(site+"\n")


def getDis(args={}, upstreamLen=5, downstreamLen=5):
    bam = ARGS["input"]
    ms = ARGS["Microsatellite"]
    dis = ARGS["workspace"] + "test.dis"
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
                              # alignmentList=alignmentList,
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


def main():
    # global caseInfo, Result, arguments
    global ARGS
    if not argumentProcress():
        return 1
    getDis(args=ARGS)
    return 0


if __name__ == "__main__":
    if main() > 1:
        os.system("rm -r " + ARGS["workspace"])
