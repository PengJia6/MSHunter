# =============================================================================
# Project : MShunter0.0.1
# Py Name: bam2dis.py
# Author :
# Date : 20-05-25
# Email : pengjia@stu.xjtu.edu.cn
# Description : ''
# =============================================================================

import multiprocessing
import pysam
import pandas as pd
from src.global_dict import *


class MSDeail:
    repeatDict = {}
    depth = 0
    depthCall = 0
    disStat = False
    lowSupport = True
    p = 0
    q = 0

    def __init__(self, chr_id, posStart, posEnd, queryStart, queryEnd, motif, motifLen, repeat_times, prefix, suffix,
                 bamfile, min_support_reads, min_mapping_qual, input_format, reference):
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
        self.input_format = input_format
        self.reference = reference

    def get_dis(self):

        if self.input_format == "bam":
            bamfile = pysam.AlignmentFile(self.bamfile, "rb")
        else:
            bamfile = pysam.AlignmentFile(self.bamfile, mode="rb", reference_filename=self.reference)

        alignmentList = [alignment for alignment in bamfile.fetch(self.chrId, self.queryStart, self.queryEnd)]
        bamfile.close()
        depth = len(alignmentList)
        if depth < self.min_support_reads:
            self.lowSupport = True
        repeatTimesDict = {}
        for alignment in alignmentList:
            if alignment.is_unmapped: continue
            # add other condition for read selection
            thisRepeatTimes = self.getRepeatTimes(alignment, self.motif, self.motifLen, self.prefix, self.suffix,
                                                  min_mapping_qual=self.min_mapping_qual)
            if thisRepeatTimes < 0: continue
            if thisRepeatTimes not in repeatTimesDict: repeatTimesDict[thisRepeatTimes] = 0
            repeatTimesDict[thisRepeatTimes] += 1
        if sum(repeatTimesDict.values()) < self.min_support_reads:
            self.lowSupport = True
            # return False
        else:

            self.lowSupport = False
        self.repeatDict = repeatTimesDict
        self.depth = depth
        self.support_reads = sum(list(repeatTimesDict.values()))
        if self.support_reads >= 1:
            self.disStat = True

    def calcuShiftProbability(self):
        """
        :param disDict:
        :param refRepeatTimes:
        :return:
        """
        # disDict=
        if not self.disStat:
            self.p = None
            self.q = None
            return False
        refRepeatTimes = self.repeatTimes
        insShfit = 0
        delShfit = 0
        normal = 0
        for rpt in self.repeatDict:
            if rpt - refRepeatTimes > 0:
                insShfit = insShfit + (rpt - refRepeatTimes) * self.repeatDict[rpt]
                normal = normal + rpt * self.repeatDict[rpt] - (rpt - refRepeatTimes) * self.repeatDict[rpt]
            else:
                delShfit = delShfit + (refRepeatTimes - rpt) * self.repeatDict[rpt]
                normal = normal + rpt * self.repeatDict[rpt] - (refRepeatTimes - rpt) * self.repeatDict[rpt]
        # print()
        self.p = round(delShfit / (insShfit + delShfit + normal), 4)
        self.q = round(insShfit / (insShfit + delShfit + normal), 4)
        # print(self.p,self.q)
        return True

    def getRepeatTimes(self, alignment, motif, motifLen, prefix, suffix, min_mapping_qual=0):
        """
        :param alignment:
        :param motif:
        :param motifLen:
        :param prefix:
        :param suffix:
        :return:
        """
        if alignment.mapping_quality < min_mapping_qual:
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
                    return count
            prefixState = readString.find(prefix, prefixState + 1)
        return -4


def loadMicroSatellite(args):
    """
    :return:
    """
    print("[INFO] Loading microsatellite file...")
    ms = args["microsatellite"]
    separator = args["separator"]
    # ID,chr,pos,motif,motifLen,repeatTimes,prefix,suffix
    if separator == "comma":
        dfMicroSatellites = pd.read_csv(ms, index_col=0)
        # print(dfMicroSatellites.columns)
    elif separator == "tab":
        dfMicroSatellites = pd.read_table(ms, header=0)
        columns = dfMicroSatellites.columns
        if "chromosome" in columns:
            dfMicroSatellites.rename(columns={"chromosome":"chr"}, inplace = True)
            # dfMicroSatellites["chr"] = dfMicroSatellites["chromosome"]
            # del dfMicroSatellites["chromosome"]
        if "location" in columns:
            dfMicroSatellites.rename(columns={"location":"pos"}, inplace = True)
            # dfMicroSatellites["pos"] = dfMicroSatellites["location"]
            # del dfMicroSatellites["location"]
        if "repeat_unit_bases" in columns:
            dfMicroSatellites.rename(columns={"repeat_unit_bases":"motif"}, inplace = True)
            # dfMicroSatellites["motif"] = dfMicroSatellites["repeat_unit_bases"]
            # del dfMicroSatellites["repeat_unit_bases"]
        if "repeat_unit_length" in columns:
            dfMicroSatellites.rename(columns={"repeat_unit_length":"motifLen"}, inplace = True)
            # dfMicroSatellites["motifLen"] = dfMicroSatellites["repeat_unit_length"]
            # del dfMicroSatellites["repeat_unit_length"]
        if "repeat_times" in columns:
            dfMicroSatellites.rename(columns={"repeat_times":"repeatTimes"}, inplace = True)
            # dfMicroSatellites["repeatTimes"] = dfMicroSatellites["repeat_times"]
            # del dfMicroSatellites["repeat_times"]
        if "left_flank_bases" in columns:
            dfMicroSatellites.rename(columns={"left_flank_bases":"prefix"}, inplace = True)
            # dfMicroSatellites["prefix"] = dfMicroSatellites["left_flank_bases"]
            # del dfMicroSatellites["left_flank_bases"]
        if "right_flank_bases" in columns:
            dfMicroSatellites.rename(columns={"right_flank_bases":"suffix"}, inplace = True)
            # dfMicroSatellites["suffix"] = dfMicroSatellites["right_flank_bases"]
            # del dfMicroSatellites["right_flank_bases"]
        dfMicroSatellites.index = dfMicroSatellites["chr"] + "_" + dfMicroSatellites["pos"].astype(str)
    elif separator == "space":
        dfMicroSatellites = pd.read_table(ms, header=0, sep=" ")
    if args["only_homopolymer"]:
        dfMicroSatellites = dfMicroSatellites[dfMicroSatellites['motifLen'] == 1]
    if args["debug"]:
        # dfMicroSatellites = dfMicroSatellites[dfMicroSatellites["motifLen"] == 1]
        if len(dfMicroSatellites) > 3000000:
            dfMicroSatellites = dfMicroSatellites.sample(30000)
    chromList = get_value("chrom_list")
    dfMicroSatellites = dfMicroSatellites[dfMicroSatellites['chr'].isin(chromList)]
    repeatRange = args["ranges_of_repeat_times"]
    repeatUnitList = sorted(repeatRange.keys())

    newDf = pd.DataFrame()
    for ul in repeatUnitList:
        minr = repeatRange[ul]["min"]
        maxr = repeatRange[ul]["max"]
        newDf = pd.concat([newDf, dfMicroSatellites[(dfMicroSatellites["motifLen"] == ul) &
                                                    (dfMicroSatellites["repeatTimes"] >= minr) &
                                                    (dfMicroSatellites["repeatTimes"] <= maxr)
                                                    ]])

    print("[INFO] There are total", len(newDf), "microsatellites.")
    set_value("ms_number", len(newDf))
    set_value("motifList", set(newDf["motif"]))
    # print(set(newDf["motif"]))
    return newDf


def processOneMs(msDetail):
    msDetail.get_dis()
    msDetail.calcuShiftProbability()
    return msDetail


def multiRun(thread, datalist):
    pool = multiprocessing.Pool(processes=thread)
    result_list = pool.map(processOneMs, datalist)
    pool.close()
    pool.join()
    return result_list


def write_vcf_init(outputpath, sampleNameList):
    outputfile = pysam.VariantFile(outputpath, "wb")
    outputfile.header.add_line('##INFO=<ID=chrom,Number=1,Type=String,Description="Chromosome">')
    outputfile.header.add_line('##INFO=<ID=pos,Number=1,Type=Integer,Description="Position">')
    outputfile.header.add_line('##INFO=<ID=Start,Number=1,Type=Integer,Description="Position start">')
    outputfile.header.add_line('##INFO=<ID=End,Number=1,Type=Integer,Description="Position End">')
    outputfile.header.add_line('##INFO=<ID=motif,Number=1,Type=String,Description="Repeat unit">')
    outputfile.header.add_line('##INFO=<ID=repeatTimes,Number=1,Type=Integer,Description="Repeat imes">')
    outputfile.header.add_line('##INFO=<ID=prefix,Number=1,Type=String,Description="Prefix of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=suffix,Number=1,Type=String,Description="Suffix of microsatellite">')
    outputfile.header.add_line(
        '##INFO=<ID=depth,Number=1,Type=Integer,Description="Number of reads associated with the position">')
    outputfile.header.add_line(
        '##INFO=<ID=support_reads,Number=1,Type=Integer,Description="Reads covered the microsatellite">')
    outputfile.header.add_line('##INFO=<ID=dis,Number=1,Type=String,Description="Distribution of repeat times">')
    outputfile.header.add_line('##INFO=<ID=proD,Number=1,Type=Float,Description="Probability of deletion">')
    outputfile.header.add_line('##INFO=<ID=proI,Number=1,Type=Float,Description="Probability of insertion">')
    outputfile.header.add_line('##INFO=<ID=lowSupport,Number=1,Type=String,Description="Low support reads">')
    outputfile.header.add_line('##INFO=<ID=disStat,Number=1,Type=String,Description="Distribution Stat">')
    outputfile.header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    outputfile.header.add_line('##FORMAT=<ID=PS,Number=1,Type=String,Description="Genotype">')
    outputfile.header.add_sample(get_value("case"))
    # print("header",outputfile.header.samples)
    # print("header",outputfile.header.formats)
    return outputfile


def write_vcf_close(outputfile):
    outputfile.close()


def write_vcf(outputfile, dataList):
    # print(header.contigs)
    contigs = []
    for msDetail in dataList:
        chrom = msDetail.chrId
        pos = int(msDetail.posStart)
        ref = str(msDetail.repeatTimes) + "[" + msDetail.motif + "]"
        if chrom not in contigs:
            outputfile.header.add_line("##contig=<ID={chrom}>".format(chrom=chrom))
        vcfrec = outputfile.new_record()
        # print("infoKey",vcfrec.info.keys())
        vcfrec.contig = chrom
        # vcfrec.stop = pos + msDetail.repeatTimes * len(msDetail.motif)
        vcfrec.pos = pos
        vcfrec.ref = ref
        vcfrec.info["chrom"] = chrom
        vcfrec.info["pos"] = pos
        vcfrec.info["Start"] = pos
        vcfrec.info["End"] = pos + msDetail.repeatTimes * len(msDetail.motif)
        vcfrec.info["motif"] = msDetail.motif
        vcfrec.info["repeatTimes"] = msDetail.repeatTimes
        vcfrec.info["prefix"] = msDetail.prefix
        vcfrec.info["depth"] = msDetail.depth
        vcfrec.info["support_reads"] = msDetail.support_reads
        vcfrec.info["dis"] = "|".join([str(key) + ":" + str(value) for key, value in msDetail.repeatDict.items()])
        vcfrec.info["proD"] = msDetail.p
        vcfrec.info["proI"] = msDetail.q
        vcfrec.info["lowSupport"] = str(msDetail.lowSupport)
        vcfrec.info["disStat"] = str(msDetail.disStat)
        vcfrec.samples[get_value("case")]["GT"] = ()
        vcfrec.samples[get_value("case")].phased = True
        outputfile.write(vcfrec)


def getDis(args={}, upstreamLen=5, downstreamLen=5):
    chromList = get_value("chrom_list")
    dis = args["output_dis"]
    input_format = args["input_format"]
    reference = args["reference"]
    thread = args["threads"]
    batch = args["batch"]
    outputfile = write_vcf_init(dis, [get_value("case")])
    dfMicroSatellites = loadMicroSatellite(args)
    curentMSNum = 0
    tmpWindow = []
    ms_number = get_value("ms_number")
    for ms_id, info in dfMicroSatellites.iterrows():
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
                              min_mapping_qual=args["minimum_mapping_quality"],
                              input_format=input_format,
                              reference=reference
                              )
        tmpWindow.append(thisMSDeail)
        curentMSNum += 1
        if curentMSNum > 10000 and args["debug"]:
            break
        if curentMSNum % (batch * thread) == 0:
            print("[INFO] Bam2dis: Total", ms_number, "microsatelite, processing:", curentMSNum - batch * thread + 1,
                  "-", curentMSNum, "(" + str(round(curentMSNum / ms_number * 100, 2)) + "%)")
            result_list = multiRun(thread=thread, datalist=tmpWindow)
            write_vcf(outputfile, result_list)
            tmpWindow = []
    result_list = multiRun(thread=thread, datalist=tmpWindow)
    write_vcf(outputfile, result_list)
    write_vcf_close(outputfile)
    print("[INFO] Bam2dis: Total", ms_number, "microsatelite, finish all!")


def bam2dis():
    args = get_value("paras")
    getDis(args=args)


if __name__ == "__main__":
    ""
