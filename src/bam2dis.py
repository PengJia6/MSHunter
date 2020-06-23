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
import gzip
from src.global_dict import *


class MSDeail:
    repeat_dis = {"all": {}, "hap1": {}, "hap2": {}, "unphased": {}}
    depth = 0
    depth_hap1 = 0
    depth_hap2 = 0
    depthCall = 0
    disStat = False
    lowSupport = True
    support_reads = 0
    support_reads_hap1 = 0
    support_reads_hap2 = 0
    p = {"all": 0, "hap1": 0, "hap2": 0, "unphased": 0}
    q = {"all": 0, "hap1": 0, "hap2": 0, "unphased": 0}
    thishap = False  # this microsatellite is phased

    def __init__(self, chr_id, posStart, posEnd, queryStart, queryEnd, motif, motifLen, repeat_times,
                 bamfile, min_support_reads, min_mapping_qual, input_format, reference, prefix, suffix,
                 prefix_str, suffix_str, fapath, hap, tech):
        self.chrId = chr_id
        self.posStart = posStart
        self.motif = motif
        self.motifLen = motifLen
        self.repeatTimes = repeat_times
        self.posEnd = posEnd
        self.queryStart = queryStart
        self.queryEnd = queryEnd
        self.bamfile = bamfile
        self.min_support_reads = min_support_reads
        self.min_mapping_qual = min_mapping_qual
        self.input_format = input_format
        self.reference = reference
        self.prefix = prefix
        self.suffix = suffix
        self.prefix_str = prefix_str
        self.suffix_str = suffix_str
        self.fapath = fapath
        self.hap = hap
        self.tech = tech

    def get_reads_alignment(self):
        if self.input_format == "bam":
            bamfile = pysam.AlignmentFile(self.bamfile, "rb")
        else:
            bamfile = pysam.AlignmentFile(self.bamfile, mode="rb", reference_filename=self.reference)
        alignmentList = [alignment for alignment in bamfile.fetch(self.chrId, self.queryStart, self.queryEnd)]
        bamfile.close()
        depth = len(alignmentList)
        if depth < self.min_support_reads:
            self.lowSupport = True
        reads_com = []
        reads_hap1 = []
        reads_hap2 = []
        if self.hap:
            for alignment in alignmentList:
                # remove unmaped reads and duplicated reads
                if alignment.is_unmapped or alignment.is_duplicate:
                    continue
                # remove reads outside the boundary
                if alignment.query_alignment_start > self.posStart and alignment.cigartuples[0][0] == 0:
                    continue
                if alignment.query_alignment_end < self.posEnd and alignment.cigartuples[0][0] == 0:
                    continue

                if not alignment.has_tag("HP"):
                    reads_com.append(alignment)
                else:
                    if alignment.get_tag("HP") == 1:
                        reads_hap1.append(alignment)
                    else:
                        reads_hap2.append(alignment)
            if len(reads_hap1) > 3 and len(reads_hap2) > 3:
                self.thishap = True
            else:
                reads_com = reads_com + reads_hap1 + reads_hap2
            self.depth_hap1 = len(reads_hap1)
            self.depth_hap2 = len(reads_hap2)
            self.depth = len(reads_com) + self.depth_hap1 + self.depth_hap2
        else:
            for alignment in alignmentList:
                if alignment.is_unmapped or alignment.is_duplicate:
                    continue  # remove unmaped reads and duplicated reads
                if alignment.query_alignment_start > self.posStart and alignment.cigartuples[0][0] == 0:
                    continue  # remove reads outside the boundary
                if alignment.query_alignment_end < self.posEnd and alignment.cigartuples[0][0] == 0:
                    continue  # remove reads outside the boundary
                reads_com.append(alignment)

            self.depth = len(reads_com)

        return {"unphased": reads_com, "hap1": reads_hap1, "hap2": reads_hap2}

    def get_dis(self):
        reads = self.get_reads_alignment()
        repeat_times_dict = {"unphased": {}, "hap1": {}, "hap2": {}}

        if not self.thishap:
            repeatTimesDict = {}
            for alignment in reads["unphased"]:
                # add other condition for r ead selection
                thisRepeatTimes = self.getRepeatTimes(alignment, self.motif, self.motifLen, self.prefix, self.suffix,
                                                      min_mapping_qual=self.min_mapping_qual)
                # self.getRepeatTimes2(alignment)
                if thisRepeatTimes < 0: continue
                if thisRepeatTimes not in repeatTimesDict: repeatTimesDict[thisRepeatTimes] = 0
                repeatTimesDict[thisRepeatTimes] += 1
            repeat_times_dict["unphased"] = repeatTimesDict
        else:
            repeatTimesDict = {}
            for alignment in reads["hap1"]:
                # add other condition for read selection
                thisRepeatTimes = self.getRepeatTimes(alignment, self.motif, self.motifLen, self.prefix, self.suffix,
                                                      min_mapping_qual=self.min_mapping_qual)
                # self.getRepeatTimes2(alignment)
                if thisRepeatTimes < 0: continue
                if thisRepeatTimes not in repeatTimesDict: repeatTimesDict[thisRepeatTimes] = 0
                repeatTimesDict[thisRepeatTimes] += 1
            repeat_times_dict["hap1"] = repeatTimesDict
            self.repeatDict = repeatTimesDict
            repeatTimesDict = {}
            for alignment in reads["hap2"]:
                # add other condition for read selection
                thisRepeatTimes = self.getRepeatTimes(alignment, self.motif, self.motifLen, self.prefix, self.suffix,
                                                      min_mapping_qual=self.min_mapping_qual)
                # self.getRepeatTimes2(alignment)
                if thisRepeatTimes < 0: continue
                if thisRepeatTimes not in repeatTimesDict: repeatTimesDict[thisRepeatTimes] = 0
                repeatTimesDict[thisRepeatTimes] += 1
            repeat_times_dict["unphased"] = repeatTimesDict
            self.repeatDict = repeatTimesDict
            repeatTimesDict = {}
            for alignment in reads["hap2"]:
                # add other condition for read selection
                thisRepeatTimes = self.getRepeatTimes(alignment, self.motif, self.motifLen, self.prefix, self.suffix,
                                                      min_mapping_qual=self.min_mapping_qual)
                self.getRepeatTimes2(alignment)
                if thisRepeatTimes < 0: continue
                if thisRepeatTimes not in repeatTimesDict: repeatTimesDict[thisRepeatTimes] = 0
                repeatTimesDict[thisRepeatTimes] += 1
            repeat_times_dict["unphased"] = repeatTimesDict
            self.repeatDict = repeatTimesDict
        repeat_times_dict["all"] = self.dis_sum(
            [repeat_times_dict["unphased"],
             repeat_times_dict["hap1"],
             repeat_times_dict["hap2"]])
        self.repeat_dis = repeat_times_dict
        if sum(repeat_times_dict["all"].values()) < self.min_support_reads:
            self.lowSupport = True
        else:
            self.lowSupport = False
        self.support_reads = sum(list(repeat_times_dict["all"].values()))
        self.support_reads_hap1 = sum(list(repeat_times_dict["hap1"].values()))
        self.support_reads_hap2 = sum(list(repeat_times_dict["hap2"].values()))
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
            # self.p = None
            # self.q = None
            return False
        refRepeatTimes = self.repeatTimes
        insShfit = 0
        delShfit = 0
        normal = 0
        p_dict = {}
        q_dict = {}
        for reads_type in self.repeat_dis:
            one_repeat_dict = self.repeat_dis[reads_type]
            for rpt in one_repeat_dict:
                if rpt - refRepeatTimes > 0:
                    insShfit = insShfit + (rpt - refRepeatTimes) * one_repeat_dict[rpt]
                    normal = normal + rpt * one_repeat_dict[rpt] - (rpt - refRepeatTimes) * one_repeat_dict[rpt]
                else:
                    delShfit = delShfit + (refRepeatTimes - rpt) * one_repeat_dict[rpt]
                    normal = normal + rpt * one_repeat_dict[rpt] - (refRepeatTimes - rpt) * one_repeat_dict[rpt]
            # print()
            if insShfit + delShfit + normal > 0:
                p_dict[reads_type] = round(delShfit / (insShfit + delShfit + normal), 4)
                q_dict[reads_type] = round(insShfit / (insShfit + delShfit + normal), 4)
            else:
                p_dict[reads_type] = 0
                q_dict[reads_type] = 0
        self.p = p_dict
        self.q = q_dict
        # print(self.p,self.q)
        return True

    # ref                     -----=======-----
    # read                       ---------->              spaning
    # read                       <-----------

    # ref                     -----=======-----
    # read                ---->  <---|                    remove
    # read                            |--->  <----        remove

    # ref                     -----=======-----
    # read                 ---->   <---/                  probably remove
    # read                            \--->  <----        probably remove

    # ref                     -----=======-----
    # read                 ---->         <---             read pair info
    # read                      ---->       <----         read pair info
    # read                       ---->   <----            read pair info

    def getRepeatTimes2(self, alignment):

        align_start = alignment.reference_start
        align_end = alignment.reference_end
        query_start = alignment.query_alignment_start
        query_end = alignment.query_alignment_end
        query_length = alignment.query_length
        align_length = alignment.query_alignment_length
        # print(self.posStart)
        ref_block = []
        read_block = []
        read_pos = 0
        ref_pos = align_start
        # print(alignment.cigarstring,alignment.cigartuples,)
        for cigartupe in alignment.cigartuples:
            if cigartupe[0] in [0, 7, 8]:  # 0 : M : match or mishmatch ; 7: :=:match; 8:X:mismatch
                ref_block.append((ref_pos, ref_pos + cigartupe[1], 0))
                read_block.append((read_pos, read_pos + cigartupe[1], 0))
                read_pos += cigartupe[1]
                ref_pos += cigartupe[1]
            elif cigartupe[0] in [1, 4, 5]:  # 1:I:inserion ;4:S:soft clip 5:H:hardclip
                ref_block.append((ref_pos, ref_pos + 0, 1))
                read_block.append((read_pos, read_pos + cigartupe[1], 1))
                read_pos += cigartupe[1]
            elif cigartupe[0] in [2, ]:  # 2:D; 3:N: skip region of reference
                ref_block.append((ref_pos, ref_pos + cigartupe[1], 2))
                read_block.append((read_pos, read_pos, 2))
                ref_pos += cigartupe[1]
            else:
                print(alignment.cigarstring)
        if ref_pos != align_end:
            print(alignment)

        if read_pos != (query_length):
            print(read_pos, query_length, query_end)
        # print(alignment)
        # read_start, read_end = 0, alignment.query_length
        # for ref_b, read_b in zip(ref_block, read_block):
        #     if self.posStart >= ref_b[0] and self.posStart < ref_b[1]:
        #         read_start = read_b[0] + (self.posStart - ref_b[0]) + 1
        #     if self.posEnd >= ref_b[0] and self.posEnd < ref_b[1]:
        #         read_end = read_b[0] + (self.posEnd - ref_b[0]) - 1
        # print(read_start, read_end, cigartupe)
        # fafile = pysam.FastaFile(self.fapath)
        # # print(fafile.fetch(self.chrId,self.posStart,self.posStart+self.repeatTimes*self.motifLen))
        # print(alignment.query[read_start - 6:read_start])
        # print(alignment.query)
        #
        # print(self.prefix)
        # print()
        # print(alignment.query[read_end:read_end + 6])
        # # print(self.suffix)
        # print()
        # print()
        # fafile.close()
        # # elif cigartupe[0] in [5]:
        #     continue
        # else:
        #
        # if align_end-align_start!=align_length:
        #
        #     print(len(alignment.get_reference_positions()),len(alignment.get_reference_positions(full_length=True)))
        #     print(alignment.pos,alignment.reference_start, alignment.reference_end,align_end-align_start,align_length)
        #     print(query_start, query_end, query_length,query_end-query_start)
        #     print("")

        # if query_length!=query_end-query_start:
        #     print(alignment.cigarstring,alignment.cigartuples)
        #     print(alignment.pos,alignment.reference_start, alignment.reference_end,align_end-align_start,align_length)
        #     print(query_start, query_end, query_length,query_end-query_start)
        #     print("")

        # print(alignment.to_dict())

    def dis_sum(self, dict_list):
        keylist = []
        for item in dict_list:
            for key in item:
                if key not in keylist:
                    keylist.append(key)
        res_dict = {}
        for key in keylist:
            res_dict[key] = 0
        for item in dict_list:
            for key in item:
                res_dict[key] += item[key]
        return res_dict

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
            dfMicroSatellites.rename(columns={"chromosome": "chr"}, inplace=True)
            # dfMicroSatellites["chr"] = dfMicroSatellites["chromosome"]
            # del dfMicroSatellites["chromosome"]
        if "location" in columns:
            dfMicroSatellites.rename(columns={"location": "pos"}, inplace=True)
            # dfMicroSatellites["pos"] = dfMicroSatellites["location"]
            # del dfMicroSatellites["location"]
        if "repeat_unit_bases" in columns:
            dfMicroSatellites.rename(columns={"repeat_unit_bases": "motif"}, inplace=True)
            # dfMicroSatellites["motif"] = dfMicroSatellites["repeat_unit_bases"]
            # del dfMicroSatellites["repeat_unit_bases"]
        if "repeat_unit_length" in columns:
            dfMicroSatellites.rename(columns={"repeat_unit_length": "motifLen"}, inplace=True)
            # dfMicroSatellites["motifLen"] = dfMicroSatellites["repeat_unit_length"]
            # del dfMicroSatellites["repeat_unit_length"]
        if "repeat_times" in columns:
            dfMicroSatellites.rename(columns={"repeat_times": "repeatTimes"}, inplace=True)
            # dfMicroSatellites["repeatTimes"] = dfMicroSatellites["repeat_times"]
            # del dfMicroSatellites["repeat_times"]
        if "left_flank_bases" in columns:
            dfMicroSatellites.rename(columns={"left_flank_bases": "prefix"}, inplace=True)
            # dfMicroSatellites["prefix"] = dfMicroSatellites["left_flank_bases"]
            # del dfMicroSatellites["left_flank_bases"]
        if "right_flank_bases" in columns:
            dfMicroSatellites.rename(columns={"right_flank_bases": "suffix"}, inplace=True)
            # dfMicroSatellites["suffix"] = dfMicroSatellites["right_flank_bases"]
            # del dfMicroSatellites["right_flank_bases"]
        dfMicroSatellites.index = dfMicroSatellites["chr"] + "_" + dfMicroSatellites["pos"].astype(str)
    elif separator == "space":
        dfMicroSatellites = pd.read_table(ms, header=0, sep=" ")
    chromList = get_value("chrom_list")
    dfMicroSatellites = dfMicroSatellites[dfMicroSatellites['chr'].isin(chromList)]
    if args["only_homopolymer"]:
        dfMicroSatellites = dfMicroSatellites[dfMicroSatellites['motifLen'] == 1]

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
    if args["debug"]:
        locis_num = 80000
        newDf = newDf.iloc[100:locis_num + 100, :]
        # dfMicroSatellites = dfMicroSatellites[dfMicroSatellites["motifLen"] == 1]
        # if len(newDf) > locis_num:
        #     newDf = newDf.sample(locis_num)

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
    # print("input",len(datalist))
    # print("output",len(result_list))
    return result_list


def write_vcf_init(outputpath, sampleNameList):
    outputfile = pysam.VariantFile(outputpath, "w")
    bamfile = pysam.AlignmentFile(get_value("paras")["input"], "rb")
    contigs = bamfile.references
    contigsLen = bamfile.lengths
    chromList = get_value("chrom_list")
    contigs_len_dict = {}
    sortedContig = []
    for contig, length in zip(contigs, contigsLen):
        if contig in chromList:
            sortedContig.append(contig)
            contigs_len_dict[contig] = length
    for contig in sortedContig:
        outputfile.header.add_line(
            "##contig=<ID={chrom},length={length}>".format(chrom=contig, length=contigs_len_dict[contig]))
    set_value("contigsInfo", contigs_len_dict)
    outputfile.header.add_line('##INFO=<ID=chrom,Number=1,Type=String,Description="Chromosome">')
    outputfile.header.add_line('##INFO=<ID=pos,Number=1,Type=Integer,Description="Position">')
    outputfile.header.add_line('##INFO=<ID=Start,Number=1,Type=Integer,Description="Position start">')
    outputfile.header.add_line('##INFO=<ID=End,Number=1,Type=Integer,Description="Position End">')
    outputfile.header.add_line('##INFO=<ID=motif,Number=1,Type=String,Description="Repeat unit">')
    outputfile.header.add_line('##INFO=<ID=repeatTimes,Number=1,Type=Integer,Description="Repeat imes">')
    outputfile.header.add_line('##INFO=<ID=prefix,Number=1,Type=String,Description="Prefix of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=suffix,Number=1,Type=String,Description="Suffix of microsatellite">')
    outputfile.header.add_line('##INFO=<ID=depth,Number=1,Type=String,Description='
                               '"Number of reads associated with the position, all|hap1|hap2">')
    outputfile.header.add_line('##INFO=<ID=support_reads,Number=1,Type=String,Description='
                               '"Reads covered the microsatellite, all|hap1|hap2">')
    outputfile.header.add_line('##INFO=<ID=dis,Number=1,Type=String,Description='
                               '"Distribution of repeat times, all|hap1|hap2|unphased">')
    outputfile.header.add_line('##INFO=<ID=proD,Number=1,Type=String,Description='
                               '"Probability of deletion,all|hap1|hap2|unphased">')
    outputfile.header.add_line('##INFO=<ID=proI,Number=1,Type=String,Description='
                               '"Probability of insertion, all|hap1|hap2|unphased">')
    outputfile.header.add_line('##INFO=<ID=lowSupport,Number=1,Type=String,Description="Low support reads">')
    outputfile.header.add_line('##INFO=<ID=Phased,Number=1,Type=String,Description="This site is phased">')
    outputfile.header.add_line('##INFO=<ID=disStat,Number=1,Type=String,Description="Distribution Stat">')
    outputfile.header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    outputfile.header.add_line('##FORMAT=<ID=PS,Number=1,Type=String,Description="Genotype">')
    outputfile.header.add_sample(get_value("case"))
    # print("header",outputfile.header.samples)
    # print("header",outputfile.header.formats)
    return outputfile


def write_vcf_close(outputfile):
    outputfile.close()
    pysam.tabix_index(get_value("paras")["output_dis"], force=True, preset="vcf")


def write_vcf(outputfile, dataList):
    # print(header.contigs)
    # print("write", len(dataList))
    for msDetail in dataList:
        chrom = msDetail.chrId
        pos = int(msDetail.posStart)
        ref = str(msDetail.repeatTimes) + "[" + msDetail.motif + "]"
        vcfrec = outputfile.new_record()
        # print("infoKey",vcfrec.info.keys())
        vcfrec.contig = chrom
        # vcfrec.stop = pos + msDetail.repeatTimes * len(msDetail.motif)
        vcfrec.pos = pos
        vcfrec.ref = ref
        vcfrec.info["chrom"] = chrom
        vcfrec.info["pos"] = pos
        vcfrec.info["Start"] = pos
        vcfrec.stop = pos + msDetail.repeatTimes * len(msDetail.motif)
        vcfrec.info["End"] = pos + msDetail.repeatTimes * len(msDetail.motif)
        vcfrec.info["motif"] = msDetail.motif
        vcfrec.info["repeatTimes"] = msDetail.repeatTimes
        vcfrec.info["prefix"] = msDetail.prefix
        vcfrec.info["suffix"] = msDetail.suffix
        vcfrec.info["depth"] = "|".join(list(map(str, [msDetail.depth, msDetail.depth_hap1, msDetail.depth_hap2])))
        vcfrec.info["support_reads"] = "|".join(list(map(str, [msDetail.support_reads,
                                                               msDetail.depth_hap1, msDetail.depth_hap2])))
        vcfrec.info["dis"] = "|".join([
            ":".join([str(key) + "-" + str(value) for key, value in msDetail.repeat_dis["all"].items()]),
            ":".join([str(key) + "-" + str(value) for key, value in msDetail.repeat_dis["hap1"].items()]),
            ":".join([str(key) + "-" + str(value) for key, value in msDetail.repeat_dis["hap2"].items()]),
            ":".join([str(key) + "-" + str(value) for key, value in msDetail.repeat_dis["unphased"].items()])

        ])

        vcfrec.info["proD"] = "|".join(list(map(str, [
            msDetail.p["all"], msDetail.p["hap1"], msDetail.p["hap2"], msDetail.p["unphased"]
        ])))
        vcfrec.info["proI"] = "|".join(list(map(str, [
            msDetail.q["all"], msDetail.q["hap1"], msDetail.q["hap2"], msDetail.q["unphased"]
        ])))
        vcfrec.info["lowSupport"] = str(msDetail.lowSupport)
        vcfrec.info["disStat"] = str(msDetail.disStat)
        vcfrec.samples[get_value("case")]["GT"] = ()
        vcfrec.samples[get_value("case")].phased = True
        outputfile.write(vcfrec)


def getDis(args={}, upstreamLen=5, downstreamLen=5):
    dis = args["output_dis"]
    input_format = args["input_format"]
    reference = args["reference"]
    thread = args["threads"]
    batch = args["batch"]
    outputfile = write_vcf_init(dis, [get_value("case")])
    contigs_info = get_value("contigsInfo")
    dfMicroSatellites = loadMicroSatellite(args)
    curentMSNum = 0
    tmpWindow = []
    ms_number = get_value("ms_number")
    tech = args["tech"]
    hap = args["hap"]
    fafile = pysam.FastaFile(args["reference"])
    prefix_len = args["prefix_len"]
    suffix_len = args["suffix_len"]
    # print(len(dfMicroSatellites))
    for ms_id, info in dfMicroSatellites.iterrows():
        chr_id = info["chr"]
        if chr_id not in contigs_info:
            continue
        posStart = int(info["pos"])
        repeat_times = int(info["repeatTimes"])
        motif = info["motif"]
        motifLen = len(motif)
        posEnd = posStart + motifLen * repeat_times
        queryStart = posStart - upstreamLen
        queryEnd = posEnd + downstreamLen
        prefix_str = fafile.fetch(chr_id, posStart - prefix_len, posStart)
        suffix_str = fafile.fetch(chr_id, posEnd, posEnd + suffix_len)
        # print(prefix_str,info["prefix"],motif)
        # print(suffix_str,info["suffix"],motif)
        thisMSDeail = MSDeail(chr_id=chr_id,
                              posStart=posStart,
                              posEnd=posEnd,
                              motif=info["motif"],
                              motifLen=motifLen,
                              repeat_times=int(info["repeatTimes"]),
                              queryStart=queryStart,
                              queryEnd=queryEnd,
                              bamfile=args["input"],
                              min_support_reads=args["minimum_support_reads"],
                              min_mapping_qual=args["minimum_mapping_quality"],
                              input_format=input_format,
                              reference=reference,
                              prefix=info['prefix'],
                              suffix=info['suffix'],
                              prefix_str=prefix_str,
                              suffix_str=suffix_str,
                              fapath=args["reference"],
                              hap=hap,
                              tech=tech,
                              )
        tmpWindow.append(thisMSDeail)
        curentMSNum += 1
        # if curentMSNum > 10000 and args["debug"]:
        #     break
        if curentMSNum % (batch * thread) == 0:
            print("[INFO] Bam2dis: Total", ms_number, "microsatelite, processing:", curentMSNum - batch * thread + 1,
                  "-", curentMSNum, "(" + str(round(curentMSNum / ms_number * 100, 2)) + "%)")
            result_list = multiRun(thread=thread, datalist=tmpWindow)
            write_vcf(outputfile, result_list)
            tmpWindow = []
    result_list = multiRun(thread=thread, datalist=tmpWindow)
    fafile.close()
    write_vcf(outputfile, result_list)
    write_vcf_close(outputfile)
    print("[INFO] Bam2dis: Total", ms_number, "microsatelite, finish all!")


def bam2dis():
    args = get_value("paras")
    getDis(args=args)


if __name__ == "__main__":
    ""
