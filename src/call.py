from src.units import *
from src.global_dict import *
import pysam
# import yaml
import multiprocessing
from functools import partial
import copy

Lock = multiprocessing.Lock()
Model = {}


class MSCall:
    chr_id = ""
    pos = ""
    ref = ""
    info = {}
    dis = {}
    dis_norm = {}
    minAllele = 1
    maxAllele = 100
    modelStat = False
    model = {}
    distance = 0
    distance_dict={}
    firstAllels = ""
    secondAllels = ""
    qual = 0
    filter = ""
    precision = "NoReadSpan"  # LowCov, Fuzzy, High
    GT = (0, 0)
    alleles = "0:0"
    alt = (".",)
    firsttwoDistance = 0

    def __init__(self, chr_id, pos, info, sample):
        self.info = info
        self.chr_id = chr_id
        self.pos = pos
        self.sample = sample
        self.disStat = False if info["disStat"] == "False" else True
        self.ref = str(info["repeatTimes"]) + "[" + info["motif"] + "]"

    def getDisdistance2(self, dict1, dict2):
        dictkey = list(dict1.keys()) + list(dict2.keys())
        for key in dictkey:
            if key not in dict1:
                dict1[key] = 0
            if key not in dict2:
                dict2[key] = 0
        sum = 0
        for key in dictkey:
            err = dict1[key] - dict2[key]
            sum += (err * err)
        return round(np.sqrt(sum), 6)

    def getDisdistance(self, dict1, dict2):
        dictkey = list(dict1.keys()) + list(dict2.keys())
        list1 = []
        list2 = []
        for key in dictkey:
            if key not in dict1:
                list1.append(0)
            else:
                list1.append(dict1[key])
            if key not in dict2:
                list2.append(0)
            else:
                list2.append(dict2[key])
        # print(np.linalg.norm(np.array(list1)-np.array(list2)))
        # print(np.sqrt(np.sum(np.square(np.array(list1)-np.array(list2)))))
        return np.linalg.norm(np.array(list1) - np.array(list2))

    def modelpre(self):
        # global Model
        if not self.disStat:
            return
        else:
            model = get_value("model")
            # lock.release()
            disnormal = {}
            self.dis = [list(map(int, i.split(":"))) for i in self.info["dis"].split("|")]
            for i in self.dis:
                disnormal[i[0]] = i[1] / self.info["support_reads"]
            self.dis_norm = disnormal
            motif = self.info["motif"]
            maxRepeat = model[motif]["maxRepeat"]
            self.minAllele = max([min(disnormal.keys()) - 1, 1])
            self.maxAllele = min([max(disnormal.keys()) + 1, maxRepeat])

            model_id_list = []
            for first in range(self.minAllele, self.maxAllele + 1):
                for second in range(self.minAllele, self.maxAllele + 1):
                    if first <= second:
                        model_id_list.append(first * 1000 + second)
            # self.model_id_list = model_id_list
            thismodel = {}
            if motif not in model:
                return
            else:
                self.modelStat = True
                # print(model_id_list)
                for model_id in model_id_list:
                    thismodel[model_id] = model[motif]["maxture"][model_id]
            self.model = thismodel

    def mscall(self):
        qual = 0
        if not self.disStat:
            self.qual = qual
            self.precision = "NoReadSpan"
            self.filter = "NoReadSpan"
            return
        if len(self.model) < 2:
            self.qual = qual
            self.precision = "NoAvailableModel"
            self.filter = "NoAvailableModel"
            return
        distance_dict = {}
        for model_id in self.model:
            distance_dict[model_id] = getDisdistance(self.dis_norm, self.model[model_id])
        self.model = {}
        distance_tuple = sorted(distance_dict.items(), key=lambda kv: (kv[1], kv[0]))
        first_1 = distance_tuple[0][0] // 1000
        first_2 = distance_tuple[0][0] % 1000
        firsttwoDistanceRatio = (distance_tuple[1][1] - distance_tuple[0][1]) / (distance_tuple[0][1] + 0.000001)
        self.firsttwoDistance = distance_tuple[1][1] - distance_tuple[0][1]
        self.distance = distance_tuple[0][1]
        qual = firsttwoDistanceRatio * self.info["support_reads"]
        if self.info["lowSupport"] == "True":
            self.precision = "LowCov"
            self.filter = "LowCov"
        elif firsttwoDistanceRatio < 0.3:
            self.precision = "Fuzzy"
            self.filter = "Fuzzy"
        else:
            self.precision = "High"
            self.filter = "PASS"
        self.firstAllels = first_1
        self.secondAllels = first_2
        self.qual = qual

        if first_1 == first_2:
            if first_1 == self.info["repeatTimes"]:
                self.GT = (0, 0)
                # self.alt = (".")
            else:
                self.GT = (1, 1)
                self.alt = (str(first_1) + "[" + self.info["motif"] + "]",)

                # self.allels=(first_1,first_2)
        else:
            if self.info["repeatTimes"] in [first_1, first_2]:
                self.GT = (0, 1)
                if first_1 == self.info["repeatTimes"]:
                    self.alt = (str(first_2) + "[" + self.info["motif"] + "]",)
                else:
                    self.alt = (str(first_1) + "[" + self.info["motif"] + "]",)
            else:
                self.GT = (1, 2)
                if first_1 < first_2:
                    self.alt = (str(first_1) + "[" + self.info["motif"] + "]",
                                str(first_2) + "[" + self.info["motif"] + "]")
                else:
                    self.alt = (str(first_2) + "[" + self.info["motif"] + "]",
                                str(first_1) + "[" + self.info["motif"] + "]")

        if first_1 <= first_2:
            self.alleles = ":".join(list(map(str, [first_1, first_2])))
        else:
            self.alleles = ":".join(list(map(str, [first_2, first_1])))


def call_one_ms(msCall):
    msCall.modelpre()
    msCall.mscall()
    return msCall


def multicallMS(mscall_list, outputfile, thread=4):
    pool = multiprocessing.Pool(thread)
    result_list = pool.map(call_one_ms, mscall_list)
    pool.close()
    pool.join()
    write_call2vcf(result_list, outputfile)

    return result_list


def write_call2vcf_init():
    paras = get_value("paras")
    inputpath = paras["output_dis"]
    outputpath = paras["output_call"]
    inputfile = pysam.VariantFile(inputpath, "rb")
    outputfile = pysam.VariantFile(outputpath, "wb", header=inputfile.header)
    outputfile.header.add_line(
        '##INFO=<ID=FirstAlleles,Number=1,Type=String,Description="The first allele type of this point">')
    outputfile.header.add_line(
        '##INFO=<ID=SecondAlleles,Number=1,Type=String,Description="The second allele type of this point">')
    outputfile.header.add_line('##INFO=<ID=Quality,Number=1,Type=Float,Description="Genotype quality">')
    outputfile.header.add_line(
        '##INFO=<ID=Distance,Number=1,Type=Float,Description="Distance between two distribution.">')
    outputfile.header.add_line(
        '##INFO=<ID=FirstTwoDistance,Number=1,Type=Float,Description="Distance between two distribution.">')
    outputfile.header.add_line('##INFO=<ID=Precision,Number=1,Type=String,Description="Genotype quality level">')
    outputfile.header.add_line('##INFO=<ID=Alleles,Number=1,Type=String,Description="Alleles">')
    outputfile.header.add_line('##INFO=<ID=Mutation,Number=1,Type=String,Description="Mutation">')
    outputfile.header.add_line('##INFO=<ID=MutationType,Number=1,Type=String,Description="Mutation type">')
    outputfile.header.add_line('##FORMAT=<ID=DP,Number=1,Type=String,Description="Allele Depth">')
    outputfile.header.add_line('##FORMAT=<ID=AF,Number=1,Type=String,Description="Allele Frequency">')
    outputfile.header.add_line('##INFO=<ID=Type,Number=1,Type=String,Description="Mutation type">')
    return outputfile


def write_call2vcf(mscall_list, outputfile):
    for mscall in mscall_list:
        vcfrec = outputfile.new_record()
        # print("infoKey",vcfrec.info.keys())
        vcfrec.contig = mscall.chr_id
        # vcfrec.stop = pos + msDetail.repeatTimes * len(msDetail.motif)
        vcfrec.pos = mscall.pos
        vcfrec.ref = mscall.ref
        vcfrec.alts = mscall.alt
        vcfrec.info["chrom"] = mscall.chr_id
        vcfrec.info["pos"] = mscall.info["pos"]
        vcfrec.info["Start"] = mscall.info["Start"]
        vcfrec.info["End"] = mscall.info["End"]
        vcfrec.info["motif"] = mscall.info["motif"]
        vcfrec.info["repeatTimes"] = mscall.info["repeatTimes"]
        vcfrec.info["prefix"] = mscall.info["prefix"]
        vcfrec.info["depth"] = mscall.info["depth"]
        vcfrec.info["support_reads"] = mscall.info["support_reads"]
        vcfrec.info["dis"] = mscall.info["dis"]
        vcfrec.info["proD"] = mscall.info["proD"]
        vcfrec.info["proI"] = mscall.info["proI"]
        vcfrec.info["lowSupport"] = str(mscall.info["lowSupport"])
        vcfrec.info["disStat"] = str(mscall.info["disStat"])
        vcfrec.info["Precision"] = mscall.precision
        vcfrec.info["Quality"] = round(mscall.qual, 6)
        vcfrec.info["FirstTwoDistance"] = mscall.firsttwoDistance
        vcfrec.info["Distance"] = mscall.distance
        vcfrec.qual = round(mscall.qual, 6)
        # vcfrec.filter = mscall.filter
        vcfrec.info["Alleles"] = mscall.alleles
        vcfrec.samples[get_value("case")]["GT"] = mscall.GT
        vcfrec.samples[get_value("case")].phased = False
        outputfile.write(vcfrec)

    return outputfile


def write_call2vcf_close(outputfile):
    outputfile.close()


def call():
    paras = get_value("paras")
    path_dis_parameter = paras["output"]
    thread = paras["threads"]

    batch = paras["batch"]
    disvcfpath = paras["output_dis"]
    # modelfile = open(path_dis_parameter + "tmp_motif_" + motif + ".model", "r")
    # model = yaml.load(model)
    vcffile = pysam.VariantFile(disvcfpath, "rb")
    # min_support_reads=paras["minimum_support_reads"]
    # print(model.keys())
    curentMSNum = 0
    tmpWindow = []
    ms_number = get_value("ms_number")
    outputfile = write_call2vcf_init()
    sampleID = get_value("case")
    for rec in vcffile.fetch():
        # rec_info = dict(rec.info)
        # rec_format = dict(rec.format)
        # rec_sample = dict(rec.samples)
        # print(rec_info, rec_format, rec_sample)

        msCall = MSCall(chr_id=rec.chrom, pos=rec.pos, info=dict(rec.info), sample=sampleID)
        # print(id(tmpWindow))
        tmpWindow.append(msCall)
        curentMSNum += 1
        if curentMSNum > 10000 and paras["debug"]:
            break
        if curentMSNum % (batch * thread) == 0:
            print("[Info] MS calling: Total", ms_number, "microsatelite, processing:", curentMSNum - batch * thread + 1,
                  "-", curentMSNum, "(" + str(round(curentMSNum / (0.1 + ms_number) * 100, 2)) + "%)")
            multicallMS(tmpWindow, outputfile=outputfile, thread=thread)
            tmpWindow = []
    multicallMS(tmpWindow, outputfile=outputfile, thread=thread)
    print("[Info] MS calling: Total", ms_number, "microsatelite, finish all!")
