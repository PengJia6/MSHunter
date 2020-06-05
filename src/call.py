from src.units import *
from src.global_dict import *
import pysam
import yaml


class MSCall:
    dis = {}
    dis_norm = {}
    minAllele = 1
    maxAllele = 100
    modelStat = False
    model = {}
    # model_id_list = []
    distance = {}
    firstAllels = ""
    secondAllels = ""
    qual = 0
    precision = "NoReadSpan"  # LowCov, Fuzzy, High

    def __init__(self, vcf_dis_rec):
        self.vcf_dis_rec = vcf_dis_rec
        self.disStat = False if vcf_dis_rec.info["disStat"] == "False" else True

    def getDisdistance(self, dict1, dict2):
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

    def modelpre(self, model):
        if not self.disStat:
            return
        else:
            disnormal = {}
            self.dis = [list(map(int, i.split(":"))) for i in self.vcf_dis_rec.info["dis"].split("|")]
            for i in self.dis:
                disnormal[i[0]] = i[1] / self.vcf_dis_rec.info["support_reads"]
            motif = self.vcf_dis_rec.info["motif"]
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
            # print("ldldl")

    def mscall(self):
        qual = 0
        if not self.disStat:
            self.qual = qual
            self.precision = "NoReadSpan"
            return
        if len(self.model) < 2:
            self.qual = qual
            self.precision = "NoAvailableModel"
            return
        distance = {}
        # print("gfjfjffj")
        for model_id in self.model:
            distance[model_id] = getDisdistance(self.dis_norm, self.model[model_id])
        self.model = {}
        distance_tuple = sorted(distance.items(), key=lambda kv: (kv[1], kv[0]))
        # qual=distance_tuple()

        if self.vcf_dis_rec.info["lowSupport"] == "True":
            self.qual = qual
        # rec.info["FirstAlleles"] = distance_tuple[0][0]
        print(distance_tuple)

    """
    if recordInfo["disStat"] == "False":
            continue
        motif = recordInfo["motif"]
        disList = [list(map(int, i.split(":"))) for i in recordInfo["dis"].split("|")]
        # thismaxRepeat = max([i[0] for i in disList])
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
    """


def call(model):
    paras = get_value("paras")
    path_dis_parameter = paras["output"]
    disvcfpath = paras["output_dis"]
    # modelfile = open(path_dis_parameter + "tmp_motif_" + motif + ".model", "r")
    # model = yaml.load(model)
    vcffile = pysam.VariantFile(disvcfpath, "rb")
    # min_support_reads=paras["minimum_support_reads"]
    print(model.keys())
    for rec in vcffile.fetch():
        msCall = MSCall(rec)
        # msCall.modelpre()
        msCall.modelpre(model)
        msCall.mscall()
