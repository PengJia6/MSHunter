import matplotlib.pyplot as plt
import pysam

path_project = "/mnt/project/ProjectMSI/MSCalling/"
filePathList = {"CCS": path_project + "data/GIAB/dis/HG001_CCS_GRCh38.dis.bcf",  # 30fold
                "ILM": path_project + "data/GIAB/dis/HG001_NGS_GRCh38.dis.bcf"  # 30fold
                }

# NGS 248945271
# NGS 248946344
endPos = 248940000
lable_list = ["ILM", "CCS"]


def getDetectionNum(type="ILM", a_repeatLen=1, repeatRegion=[1, 30]):
    ilm = pysam.VariantFile(filePathList[type], "rb")
    Homo = {}
    for rec in ilm.fetch("chr1", 0, endPos):
        motif = rec.info["motif"]
        motifLen = len(motif)
        repetTImes = rec.info["repeatTimes"]
        # repeatLen = rec.info["repeatTimes"] * motifLen
        if motifLen == a_repeatLen:
            if repetTImes < repeatRegion[0] or repetTImes > repeatRegion[1]:
                continue
            if repetTImes not in Homo:
                Homo[repetTImes] = 0
            Homo[repetTImes] += 1
    return Homo


def getRefNum(a_repeatLen=1, repeatRegion=[1, 30]):
    path_MS = "/mnt/project/REF/microsatellite/GRCh38_l5_m100_20200528.list"
    Homo = {}
    for line in open(path_MS):
        lineinfo = line[:-1].split("\t")
        chrom = lineinfo[0]
        if chrom != "chr1":
            if chrom == "chr2":
                break
            continue
        pos = int(lineinfo[1])
        if pos > endPos:
            break
        motifLen = int(lineinfo[2])
        repetTImes = int(lineinfo[4])
        # repeatLen = int(lineinfo[4]) * motifLen
        if motifLen == a_repeatLen:
            if repetTImes < repeatRegion[0] or repetTImes > repeatRegion[1]:
                continue
            if repetTImes not in Homo:
                Homo[repetTImes] = 0
            Homo[repetTImes] += 1
    return Homo


def analysis(motifLen=1, repeatRegion=[1, 30]):
    # lable_list = ["ILM", "CCS"]
    # repeatLen=motifLen
    microRef = getRefNum(a_repeatLen=motifLen, repeatRegion=repeatRegion)
    microCCS = getDetectionNum(type="CCS", a_repeatLen=motifLen, repeatRegion=repeatRegion)
    microILM = getDetectionNum(type="ILM", a_repeatLen=motifLen, repeatRegion=repeatRegion)
    lable = ["Ref", "ILM", "CCS"]
    markers = ["*", ".", "+"]
    plt.figure(figsize=(10, 7))
    plt.axes(yscale="log")
    datalist = [microRef, microILM, microCCS]
    # print(len(microILM))
    typenum = 0
    for numdict in datalist:
        typenum += 1
        x = []
        y = []
        for repeat in sorted(numdict.keys()):
            # if repeat < 30:
            x.append(repeat)
            y.append(numdict[repeat])
        if typenum == 1:
            color = "gray"
        elif typenum == 2:
            color = "red"
        if typenum == 3:
            color = "blue"

        plt.plot(x, y,color=color)
        plt.scatter(x, y, marker=markers[typenum - 1],color = color)
    fontsize = 19
    plt.legend(lable, fontsize=fontsize - 2)
    plt.xticks(fontsize=fontsize-2)
    plt.yticks(fontsize=fontsize-2)
    plt.xlabel("Repeat times", fontsize=fontsize)
    plt.ylabel("No. of microsatellite", fontsize=fontsize)
    plt.title("MotifLen:" + str(motifLen) + "   " +
              str(repeatRegion), fontsize=fontsize )

    # plt.show()
    plt.savefig(path_project +
                "data/GIAB/dis/fig/detection_ability_m" +
                str(motifLen) + "repeat" + str(repeatRegion[0]) + "-" + str(repeatRegion[1]),
                dpi=400)

    plt.close()

for motifLen in range(1,6):
    small=20//motifLen
    repeatRegion=[small,1000]
    analysis(motifLen=motifLen, repeatRegion=repeatRegion)
    repeatRegion = [5, 1000]
    analysis(motifLen=motifLen, repeatRegion=repeatRegion)
# analysis(type="micro")
