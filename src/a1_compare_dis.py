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

def getHomoNormalDis(motifDis_tmp, maxRepeat):
    homList = {}

    for homo in motifDis_tmp:
        if len(motifDis_tmp[homo])<5:
            continue
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
            continue
        x, y = [], []
        for i in homList[homo]:
            x.append(i)
            y.append(homList[homo][i])
        # print(x,y)
        # # print(homo,sumdis)
        # if homo>9:
        #     plt.plot(x,y)
        #     plt.scatter(x,y)
        #     plt.vlines(homo,0,1)
        #     plt.show()
        #
        #     plt.close()
        #
        #     print(homList[homo])
        # print(homo,homList[homo])
    return homList

def getDis(type="ILM", a_repeatLen=1, repeatRegion=[1, 30]):
    ilm = pysam.VariantFile(filePathList[type], "rb")

    repeatList=[]
    motifDis_tmp={}
    maxRepeat=0
    for rec in ilm.fetch("chr1", 0, endPos):
        motif = rec.info["motif"]
        motifLen = len(motif)
        repeatTimes = rec.info["repeatTimes"]
        # repeatLen = rec.info["repeatTimes"] * motifLen

        if motifLen == a_repeatLen:
            if repeatTimes < repeatRegion[0] or repeatTimes > repeatRegion[1]:
                continue
            disList = [list(map(int, i.split(":"))) for i in rec.info["dis"].split("|")]
            thismaxRepeat = max([i[0] for i in disList])
            maxRepeat = maxRepeat if maxRepeat >= thismaxRepeat else thismaxRepeat
            support_reads = rec.info["support_reads"]
            disListnormal = {}
            for i in disList:
                disListnormal[i[0]] = i[1] / support_reads
            repeatList.append(repeatTimes)
            if repeatTimes not in motifDis_tmp:
                motifDis_tmp[repeatTimes] = []
            thistmp = motifDis_tmp[repeatTimes]
            thistmp.append(disListnormal)
            motifDis_tmp[repeatTimes] = thistmp
    homList = getHomoNormalDis(motifDis_tmp, maxRepeat)
    return homList

def analysis(motifLen=1, repeatRegion=[5,19]):
    # lable_list = ["ILM", "CCS"]
    # repeatLen=motifLen
    # microRef = getRefNum(a_repeatLen=motifLen, repeatRegion=repeatRegion)
    microCCS = getDis(type="CCS", a_repeatLen=motifLen, repeatRegion=repeatRegion)
    microILM = getDis(type="ILM", a_repeatLen=motifLen, repeatRegion=repeatRegion)
    print("Finished data")
    lable = ["ILM", "CCS"]
    markers = ["*", ".", "+"]
    datalist = [microILM, microCCS]
    repeats=list(set([i for i in microILM]+[j for j in microCCS]))
    repeats.sort()
    row =5
    col=len(repeats)//row
    fig=plt.figure(figsize=(15, 8))
    # plt.axes(yscale="log")
    # fig, ax = plt.subplots(row, col, sharex='col', sharey='row')
    plotnum=0
    for rp in repeats:
        # thiscol=plotnum//col
        # thisrow=plotnum%col
        plotnum+=1
        plt.subplot(col,row,plotnum)
        ccs=microCCS[rp]
        x_ccs=list(ccs.keys())
        x_ccs.sort()
        y_ccs=[ccs[i] for i in ccs]
        ilm = microILM[rp]
        x_ilm = list(ilm.keys())
        x_ilm.sort()
        y_ilm = [ilm[i] for i in ilm]
        if plotnum<11:
            plt.xticks([])
        if (plotnum%row)!=1:
            plt.yticks([])
        else:
            plt.yticks([0,0.5,1])
        plt.plot(x_ilm,y_ilm,color="r")
        plt.plot(x_ccs,y_ccs,color ="blue")
        plt.text(repeatRegion[1]+2,0.85,str(rp))
        plt.xlim(repeatRegion[0]-5,repeatRegion[1]+5)
        plt.vlines(rp,0,1,color="gray")

        # plt.title(str(rp))

    plt.subplots_adjust(wspace=0.1, hspace=0.1) # 调整子图空白距离
    # fontsize = 19
    fig.legend(lable,loc ="right")
    # fig.set_xlabel('Month')
    # plt.xticks(fontsize=fontsize - 2)
    # plt.yticks(fontsize=fontsize - 2)


    # plt.set_xlabel()
    plt.suptitle("MotifLen:" + str(motifLen) + "   " +
              str(repeatRegion))

    #
    # plt.show()
    plt.savefig(path_project +
                "data/GIAB/dis/fig/detection_dis_m" +
                str(motifLen) + "repeat" + str(repeatRegion[0]) + "-" + str(repeatRegion[1]),
                dpi=400)

    plt.close()

analysis(motifLen=1, repeatRegion=[5,19])
analysis(motifLen=1, repeatRegion=[20,34])
analysis(motifLen=2, repeatRegion=[5,19])


#
# def getRefNum(a_repeatLen=1, repeatRegion=[1, 30]):
#     path_MS = "/mnt/project/REF/microsatellite/GRCh38_l5_m100_20200528.list"
#     Homo = {}
#     for line in open(path_MS):
#         lineinfo = line[:-1].split("\t")
#         chrom = lineinfo[0]
#         if chrom != "chr1":
#             if chrom == "chr2":
#                 break
#             continue
#         pos = int(lineinfo[1])
#         if pos > endPos:
#             break
#         motifLen = int(lineinfo[2])
#         repetTImes = int(lineinfo[4])
#         repeatLen = int(lineinfo[4]) * motifLen
#         if motifLen == a_repeatLen:
#             if repeatLen < repeatRegion[0] or repeatLen > repeatRegion[1]:
#                 continue
#             if repetTImes not in Homo:
#                 Homo[repetTImes] = 0
#             Homo[repetTImes] += 1
#     return Homo

