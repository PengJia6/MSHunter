import os
from src.bam2dis import *
from src.call import *
from src.errEval import *


class MSHAP:
    prefix_len = 10
    suffix_len = 10
    repeat_length_dis = {}
    depth = 0
    depth_hap1 = 0
    depth_hap2 = 0
    depthCall = 0
    dis_stat = False
    support_reads = 0
    support_reads_hap1 = 0
    support_reads_hap2 = 0
    more_than_one_alleles = False
    start_pre = 0
    end_suf = 0

    def __init__(self, chr_id, posStart, posEnd, queryStart, queryEnd, motif, motifLen, repeat_times,
                 bamfile, input_format, reference, prefix, suffix,
                 prefix_str, suffix_str, fapath):
        self.chrId = chr_id
        self.posStart = posStart
        self.motif = motif
        self.motifLen = motifLen
        self.repeatTimes = repeat_times
        self.posEnd = posEnd
        self.queryStart = queryStart
        self.queryEnd = queryEnd
        self.bamfile = bamfile
        self.input_format = input_format
        self.reference = reference
        self.prefix = prefix
        self.suffix = suffix
        self.prefix_str = prefix_str
        self.suffix_str = suffix_str
        self.fapath = fapath
        self.start_pre = posStart - self.prefix_len
        self.end_suf = posEnd + self.suffix_len

    def get_reads_alignment(self):
        if self.input_format == "bam":
            bamfile = pysam.AlignmentFile(self.bamfile, "rb")
        else:
            bamfile = pysam.AlignmentFile(self.bamfile, mode="rb", reference_filename=self.reference)
        alignmentList = [alignment for alignment in bamfile.fetch(self.chrId, self.queryStart, self.queryEnd)]
        bamfile.close()
        depth = len(alignmentList)
        if depth < 1:
            self.lowSupport = True
        reads_com = []
        # print(len(alignmentList))

        for alignment in alignmentList:
            if alignment.is_unmapped or alignment.is_duplicate or alignment.is_secondary:
                # print("llflflfl")
                continue  # remove unmaped reads and duplicated reads
            # if alignment.query_alignment_start > self.posStart and alignment.cigartuples[0][0] == 0:
            #     continue  # remove reads outside the boundary
            # if alignment.query_alignment_end < self.posEnd and alignment.cigartuples[0][0] == 0:
            #     continue  # remove reads outside the boundary
            if alignment.reference_start > self.posStart or alignment.reference_end < self.posEnd:
                # print("lflflfl")
                # print("--------------")
                # print(self.posStart,self.posEnd)
                # print(alignment.query_alignment_start,alignment.query_alignment_end)
                # print(alignment.get_reference_positions())
                continue
            reads_com.append(alignment)

        self.depth = len(reads_com)

        return reads_com

    def get_dis(self):

        repeat_length_dict = {}
        reads = self.get_reads_alignment()
        for alignment in reads:
            # add other condition for read selection
            repeat_length = self.get_repeat_length(alignment)
            if repeat_length not in repeat_length_dict:
                repeat_length_dict[repeat_length] = 0
            repeat_length_dict[repeat_length] += 1
        self.repeat_length_dis = repeat_length_dict
        if len(repeat_length_dict) > 1:
            self.more_than_one_alleles = True
            # print(repeat_length_dict)
            # print("ldflldl")
        if len(repeat_length_dict) == 1:
            self.dis_stat = True

    def get_alignment_detail(self):


        print()

    def pos_convert_ref2read(self, ref_block: list, read_block: list, pos: int, direction="start") -> tuple:
        """
        @param direction:  start of end
        @param ref_block:
        @param read_block:
        @param start:
        @param end:
        @return:
        """
        # print("=================================")

        if direction == "start":
            block_index = 0
            for sub_ref_block in ref_block:
                if sub_ref_block[0] <= pos <= sub_ref_block[1]:
                    # print(sub_ref_block,pos)
                    break
                block_index += 1

            if pos == ref_block[block_index][1]:  # M|D M|I D|M
                read_pos = read_block[block_index][1]
                if ref_block[block_index][2] == 2:
                    self.start_pre = min([self.start_pre, ref_block[block_index][0]])
                return pos, read_pos

            if ref_block[block_index][2] == 0:  # match
                read_pos = read_block[block_index][0] + (pos - ref_block[block_index][0])
                return pos, read_pos
            elif ref_block[block_index][2] == 2:  # deletion
                pos = pos
                read_pos = read_block[block_index][1]
                self.start_pre = min([read_block[block_index][0], self.start_pre])

                # pos = ref_block[block_index - 1][1] - 1
                # read_pos = read_block[block_index - 1][1] - 1
                return pos, read_pos
            else:
                print("gjjgjdgffkfkfkfpooooo")
                pos = ref_block[block_index - 1][1] - 1
                read_pos = read_block[block_index - 1][1] - 1
                return pos, read_pos

        if direction == "end":
            block_index = len(ref_block) - 1
            for sub_ref_block in ref_block[::-1]:
                if sub_ref_block[0] <= pos <= sub_ref_block[1]:
                    # print(sub_ref_block,pos)
                    break
                block_index = block_index - 1
            # if pos==ref_block[block_index][1]:
            #     print(block_index,len(ref_block))
            #     print( pos,ref_block[block_index-1],ref_block[block_index],)
            #     print("end")
            if pos == ref_block[block_index][0]:  # D|M I|M M|D
                read_pos = read_block[block_index][0]
                if ref_block[block_index][2] == 2:
                    self.end_suf = max([self.end_suf, ref_block[block_index][1]])
                return pos, read_pos

            if ref_block[block_index][2] == 0:
                read_pos = read_block[block_index][0] + (pos - ref_block[block_index][0])
                return pos, read_pos
            elif ref_block[block_index][2] == 2:
                pos = pos
                read_pos = read_block[block_index][0]
                self.end_suf = max([self.end_suf, read_block[block_index][1]])
                # pos = ref_block[block_index + 1][0] + 1
                # read_pos = read_block[block_index + 1][0] + 1
                return pos, read_pos
            else:
                print("lfllflflfllflflffkfjfhvjfnjgffjfjjjjjjjj")
                # ref_block[block_index][0] == pos:
                # print("llff")
                pos = pos + 1
                read_pos = read_block[block_index - 1][0] - 1
                return pos, read_pos

    def get_repeat_length(self, alignment):
        align_start = alignment.reference_start
        ref_block = []
        read_block = []
        read_pos = 0
        ref_pos = align_start
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
                return -1

        if self.posStart >= ref_block[-1][1] or self.posStart <= ref_block[0][0]:
            print("pos2")
            return -1
        if self.posEnd >= ref_block[-1][1] or self.posEnd <= ref_block[0][0]:
            print("pos3")
            return -1

        ref_start, read_start = self.pos_convert_ref2read(ref_block, read_block, self.posStart, direction="start")
        ref_end, read_end = self.pos_convert_ref2read(ref_block, read_block, self.posEnd, direction="end")
        rpt = self.repeatTimes * self.motifLen + ((read_end - read_start) - (ref_end - ref_start))

        return rpt

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

    def get_snp_info(self):
        pre_content = 5
        suf_content = 5
        start = self.posStart
        end = self.posEnd
        fafile = pysam.FastaFile(self.reference)
        start_pos = start - pre_content
        end_pos = end + suf_content
        ref_seq = fafile.fetch(self.chrId, start_pos, end_pos)
        bamfile = pysam.AlignmentFile(self.bamfile)
        print('++++++++++++++++++')
        outfile = pysam.AlignmentFile("-", "w", template=bamfile, index_filename=self.bamfile + ".bai")
        for pot in bamfile.fetch(self.chrId, start_pos, end_pos):
            outfile.write(pot)
        for pot in outfile.pileup(self.chrId, start_pos, end_pos, truncate=True, index_filename=self.bamfile + ".bai"):
            print(pot)

        return

    def getRepeatTimes2(self, alignment):

        """
        :param alignment:
        :param motif:
        :param motifLen:
        :param prefix:
        :param suffix:
        :return:
        """

        # self.getRepeatTimes2(alignment)
        if alignment.mapping_quality < self.min_mapping_qual:
            return -1
        readString = alignment.query
        prefixState = readString.find(self.prefix)
        if prefixState < 0: return -1
        suffixState = readString.rfind(self.suffix)
        if suffixState < 0: return -3
        if prefixState + 5 >= suffixState: return -2
        while prefixState >= 0:
            count = 0
            start = prefixState + 5
            while start == readString.find(self.motif, start):
                count += 1
                start = readString.find(self.motif, start) + self.motifLen
            if (self.motifLen == 1 and count >= 1) or (self.motifLen > 1 and count >= 1):
                if start == readString.find(self.suffix, start):
                    return count
            prefixState = readString.find(self.prefix, prefixState + 1)
        return -4


def benchmark_init(args):
    """
    argument procress
    """
    paras = {}
    paras["input"] = args.input[0]
    paras["output"] = args.output[0]
    paras["microsatellite"] = args.microsatellite[0]
    paras["reference"] = args.reference[0]
    paras["separator"] = args.separator[0]
    paras["prefix_len"] = args.prefix_len[0]
    paras["suffix_len"] = args.suffix_len[0]
    paras["debug"] = args.debug[0]
    paras["only_homopolymer"] = args.only_homopolymers[0]
    paras["minimum_support_reads"] = args.minimum_support_reads[0]
    paras["threads"] = args.threads[0]
    paras["batch"] = args.batch[0]
    paras["ranges_of_repeat_times"] = {}

    for i in args.minimum_repeat_times[0].split(";"):
        unitRange, repeatRange = i.split(":")
        if "-" in unitRange:
            unitStart, unitEnd = tuple(map(int, unitRange.split("-")))
        else:
            unitStart = int(unitRange)
            unitEnd = unitStart
        repeatStart = int(repeatRange)
        # print(unitStart,unitEnd,"  ",repeatStart, repeatEnd)
        for ur in range(unitStart, unitEnd + 1):
            if ur not in paras["ranges_of_repeat_times"]:
                paras["ranges_of_repeat_times"][ur] = {}
            paras["ranges_of_repeat_times"][ur]["min"] = repeatStart
        for i in args.maximum_repeat_times[0].split(";"):
            # print(i)
            unitRange, repeatRange = i.split(":")
            if "-" in unitRange:
                unitStart, unitEnd = tuple(map(int, unitRange.split("-")))
            else:
                unitStart = int(unitRange)
                unitEnd = unitStart
            repeatStart = int(repeatRange)
            # print(unitStart,unitEnd,"  ",repeatStart, repeatEnd)
            for ur in range(unitStart, unitEnd + 1):
                if ur not in paras["ranges_of_repeat_times"]:
                    paras["ranges_of_repeat_times"][ur] = {}
                paras["ranges_of_repeat_times"][ur]["max"] = repeatStart
    error_stat = False
    if os.path.exists(paras["input"]):
        print("[INFO] The input file is : '" + paras["input"] + "'.")
    else:
        print('[ERROR] The input file '
              + paras["input"] + ' is not exist, please check again')
        error_stat = True

    if os.path.isfile(paras["microsatellite"]):
        print("[INFO] The microsatellites file  is : " + paras["microsatellite"])
    else:
        print('[ERROR] The microsatellites file '
              + paras["microsatellite"] + ' is not exist, please check again')
        error_stat = True
    if os.path.isfile(paras["reference"]):
        print("[INFO] The reference file is : '" + paras["reference"] + "'.")
    else:
        paras["reference"] = "" if paras["reference"] == "." else paras["reference"]
        print('[ERROR] The reference file ' + paras["reference"] + ' is not exist, please check again')
        error_stat = True
    if paras["input"][-4:] == "cram":
        paras["input_format"] = "cram"
        cramfile = pysam.AlignmentFile(paras["input"], mode="rb", reference_filename=paras["reference"])
        if not cramfile.has_index():
            print("[INFO] Build index for the input cram ...")
            pysam.index(paras["input"])
        cramfile.close()
    else:
        paras["input_format"] = "hap1"
        bamfile = pysam.AlignmentFile(paras["input"], mode="rb")
        if not bamfile.has_index():
            print("[INFO] Build index for the input bam ...")
            pysam.index(paras["input"])
        bamfile.close()
    if not os.path.exists(paras["output"]):
        print("[INFO] The output is : " + paras["output"] + ".")
    else:
        print(
            '[ERROR] The output ' + paras["output"] +
            ' is still exist! in case of overwrite files in this workspace, '
            'please check your script!')
        if not paras["debug"]:
            error_stat = True
    if error_stat: return False
    output_path = paras["output"]
    output_path = output_path if output_path[-1] == "/" else output_path + "/"
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    paras["output"] = output_path
    input_path = paras["input"]
    input_path = input_path[:-1] if input_path[-1] == "/" else input_path
    case = input_path.split("/")[-1].strip(".bam")
    case = case.strip(".cram")
    paras["output_dis"] = paras["output"] + case + ".dis.vcf.gz"
    paras["output_tmp"] = paras["output"] + case + "_tmp"
    if not os.path.exists(paras["output_tmp"]):
        os.makedirs(paras["output_tmp"])
    paras["output_model"] = paras["output"] + case + ".model"
    paras["output_call"] = paras["output"] + case + ".vcf.gz"
    set_value("case", case)
    set_value("paras", paras)
    return True


def bm_processOneMs(msDetail):
    msDetail.get_dis()
    # msDetail.calcuShiftProbability()
    # msDetail.get_snp_info()

    return msDetail


def bm_multiRun(thread, datalist):
    # pool = multiprocessing.Pool(processes=thread)
    # result_list = pool.map(bm_processOneMs, datalist)
    # pool.close()
    # pool.join()
    result_list = []
    for ms in datalist:
        result_list.append(bm_processOneMs(ms))

    # print("input",len(datalist))
    # print("output",len(result_list))
    return result_list


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
        thisMSDeail = MSHAP(chr_id=chr_id,
                            posStart=posStart,
                            posEnd=posEnd,
                            motif=info["motif"],
                            motifLen=motifLen,
                            repeat_times=int(info["repeatTimes"]),
                            queryStart=queryStart,
                            queryEnd=queryEnd,
                            bamfile=args["input"],
                            input_format=input_format,
                            reference=reference,
                            prefix=info['prefix'],
                            suffix=info['suffix'],
                            prefix_str=prefix_str,
                            suffix_str=suffix_str,
                            fapath=args["reference"],
                            )
        tmpWindow.append(thisMSDeail)
        curentMSNum += 1
        if curentMSNum < 16000 and args["debug"]:
            continue
            # break
        if curentMSNum % (batch * thread) == 0:
            print("[INFO] Bam2dis: Total", ms_number, "microsatelite, processing:", curentMSNum - batch * thread + 1,
                  "-", curentMSNum, "(" + str(round(curentMSNum / ms_number * 100, 2)) + "%)")
            result_list = bm_multiRun(thread=thread, datalist=tmpWindow)
            # write_vcf(outputfile, result_list)
            tmpWindow = []
    result_list = bm_multiRun(thread=thread, datalist=tmpWindow)
    fafile.close()
    # write_vcf(outputfile, result_list)
    # write_vcf_close(outputfile)
    print("[INFO] Bam2dis: Total", ms_number, "microsatelite, finish all!")


def benchmark(parase):
    print("1")
    if benchmark_init(parase):
        args = get_value("paras")
        getDis(args)
        print("hhh")
