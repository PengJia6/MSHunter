import pysam

vcf = pysam.VariantFile(
    "/media/pengjia/ZS6Y_Feng_1/NA12878/phasedVcf/NA12878.whatshap.strandseq-pacbio-phasing.2017-07-19.vcf.gz", "rb")
for rec in vcf.fetch():
    # print(rec)
    if rec.pos==837214:
        print(rec.id)
        print(rec.format.keys(), rec.format.values()[0])
        print(rec.samples.keys(),type(rec.samples.values()[0].values()[0]))
        print(rec.samples.values()[0].phased)


    # print(i.format)
