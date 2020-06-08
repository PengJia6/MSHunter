# =============================================================================
# Project : mstools-0.0.1
# Py Name: global_dict
# Author : Peng Jia
# Date : 20-05-26
# Email : pengjia@stu.xjtu.edu.cn
# Description : 'For global variable'
# =============================================================================


def global_init():
    global _global_dict
    _global_dict = {}
    _global_dict["ms_number"] = 0
    _global_dict["tools_version"] = "1.0"
    _global_dict["tools_name"] = "mstools"
    _global_dict["chrom_list"] = [str(i) for i in range(1, 23)] + \
                                 ["chr" + str(i) for i in range(1, 23)] + \
                                 ["X", "Y", "chrX", "chrY", "chrM", "MT"]
    _global_dict["default"] = {
        "bam2dis": {
            "threads": 4,
            "minimum_mapping_quality": 1,
            "minimum_support_reads": 2,
            "batch": 2000,
            "debug": False,
            "separator": "tab",  # comma,tab,space
            "only_homopolymers": False,
            "allow_mismatch": True,
            "minimum_repeat_times": "1:8;2-5:5",
            "maximum_repeat_times": "1-5:100",

        },
        "errEval": {
            "threads": 4,
            "batch": 2000,

        },
        "call": {
        },
        "genotype": {
            "reference": ".",
            "threads": 4,
            "minimum_mapping_quality": 1,
            "minimum_support_reads": 2,
            "batch": 2000,
            "debug": False,
            "separator": "tab",  # comma,tab,space
            "only_homopolymers": False,
            "allow_mismatch": True,
            "minimum_repeat_times": "1:8;2-5:5",
            "maximum_repeat_times": "1-5:100",

        },

    }

def set_value(name, value):
    _global_dict[name] = value

def get_value(name, defValue=None):
    try:
        return _global_dict[name]
    except KeyError:
        print(["ERROR not such variable"])
        return defValue
