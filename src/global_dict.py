# =============================================================================
# Project : mstools-0.0.1
# Py Name: global_dict
# Author :
# Date : 20-05-26
# Email : pengjia@stu.xjtu.edu.cn
# Description : 'For global variable'
# =============================================================================


def global_init():
    global _global_dict
    _global_dict = {}
    _global_dict["tools_version"]="1.0"
    _global_dict["tools_name"]="mstools"
    _global_dict["chrom_list"] = [str(i) for i in range(1, 23)] + \
                                 ["chr" + str(i) for i in range(1, 23)] + \
                                 ["chrM", "chrX", "chrY", "MT", "X", "Y"]
    _global_dict["default"]={
        "bam2dis":{
            "threads":4,
            "minimum_mapping_quality":1,
            "minimum_support_reads":5,
            "batch":2000,
            "debug":False,
            "separator":"comma"

        },
        "errEval":{

        },
        "call":{

        }

    }

def set_value(name, value):
    _global_dict[name] = value


def get_value(name, defValue=None):
    try:
        return _global_dict[name]
    except KeyError:
        print(["ERROR not such variable"])
        return defValue
