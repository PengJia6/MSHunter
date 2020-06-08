# =============================================================================
# Project : MShunter0.0.1
# Py Name: ScanMicrosatellites
# Author : units.py
# Date : 20-05-25
# Email : pengjia@stu.xjtu.edu.cn
# Description : ''
# =============================================================================
import numpy as np
import time


def removeZeroDict(dict):
    newdict = {}
    for key, value in dict.items():
        if value > 0.000001:
            newdict[key] = value
    return newdict


# the fastest
def getDisdistance(dict1, dict2):
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


def getDisdistance2(dict1, dict2):
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


def getDisdistance3(dict1, dict2):
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
    return np.sqrt(np.sum(np.square(np.array(list1) - np.array(list2))))


if __name__ == "__main__":
    ""
    dict1 = {1: 1, 2: 2, 3: 3}
    dict2 = {1: 2, 2: 4, 3: 8}
    start = time.time()
    step = 90000
    for i in range(step):
        getDisdistance(dict2, dict1)
    end = time.time()
    print(end - start)
    start = time.time()
    for i in range(step):
        getDisdistance2(dict2, dict1)
    end = time.time()
    print(end - start)
    start = time.time()
    for i in range(step):
        getDisdistance3(dict2, dict1)
    end = time.time()
    print(end - start)
    print(getDisdistance(dict1, dict2))
    print(getDisdistance2(dict1, dict2))
    print(getDisdistance3(dict1, dict2))
    # getDisdistance(dict2,dict1)
