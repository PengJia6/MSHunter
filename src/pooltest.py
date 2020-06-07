import multiprocess
import pysam
import multiprocessing as mp
import time

# class My:
#     def __init__(self, a, b):
#         self.a = a
#         # self.b = b
def run(My):
    # print(My.a)
    # print(My.b)
    # print(My)
    # print("")
    # i=0
    # while True:
    #     i+=1
    #     if i >1000000000:
    #         break
    # # time.sleep(100)
    return My


list0=[]
for i in range(5):
    if i==2:
        continue
    list0.append(My(i,i+1))

disvcfpath = "/mnt/project/ProjectMSI/MSCalling/note/py/test/NA12878.final.1000GDeep.dis.bcf"
vcffile = pysam.VariantFile(disvcfpath, "rb")
list1 = []
inum = 0
for i in vcffile.fetch():
    list1.append(i)
    inum += 1
    if inum > 100:
        break
pool=mp.Pool(3)
pool.map(run,list0)
pool.close()
pool.join()
print("Finished!")

#
# pool = multiprocess.Pool(2)
# res_l=[]
# for rec in list1:
#     res=pool.apply_async(run,args=(rec,))
#     res_l.append(res)

# res = pool.map(run, list1)
# pool.close()
# pool.join()
# for res in res_l:
#     print(res.get())
#
# results = [pool.apply_async(run, args=(name,)) for name in list1]
# results = [p.get() for p in results]
#
# end_t = datetime.datetime.now()
# elapsed_sec = (end_t - start_t).total_seconds()
#     print("多进程计算 共消耗: " + "{:.2f}".format(elapsed_sec) + " 秒")

#
# class MyThread(threading.Thread):
#
#     def __init__(self, func, args=()):
#         super(MyThread, self).__init__()
#         self.func = func
#         self.args = args
#
#     def run(self):
#         self.result = self.func(*self.args)
#
#     def get_result(self):
#         try:
#             return self.result  # 如果子线程不使用join方法，此处可能会报没有self.result的错误
#         except Exception:
#             return None
#
# threadList=[]
# for i in vcffile.fetch():
#     list1.append(i)
#     inum += 1
#     t = MyThread(run,[i])
#     threadList.append(t)
#     if inum > 100:
#         break
# print(len(threadList))
# i=1
# for t in threadList:
#     i+=1
#     print(i)
#     t.start()
# for t in threadList:
#     t.join()
