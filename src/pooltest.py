import multiprocessing

class My:
    def __init__(self,a,b):
        self.a=a
        self.b = b
def run(My):
    print(My.a)
    print(My.b)
    return My

list1=[]
for i in range(5):
    if i==2:
        continue
    list1.append(My(i,i+1))
for j in range(3):
    pool=multiprocessing.Pool(2)
    res=pool.map(run,list1)
    pool.close()
    # pool.join()
    for i in res:
        print(i.a)
