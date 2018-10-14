import re
import csv
import math
import random
import matplotlib.pyplot as plt

def read(name):
    with open(name, "r", encoding="utf-8-sig") as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        ls = []
        for col in reader:
            ls.append(col)
        # print(ls)
        for i in range(len(ls)):
            ls[i][1] = float(ls[i][1])
        ls1 = sorted(ls, key=lambda ls: ls[1], reverse=False)
    return ls1

def dealEnd(end):
    sum = 0
    for i in end:
        sum = sum + i[1]
    for i in end:
        p = round(i[1] / sum, 3) * 1000
        i.append(p)
    end.append(sum)
    # [['UAG', 0.0, 0.0], ['UGA', 0.0, 0.0], ['UAA', 0.11, 1000.0], 0.11]

def allotMin(ls, ls1):
    j = 60
    for i in range(0, len(ls1)):
        x = []
        if ls1[i][1] < ls[0][1]:
            x.insert(0, ls[j].copy())
            x[0][1] = ls1[i][1]
            x.insert(1, ls1[i][1])
            ls1[i].append(x)
            ls[j][1] = round(ls[j][1] - ls1[i][1], 2)
            j = j - 1
        else:
            return i

def backtrack(ls, p, i):
    # print(i)
    global cw, best, ps, val
    # print(best)
    if math.fabs(cw - p) == 0:
        best = ps.copy()
        val = 0
        return
    if i > len(ls):
        if math.fabs(cw - p) <= val:
            best = ps.copy()
            val = math.fabs(cw - p)
        return
    if cw + ls[i - 1][1] <= p:
        cw = cw + ls[i - 1][1]
        ps.append(ls[i - 1])
        backtrack(ls, p, i + 1)
        cw = cw - ls[i - 1][1]
        ps.remove(ls[i - 1])
    else:
        if math.fabs(cw - p) <= val:
            best = ps.copy()
            val = math.fabs(cw - p)
    if (val > 0.01):
        backtrack(ls, p, i + 1)


def out(ls):
    t = sorted(ls, key=lambda ls: ls[0], reverse=False)
    for line in t:
        print(line[0], "%.2f" % line[1], end="|  ")
        for i in range(len(line[2]) - 1):
            if line[2][0] != 0:
                print(line[2][i][0], line[2][i][1], end="  ")

        print("")
        print("sum:%.2lf" % line[2][-1])


def drawing(ls):
    t = ls
    num_list = []
    name_list = []
    num_list1 = []
    for line in t:
        name_list.append(line[0])
        num_list.append(line[1])
        num_list1.append(line[2][-1])
    x = list(range(len(num_list)))
    x1 = x.copy()
    plt.bar(x, num_list, width=0.4, label='letter', fc='r')
    for i in range(len(x)):
        x1[i] = x[i] + 0.4
    plt.bar(x1, num_list1, width=0.4, label='codon', tick_label=name_list, fc='b', align='center')

    plt.legend()
    plt.show()





CodonList = []
TerminationCodonList = []
CodonList = read("mima.csv")
LenOfCodonList = len(CodonList)

i = 0
while i < LenOfCodonList:
    if CodonList[i][0] == 'UAA' or CodonList[i][0] == 'UGA' or CodonList[i][0] == 'UAG':
        TerminationCodonList.append(CodonList.pop(i))
        i = i - 1
        LenOfCodonList = LenOfCodonList - 1
    i = i + 1

dealEnd(TerminationCodonList)

LetterList = []
LetterList = read("zimu.csv")
CodonList2 = CodonList.copy()


# 1 GAU 4.06 -- E(3.99)+Z(0.07)
# 2 GCU 3.47 -- E(3.37)+J(0.10)
# 3 GAA 3.00 -- I(2.89)+Q(0.11)
# 4 AUG 2.89 -- I(2.72)+X(0.17)
MinNumber = 0
MinNumber = allotMin(CodonList, LetterList)
for k in range(MinNumber, 26):
    cw = 0
    best = []
    ps = []
    val = 100
    sum = 0
    backtrack(CodonList2, LetterList[k][1], 1)
    for rol in best:
        sum = sum + rol[1]
        CodonList2.remove(rol)
    best.append(sum)
    LetterList[k].append(best)
out(LetterList)
t = sorted(LetterList, key=lambda LetterList: LetterList[0], reverse=False)
drawing(t)