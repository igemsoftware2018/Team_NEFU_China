import re
import csv
import math
import random
import matplotlib.pyplot as plt

#GAATTC随机序列TAGGTTGCTTCTTTTAGTGGTTTGCAAUGACU...内含子GUAUGUATGGATTCTGGTATGTTCTAGCGCTTGCACCATCCCATTTAACTGTAAGAAGAATTGCACGGTCCCAATTGCTCGAGAGATTTCTCTTTTACCTTTTTTTACTATTTTTCACTCTCCCATAACCTCCTATATTGACTGATCTGTAATAACCACGATATTATTGGAATAAATAGGGGCTTGAAATTTGGAAAAAAAAAAAAAACTGAAATATTTTCGTGATAAGTGATAGTGATATTCTTCTTTTATTTGCTACTGTTACTAAGTCTCATGTACTAACATCGATTGCTTCATTCTTTTTGTTGCTATATTATATGTTTAUACUAACCAG...ACG...内含子GUAUGUATGGGTAGAGTTAGAACCAAGACCGTCAAGCGTGCTTCTAAGGCTTTGATTGAACGTTACTATCCAAAGTTGACTTTGGATTTCCAAACCAACAAGAGACTTTGTGATGAAATCGCCACTATCCAATCCAAGAGATTGAGAAACAAGATTGCTGGTTACACCACCCATTTGATGAAGAGAATCCAAAAGGGTCCAGTTAGAGGTATCTCTTTCAAATTGCAAGAAGAAGAAAGAGAAAGAAAGGACCAATACGTCCCAGAAGTCTCTGCTTTGGACTTGTCTCGTTCTAACGGTGTTTTGAACGTTGACAACCAAACTTCTGACTTGGTTAAATCTTTGGGTTTGAAGTTGCCATTATCTGTTATCAACGTTTCTGCCCAAAGAGACAGACGTTACAGAAAGAGAGTTTAAUACUAACCAG...CGUAUAUAATTTTCGTCTCTTATTATTAAACCTTTAAAAACGCTATCCTTGACTTTATCTGTACTTTGCAATAAAAGCAGGCTCTGAGTGTTTAAATCTATTTTTCTTTCATTCGCTAGC
# 从mima.csv和zimu.csv里读出数据
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


# 处理三个终止密码子，计算相应出现概率以及总和
def dealEnd(end):

    sum = 0
    for i in end:
        sum = sum + i[1]
    for i in end:
        p = round(i[1] / sum, 3) * 1000
        i.append(p)
    end.append(sum)
    # [['UAG', 0.0, 0.0], ['UGA', 0.0, 0.0], ['UAA', 0.11, 1000.0], 0.11]
    

# 将频率最大的密码子分配给频率最小的几个字母
def allotMin(ls,ls1):
    j=60
    for i in range(0,len(ls1)):
        x = []
        if ls1[i][1]<ls[0][1]:
            x.insert(0,ls[j].copy())
            x[0][1]=ls1[i][1]
            x.insert(1,ls1[i][1])
            ls1[i].append(x)
            ls[j][1]=round(ls[j][1]-ls1[i][1],2)
            j=j-1
        else:
            return i

        
#匹配出与字母对应的密码子
def backtrack(ls, p, i):  # p:某字母的频率
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


# 输出密码子和字母的匹配结果
def out(ls):
    t = sorted(ls, key=lambda ls: ls[0], reverse=False)
    for line in t:
        print(line[0], "%.2f" % line[1], end="|  ")
        for i in range(len(line[2]) - 1):
            if line[2][0] != 0:
                print(line[2][i][0], line[2][i][1], end="  ")
        # print("方案为:", best)
        print("")
        print("sum:%.2lf" % line[2][-1])


# 把字母和密码子对应的柱状图画出来
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
    '''# x = np.arange(26)+1
    for z, y in zip(x, num_list):
        z1=str(round(z*0.35,2))
        y1=str(y+0.3)
        print(type(z))
        print(z)
        plt.text(z1, y1, '%.2f' % y, ha='center', va='bottom')

    for z, y in zip(x1, num_list1):
        #x1 = str(x)
        y1 = str(y)
        plt.text(z + 0.3, y1 + 0.05, '%.2f' % y1, ha='center', va='bottom')'''
    plt.legend()
    plt.show()


def getProbability(ls):
    for line in ls:
        for i in range(len(line[2]) - 1):
            p = round(line[2][i][1] / line[2][-1], 3) * 1000
            line[2][i].append(p)


def getCodon(ls):
    new1 = []
    for i in range(0, len(ls) - 1):
        for j in range(0, int(ls[i][2])):
            new1.append(i)
    x = random.randint(0, 1000)
    return ls[new1[x]][0]


# 字母转换成密码子
def LetterToCodon(ls, al, neihanzi1, neihanzi2, end, n, mei1, suiji, _5UTR, _3UTR, wei, mei2):
    # al:字母序列
    new = []
    flag = [0,0]
    for i in al:
        if(i == 'E'):
            flag[0] = 1
        if(i == 'I'):
            flag[1] = 1
    for i in al:
        for line in ls:
            if line[0] == i:
                if line[2][-1] == 0:
                    break
                tempCondon = getCodon(line[2])
                if(i == 'E'):
                    tempCondon = 'AAC'
                if(i == 'I'):
                    tempCondon = 'GUC'
                new.append(tempCondon)
                break
    new.append(getCodon(end))
    i = 1
    dis = []  # 存储要放进内含子的位置
    x = random.randint(0, len(al))
    dis.append(x)
    while (1):
        flag = 0
        x1 = random.randint(0, len(al))
        '''
        for j in range(0, len(dis)):
            if math.fabs(dis[j] - x1) * 3 < 2 * len(neihanzi):
                flag = 1
                break
        '''
        if flag == 0:
            dis.append(x1)
            i = i + 1
            if i >= int(n):
                break
    i = random.randint(1,2)
    if(i == 1):
        new.insert(dis[0], neihanzi1)
        new.insert(dis[1], neihanzi2)
    else:
        new.insert(dis[0], neihanzi2)
        new.insert(dis[1], neihanzi1)
    '''
    for i in dis:
        new.insert(i, neihanzi
    '''
    new.insert(0, mei1)
    new.insert(1, suiji)
    new.insert(2, _5UTR)
    new.insert(3, 'AUG')
    new.append(_3UTR)
    new.append(wei)
    new.append(mei2)
    new1 = []
    new1 = "".join(new)
    return new1
    

# 将输入的密码子序列提取出中间有用的那段序列
def dealCodon(codon, mei1, suiji, _5UTR, end, neihanzi1, neihanzi2, _3UTR, wei, mei2):
    new1 = codon[len(mei1) + len(suiji) + len(_5UTR) + 3 : -(len(_3UTR) + len(wei) + len(mei2))]
    #print('\n密码子序列为：')
    #print(new1)
    new2 = ''
    new3 = []
    # print(new1)
    x = ''
    j = 0
    i = 0
    while i < (len(new1)):
        # for i in range(0,len(new1)-len()):
        if new1[i:i + len(neihanzi1)] == neihanzi1:
            # k=i+len(neihanzi)
            i = i + len(neihanzi1)
            # print('yes1')
        elif new1[i:i + len(neihanzi2)] == neihanzi2:
            i = i + len(neihanzi2)
            # print('yes2')
        else:
            # new2 = new2 + new1[k]
            # k=k+1
            new2 = new2 + new1[i]
            i = i + 1
    for i in new2:
        x = x + i
        j = j + 1
        flag = 0
        if j == 3:
            for line in range(3):
                if x == end[line][0]:
                    flag = 1
                    break
            if flag == 1:
                break
            j = 0
            new3.append(x)
            x = ''

    return new3


# 将密码子转换成字母序列
def CodonToLetter(ls, new):
    le = []
    for i in new:
        flag = 0
        for line in ls:
            for k in range(0, len(line[2]) - 1):
                if i == line[2][k][0]:
                    le.append(line[0])
                    flag = 1
                    break
            if flag == 1:
                break
    le1 = "".join(le)
    return le1


# 检测密码子序列中是否有酶序列
def checkEnzyme(new, mei1, mei2):
    line = new[len(mei1) : -len(mei2)]
    print('\n检测密码子序列中是否有酶序列:')
    print('序列:')
    print(line)
    # line = ''
    flag_mei1 = 0
    flag_mei2 = 0
    matchobj=re.match(r"(.*)GAATTC(.*)", line)
    if(matchobj == None):
        flag_mei1 = 0
        print('没有酶1序列')
    else:
        flag_mei2 = 1
        print('有酶1序列')
    
    matchobj=re.match(r"(.*)GCTAGC(.*)", line)
    if(matchobj == None):
        flag_mei2 = 0
        print('没有酶2序列')
    else:
        flag_mei2 = 1
        print('有酶2序列')
    return flag_mei1, flag_mei2


# 检测密码子序列中是否有内含子序列
def checkIntron(new, mei1, mei2, neiahnzi1, neihanzi2):
    # line=" AAAGUAUGUCCCUACUAACCAGTTT"
    flag_long = 0
    # flag_short = 0
    line = new[len(mei1) : -len(mei2)]
    print('\n检测密码子序列中是否有内含子序列（长）:')
    # 先检测是否有长序列
    '''
    j = [0]
    k = [0]
    for i in range(0, len(line)):
        new_line = line[j[0] : i]
        #print(new_line)
        matchobj=re.match(r"(.*)GUAUGU(.*)UACUAAC(.*)CAG(.*)", new_line)
        if(matchobj == None):
            continue
        else:
            print('由内含子') 
            print(matchobj.group(2))
            
            break
        if(i == 10):
            break
    '''
        
    matchobj=re.match(r"(.*)GUAUGU(.*)UACUAAC(.*)CAG(.*)GUAUGU(.*)UACUAAC(.*)CAG(.*)",new)
    if(matchobj == None):              
        print("没有内含子(长）序列")
        flag_long = 1
    else:
        print("有内含子（长）序列")
        
        temp_neihanzi1 = '...内含子GUAUGU' + matchobj.group(2) + 'UACUAAC' + 'CAG...'
        temp_neihanzi2 = '...内含子GUAUGU' + matchobj.group(5) + 'UACUAAC' + 'CAG...'
        if((temp_neihanzi1 == neihanzi1 or temp_neihanzi1 == neihanzi2) and (temp_neihanzi2 == neihanzi2 or temp_neihanzi2 == neihanzi1)):
            flag_long = 0
            print('查询到插入的两个内含子：')
            print('其1：')
            print(matchobj.group(2))
            print('其2：')
            print(matchobj.group(5))
        else:
            flag_long = 1

        # 查询序列中是否有其他内含子，内含子歧义的情况
        if(flag_long == 0):
            temp_check = matchobj.group(1)
            temp_matchobj=re.match(r"(.*)GUAUGU(.*)", temp_check)
            if(temp_matchobj == None):
                flag_long = 0
            else:
                flag_long = 1
                print('有A-A-B-C结构序列')
                

        if(flag_long == 0):
            temp_check = matchobj.group(3)
            temp_matchobj=re.match(r"(.*)CAG(.*)", temp_check)
            if(temp_matchobj == None):
                flag_long = 0
            else:
                flag_long = 1
                print('有A-B-C-C结构序列')

        if(flag_long == 0):
            temp_check = matchobj.group(4)
            temp_matchobj=re.match(r"(.*)UACUAAC(.*)", temp_check)
            if(temp_matchobj == None):
                flag_long = 0
            else:
                flag_long = 1
                print('有A-B-C-B-A-B-C结构序列')

        if(flag_long == 0):
            temp_check = matchobj.group(4)
            temp_matchobj=re.match(r"(.*)GUAUGU(.*)", temp_check)
            if(temp_matchobj == None):
                flag_long = 0
            else:
                flag_long = 1
                print('有A-2A-B-C结构序列')
        '''
        if(flag_long == 0):
            temp_check = matchobj.group(6)
            temp_matchobj=re.match(r"(.*)CAG(.*)", temp_check)
            if(temp_matchobj == None):
                flag_long = 0
            else:
                flag_long = 1
                print('有2A-B-C-C结构序列')
        '''

        if(flag_long == 0):
            print('序列没有其他内含子')
        else:
            print('序列有其他内含子')
    return flag_long
    
    '''
    print(matchobj.group(0)) 
    print(matchobj.group(1))  
    print(matchobj.group(2))  
    print(matchobj.group(3))
    '''
    


        
CodonList = []                            # CodonList 密码子数组
TerminationCodonList = []                 # TerminationCodonList 终止密码子数组
CodonList = read("mima.csv")              
LenOfCodonList = len(CodonList)           # lenOfCodonList=64 64个密码子
        
# 密码子数组中删除终止密码子'UAA' 'UGA' 'UAG',存储于终止密码子数组
i = 0
while i < LenOfCodonList:
    if CodonList[i][0] == 'UAA' or CodonList[i][0] == 'UGA' or CodonList[i][0] == 'UAG':
        TerminationCodonList.append(CodonList.pop(i))  
        i = i - 1
        LenOfCodonList = LenOfCodonList - 1
    i = i + 1
            
# 处理三个终止密码子，计算相应出现概率以及总和    
dealEnd(TerminationCodonList)

LetterList = [] # 字母数组
LetterList = read("zimu.csv") 
CodonList2 = CodonList.copy()
        
# 将频率最大的密码子分配给频率最小的几个字母 
# 1 GAU 4.06 -- E(3.99)+Z(0.07)
# 2 GCU 3.47 -- E(3.37)+J(0.10)
# 3 GAA 3.00 -- I(2.89)+Q(0.11)
# 4 AUG 2.89 -- I(2.72)+X(0.17)
MinNumber=0
MinNumber=allotMin(CodonList,LetterList)
for k in range(MinNumber, 26):
    cw = 0     # 当前总频率
    best = []  # 最好组合
    ps = []    # 当前组合
    val = 100
    sum = 0
    backtrack(CodonList2, LetterList[k][1], 1)
    for rol in best:
        sum = sum + rol[1]
        CodonList2.remove(rol)
    '''print(ls1[k][0], "%.2f"%ls1[k][1], end="|  ")
    for line in best:
        print(line[0],line[1],end="  ")
    print("方案为:", best)
    print("")
    print("sum:%.2lf" % sum)'''
    best.append(sum)
    LetterList[k].append(best)
#out(LetterList)   
t= sorted(LetterList, key=lambda LetterList: LetterList[0], reverse=False)
#drawing(t)


# 字母转密码子
getProbability(LetterList)
while(1):
    al = input('please enter the letters:')
    # neihanzi = '0000000'  # input('please enter neihanzi:')
    #(.*)GUAUGU(.*)UACUAACCAG(.*)
    neihanzi1='...内含子GUAUGU'\
             +'ATGGGTAGAGTTAGAACCAAGACCGTCAAGCGTGCTTCTAAGGCTTTGATTGAACGTTACTATCCAAAGTTGACTTTG'\
             +'GATTTCCAAACCAACAAGAGACTTTGTGATGAAATCGCCACTATCCAATCCAAGAGATTGAGAAACAAGATTGCTGGT'\
             +'TACACCACCCATTTGATGAAGAGAATCCAAAAGGGTCCAGTTAGAGGTATCTCTTTCAAATTGCAAGAAGAAGAAAGA'\
             +'GAAAGAAAGGACCAATACGTCCCAGAAGTCTCTGCTTTGGACTTGTCTCGTTCTAACGGTGTTTTGAACGTTGACAAC'\
             +'CAAACTTCTGACTTGGTTAAATCTTTGGGTTTGAAGTTGCCATTATCTGTTATCAACGTTTCTGCCCAAAGAGACAGA'\
             +'CGTTACAGAAAGAGAGTTTAA'\
             +'UACUAACCAG...'
    neihanzi2='...内含子GUAUGU'\
             +'ATGGATTCTGGTATGTTCTAGCGCTTGCACCATCCCATTTAACTGTAAGAAGAATTGCACGGTCCCAATTGCTCGAGA'\
             +'GATTTCTCTTTTACCTTTTTTTACTATTTTTCACTCTCCCATAACCTCCTATATTGACTGATCTGTAATAACCACGAT'\
             +'ATTATTGGAATAAATAGGGGCTTGAAATTTGGAAAAAAAAAAAAAACTGAAATATTTTCGTGATAAGTGATAGTGATA'\
             +'TTCTTCTTTTATTTGCTACTGTTACTAAGTCTCATGTACTAACATCGATTGCTTCATTCTTTTTGTTGCTATATTATA'\
             +'TGTTTA'\
             +'UACUAACCAG...'

    n = 2                 # input('please enter the number of neihanzi:')
    mei1 = 'GAATTC'       # input('please enter mei1:')
    suiji = '随机序列'    # input('please enter suiji:')
    _5UTR = 'TAGGTTGCTTCTTTTAGTGGTTTGCA'
    _3UTR = 'TTTTCGTCTCTTATTATTAAACCTTTAAAAACGCTATCCTTGACTTTATCTGTACTTTGC'
    wei = 'AATAAAAGCAGGCTCTGAGTGTTTAAATCTATTTTTCTTTCATTC' # input('please enter wei:')
    mei2 = 'GCTAGC'       # input('please enter mei2:')
    
    cnt = 0
    while(1):
        cnt += 1
        print('第%d次生成序列' % cnt)
        new = []
        new = LetterToCodon(LetterList, al, neihanzi1, neihanzi2, TerminationCodonList, n, mei1, suiji, _5UTR, _3UTR, wei, mei2)
        print (new)
        flag_mei1, flag_mei2 = checkEnzyme(new, mei1, mei2)
        if(flag_mei1 == 1 or flag_mei2 == 1):
            continue
        flag_intron = checkIntron(new, mei1, mei2, neihanzi1, neihanzi2)
        if(flag_intron == 1):
            continue
        break
    print('\n-------------------------------- checking over --------------------------------')
    print('condon:')
    print(new)
    # 密码子转字母
    codon = input('please enter the codon:')
    print(codon)
    codon1 = []
    codon1 = dealCodon(codon, mei1, suiji, _5UTR, TerminationCodonList, neihanzi1, neihanzi2, _3UTR, wei, mei2)
    print (codon1)
    le1 = ''
    le1 = CodonToLetter(LetterList, codon1)
    print (le1)

    inputt = input('please enter 1 to continue, 0 to break:')
    if( inputt == '1' ):
        continue
    else:
        break
    

