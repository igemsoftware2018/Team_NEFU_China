import re
from math import log

FILE_NAME = 'dict.txt'

def getDitc(f):
    lfreq = {}
    ltotal = 0
    file = open(f, 'rb')
    for lineno, line in enumerate(file, 1):
        line = line.strip().decode('utf-8')
        #print('line '+str(lineno)+': '+line)
        regex = re.compile('\s+')
        word, freq = regex.split(line)[:2]
        freq = int(freq)
        if(word not in lfreq):
            lfreq[word] = freq
        else:
            lfreq[word] += freq
        ltotal += log(freq)
        for ch in range(len(word)):
            wfrag = word[:ch + 1]
            if wfrag not in lfreq:
                lfreq[wfrag] = 0
    file.close()
    #print(lfreq)
    return lfreq, ltotal

# DAG
def getDAG(s, freq):
    DAG = {}
    N = len(s)
    for k in range(N):
        tmplist = []
        i = k
        frag = s[k]
        while i < N and frag in freq:
            if freq[frag]:
                tmplist.append(i)
            i += 1
            frag = s[k:i + 1]
        if not tmplist:
            tmplist.append(k)
        DAG[k] = tmplist
    print(DAG)
    return DAG

# DAG
def cutAll(s, DAG):
    dog = DAG
    old_j = -1
    List = []
    for k, L in dog.items():
        if len(L) == 1 and k > old_j:
            List.append(s[k:L[0] + 1])
            old_j = L[0]
        else:
            for j in L:
                if j > k:
                    List.append(s[k:j + 1])
                    old_j = j
    return List


def getCalc(sentence, DAG, route, freq):
    N = len(sentence)
    route[N] = (0, 0)
    for idx in range(N - 1, -1, -1):
        max = -1e9
        flag = DAG[idx][0]
        for x in DAG[idx]:
            s = sentence[idx:x + 1]
            print("s=",s,",max=",max,",x=",x,",idx=",idx)
            if(s in freq):
                if(freq[s] == 0):
                    res = 0 + route[x + 1][0]
                else:
                    res = 1.0/freq[s] + route[x + 1][0]
            else:
                res = 0 + route[x + 1][0]
            if(res > max):
                max = res
                flag = x
        route[idx] = (max, flag)

    #for idx in range(N - 1, -1, -1):
        #route[idx] = max((log(freq[sentence[idx:x + 1]] or 1) + route[x + 1][0], x) for x in DAG[idx])
    print(route)
    x = 0
    buf = ''
    result = []
    #print('result:')
    while x < N:
        y = route[x][1] + 1
        l_word = sentence[x:y]
        if(y - x == 1):
            buf += l_word
        result.append(l_word)
        x = y
    #print(result)
    return result

def main():
    freq, total = getDitc(FILE_NAME)
    route = {}
    #sentence = 'Iwantogotoschool!'
    while 1:
         sentence = input('please input the sentences:')
         #sentence = 'afternoon'
         sentence = sentence.lower()
         dag = getDAG(sentence, freq)
         allResL = cutAll(sentence, dag)
         bestResL = getCalc(sentence, dag, route, freq)
         print('all results:', end = ' ')
         for word in allResL:
             print(word,  end = ' ')
         print('\nbest results:', end = ' ')
         for word in bestResL:
             print(word, end = ' ')
         end = input('\npress 1 to continue, any else to break:')
         if(end == '1'):
             continue
         else:
             break

main()
