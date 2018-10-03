from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtPrintSupport import *
from PyQt5.QtWidgets import *
from wordsegment import load, segment
import os
import sys
import uuid
import re
import csv
import math
import random
import matplotlib.pyplot as plt

FONT_SIZES = [7, 8, 9, 10, 11, 12, 13, 14, 18, 24, 36, 48, 64, 72, 96, 144, 288]
IMAGE_EXTENSIONS = ['.jpg','.png','.bmp']
HTML_EXTENSIONS = ['.htm', '.html']
CodonList = []  # CodonList 密码子数组
TerminationCodonList = []  # TerminationCodonList 终止密码子数组
LetterList = []  # 字母数组
CodonList2 = []
MinNumber = 0
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

cw = 0      # 当前总频率
best = []  # 最好组合
ps = []  # 当前组合
val = 100
sum = 0
new = []

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


# 匹配出与字母对应的密码子
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
    flag = [0, 0]
    for i in al:
        if (i == 'E'):
            flag[0] = 1
        if (i == 'I'):
            flag[1] = 1
    for i in al:
        for line in ls:
            if line[0] == i:
                if line[2][-1] == 0:
                    break
                tempCondon = getCodon(line[2])
                if (i == 'E'):
                    tempCondon = 'AAC'
                if (i == 'I'):
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
    i = random.randint(1, 2)
    if (i == 1):
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
    new1 = codon[len(mei1) + len(suiji) + len(_5UTR) + 3: -(len(_3UTR) + len(wei) + len(mei2))]
    # print('\n密码子序列为：')
    # print(new1)
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
    line = new[len(mei1): -len(mei2)]
    print('\n检测密码子序列中是否有酶序列:')
    print('序列:')
    print(line)
    # line = ''
    flag_mei1 = 0
    flag_mei2 = 0
    matchobj = re.match(r"(.*)GAATTC(.*)", line)
    if (matchobj == None):
        flag_mei1 = 0
        print('没有酶1序列')
    else:
        flag_mei2 = 1
        print('有酶1序列')

    matchobj = re.match(r"(.*)GCTAGC(.*)", line)
    if (matchobj == None):
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
    line = new[len(mei1): -len(mei2)]
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

    matchobj = re.match(r"(.*)GUAUGU(.*)UACUAAC(.*)CAG(.*)GUAUGU(.*)UACUAAC(.*)CAG(.*)", new)
    if (matchobj == None):
        print("没有内含子(长）序列")
        flag_long = 1
    else:
        print("有内含子（长）序列")

        temp_neihanzi1 = '...内含子GUAUGU' + matchobj.group(2) + 'UACUAAC' + 'CAG...'
        temp_neihanzi2 = '...内含子GUAUGU' + matchobj.group(5) + 'UACUAAC' + 'CAG...'
        if ((temp_neihanzi1 == neihanzi1 or temp_neihanzi1 == neihanzi2) and (
                temp_neihanzi2 == neihanzi2 or temp_neihanzi2 == neihanzi1)):
            flag_long = 0
            print('查询到插入的两个内含子：')
            print('其1：')
            print(matchobj.group(2))
            print('其2：')
            print(matchobj.group(5))
        else:
            flag_long = 1

        # 查询序列中是否有其他内含子，内含子歧义的情况
        if (flag_long == 0):
            temp_check = matchobj.group(1)
            temp_matchobj = re.match(r"(.*)GUAUGU(.*)", temp_check)
            if (temp_matchobj == None):
                flag_long = 0
            else:
                flag_long = 1
                print('有A-A-B-C结构序列')

        if (flag_long == 0):
            temp_check = matchobj.group(3)
            temp_matchobj = re.match(r"(.*)CAG(.*)", temp_check)
            if (temp_matchobj == None):
                flag_long = 0
            else:
                flag_long = 1
                print('有A-B-C-C结构序列')

        if (flag_long == 0):
            temp_check = matchobj.group(4)
            temp_matchobj = re.match(r"(.*)UACUAAC(.*)", temp_check)
            if (temp_matchobj == None):
                flag_long = 0
            else:
                flag_long = 1
                print('有A-B-C-B-A-B-C结构序列')

        if (flag_long == 0):
            temp_check = matchobj.group(4)
            temp_matchobj = re.match(r"(.*)GUAUGU(.*)", temp_check)
            if (temp_matchobj == None):
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

        if (flag_long == 0):
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

def hexuuid():
    return uuid.uuid4().hex

def splitext(p):
    return os.path.splitext(p)[1].lower()

class TextEdit(QTextEdit):

    def canInsertFromMimeData(self, source):

        if source.hasImage():
            return True
        else:
            return super(TextEdit, self).canInsertFromMimeData(source)

    def insertFromMimeData(self, source):

        cursor = self.textCursor()
        document = self.document()

        if source.hasUrls():

            for u in source.urls():
                file_ext = splitext(str(u.toLocalFile()))
                if u.isLocalFile() and file_ext in IMAGE_EXTENSIONS:
                    image = QImage(u.toLocalFile())
                    document.addResource(QTextDocument.ImageResource, u, image)
                    cursor.insertImage(u.toLocalFile())

                else:
                    break

            else:
                return


        elif source.hasImage():
            image = source.imageData()
            uuid = hexuuid()
            document.addResource(QTextDocument.ImageResource, uuid, image)
            cursor.insertImage(uuid)
            return

        super(TextEdit, self).insertFromMimeData(source)


class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        layout = QVBoxLayout()
        hlayout = QHBoxLayout()
        vlayout = QVBoxLayout()


        self.lettertocodon = QRadioButton("letter to codon")
        self.lettertocodon.setChecked(True)
        self.lettertocodon.toggled.connect(lambda:self.btnstate(self.lettertocodon))

        self.codontoletter = QRadioButton("codon to letter")
        self.codontoletter.setChecked(False)
        self.codontoletter.toggled.connect(lambda: self.btnstate(self.codontoletter))

        self.labl1 = QLabel("please input the letter:")

        self.editor1 = TextEdit()
        self.editor1.setAcceptRichText(False)
        self.editor1.setAutoFormatting(QTextEdit.AutoAll)
        self.editor1.selectionChanged.connect(self.update_format)
        font = QFont('Times', 12)
        self.editor1.setFont(font)
        self.editor1.setFontPointSize(12)

        self.editor2 = TextEdit()
        self.editor2.setAcceptRichText(False)
        self.editor2.setAutoFormatting(QTextEdit.AutoAll)
        self.editor2.selectionChanged.connect(self.update_format)
        font = QFont('Times', 12)
        self.editor2.setFont(font)
        self.editor2.setFontPointSize(12)

        self.button_translate = QPushButton('translate')
        self.button_translate.clicked.connect(self.onButtonClick)

        self.path = None

        hlayout.addWidget(self.lettertocodon)
        hlayout.addWidget(self.codontoletter)
        vlayout.addWidget(self.labl1)
        vlayout.addWidget(self.editor1)
        vlayout.addWidget(self.button_translate)
        vlayout.addWidget(self.editor2)

        hwg = QWidget()
        vwg = QWidget()

        hwg.setLayout(hlayout)
        vwg.setLayout(vlayout)

        layout.addWidget(hwg)
        layout.addWidget(vwg)

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

        self.status = QStatusBar()
        self.setStatusBar(self.status)

        file_toolbar = QToolBar("File")
        file_toolbar.setIconSize(QSize(14, 14))
        self.addToolBar(file_toolbar)
        file_menu = self.menuBar().addMenu("&File")

        open_file_action = QAction(QIcon(os.path.join('images', 'blue-folder-open-document.png')), "Open file...", self)
        open_file_action.setStatusTip("Open file")
        open_file_action.triggered.connect(self.file_open)
        file_menu.addAction(open_file_action)
        file_toolbar.addAction(open_file_action)

        save_file_action = QAction(QIcon(os.path.join('images', 'disk.png')), "Save", self)
        save_file_action.setStatusTip("Save current page")
        save_file_action.triggered.connect(self.file_save)
        file_menu.addAction(save_file_action)
        file_toolbar.addAction(save_file_action)

        saveas_file_action = QAction(QIcon(os.path.join('images', 'disk--pencil.png')), "Save As...", self)
        saveas_file_action.setStatusTip("Save current page to specified file")
        saveas_file_action.triggered.connect(self.file_saveas)
        file_menu.addAction(saveas_file_action)
        file_toolbar.addAction(saveas_file_action)

        print_action = QAction(QIcon(os.path.join('images', 'printer.png')), "Print...", self)
        print_action.setStatusTip("Print current page")
        print_action.triggered.connect(self.file_print)
        file_menu.addAction(print_action)
        file_toolbar.addAction(print_action)

        edit_toolbar = QToolBar("Edit")
        edit_toolbar.setIconSize(QSize(16, 16))
        self.addToolBar(edit_toolbar)
        edit_menu = self.menuBar().addMenu("&Edit")

        undo_action = QAction(QIcon(os.path.join('images', 'arrow-curve-180-left.png')), "Undo", self)
        undo_action.setStatusTip("Undo last change")
        undo_action.triggered.connect(self.editor1.undo)
        edit_menu.addAction(undo_action)

        redo_action = QAction(QIcon(os.path.join('images', 'arrow-curve.png')), "Redo", self)
        redo_action.setStatusTip("Redo last change")
        redo_action.triggered.connect(self.editor1.redo)
        edit_toolbar.addAction(redo_action)
        edit_menu.addAction(redo_action)

        edit_menu.addSeparator()

        cut_action = QAction(QIcon(os.path.join('images', 'scissors.png')), "Cut", self)
        cut_action.setStatusTip("Cut selected text")
        cut_action.setShortcut(QKeySequence.Cut)
        cut_action.triggered.connect(self.editor1.cut)
        edit_toolbar.addAction(cut_action)
        edit_menu.addAction(cut_action)

        copy_action = QAction(QIcon(os.path.join('images', 'document-copy.png')), "Copy", self)
        copy_action.setStatusTip("Copy selected text")
        cut_action.setShortcut(QKeySequence.Copy)
        copy_action.triggered.connect(self.editor1.copy)
        edit_toolbar.addAction(copy_action)
        edit_menu.addAction(copy_action)

        paste_action = QAction(QIcon(os.path.join('images', 'clipboard-paste-document-text.png')), "Paste", self)
        paste_action.setStatusTip("Paste from clipboard")
        cut_action.setShortcut(QKeySequence.Paste)
        paste_action.triggered.connect(self.editor1.paste)
        edit_toolbar.addAction(paste_action)
        edit_menu.addAction(paste_action)

        select_action = QAction(QIcon(os.path.join('images', 'selection-input.png')), "Select all", self)
        select_action.setStatusTip("Select all text")
        cut_action.setShortcut(QKeySequence.SelectAll)
        select_action.triggered.connect(self.editor1.selectAll)
        edit_menu.addAction(select_action)

        edit_menu.addSeparator()

        wrap_action = QAction(QIcon(os.path.join('images', 'arrow-continue.png')), "Wrap text to window", self)
        wrap_action.setStatusTip("Toggle wrap text to window")
        wrap_action.setCheckable(True)
        wrap_action.setChecked(True)
        wrap_action.triggered.connect(self.edit_toggle_wrap)
        edit_menu.addAction(wrap_action)

        format_toolbar = QToolBar("Format")
        format_toolbar.setIconSize(QSize(16, 16))
        self.addToolBar(format_toolbar)
        format_menu = self.menuBar().addMenu("&Format")

        self.fonts = QFontComboBox()
        self.fonts.currentFontChanged.connect(self.editor1.setCurrentFont)
        format_toolbar.addWidget(self.fonts)

        self.fontsize = QComboBox()
        self.fontsize.addItems([str(s) for s in FONT_SIZES])

        self.fontsize.currentIndexChanged[str].connect(lambda s: self.editor1.setFontPointSize(float(s)) )
        format_toolbar.addWidget(self.fontsize)

        self.bold_action = QAction(QIcon(os.path.join('images', 'edit-bold.png')), "Bold", self)
        self.bold_action.setStatusTip("Bold")
        self.bold_action.setShortcut(QKeySequence.Bold)
        self.bold_action.setCheckable(True)
        self.bold_action.toggled.connect(lambda x: self.editor1.setFontWeight(QFont.Bold if x else QFont.Normal))
        format_toolbar.addAction(self.bold_action)
        format_menu.addAction(self.bold_action)

        self.italic_action = QAction(QIcon(os.path.join('images', 'edit-italic.png')), "Italic", self)
        self.italic_action.setStatusTip("Italic")
        self.italic_action.setShortcut(QKeySequence.Italic)
        self.italic_action.setCheckable(True)
        self.italic_action.toggled.connect(self.editor1.setFontItalic)
        format_toolbar.addAction(self.italic_action)
        format_menu.addAction(self.italic_action)

        self.underline_action = QAction(QIcon(os.path.join('images', 'edit-underline.png')), "Underline", self)
        self.underline_action.setStatusTip("Underline")
        self.underline_action.setShortcut(QKeySequence.Underline)
        self.underline_action.setCheckable(True)
        self.underline_action.toggled.connect(self.editor1.setFontUnderline)
        format_toolbar.addAction(self.underline_action)
        format_menu.addAction(self.underline_action)

        format_menu.addSeparator()

        self.alignl_action = QAction(QIcon(os.path.join('images', 'edit-alignment.png')), "Align left", self)
        self.alignl_action.setStatusTip("Align text left")
        self.alignl_action.setCheckable(True)
        self.alignl_action.triggered.connect(lambda: self.editor1.setAlignment(Qt.AlignLeft))
        format_toolbar.addAction(self.alignl_action)
        format_menu.addAction(self.alignl_action)

        self.alignc_action = QAction(QIcon(os.path.join('images', 'edit-alignment-center.png')), "Align center", self)
        self.alignc_action.setStatusTip("Align text center")
        self.alignc_action.setCheckable(True)
        self.alignc_action.triggered.connect(lambda: self.editor1.setAlignment(Qt.AlignCenter))
        format_toolbar.addAction(self.alignc_action)
        format_menu.addAction(self.alignc_action)

        self.alignr_action = QAction(QIcon(os.path.join('images', 'edit-alignment-right.png')), "Align right", self)
        self.alignr_action.setStatusTip("Align text right")
        self.alignr_action.setCheckable(True)
        self.alignr_action.triggered.connect(lambda: self.editor1.setAlignment(Qt.AlignRight))
        format_toolbar.addAction(self.alignr_action)
        format_menu.addAction(self.alignr_action)

        self.alignj_action = QAction(QIcon(os.path.join('images', 'edit-alignment-justify.png')), "Justify", self)
        self.alignj_action.setStatusTip("Justify text")
        self.alignj_action.setCheckable(True)
        self.alignj_action.triggered.connect(lambda: self.editor1.setAlignment(Qt.AlignJustify))
        format_toolbar.addAction(self.alignj_action)
        format_menu.addAction(self.alignj_action)

        format_group = QActionGroup(self)
        format_group.setExclusive(True)
        format_group.addAction(self.alignl_action)
        format_group.addAction(self.alignc_action)
        format_group.addAction(self.alignr_action)
        format_group.addAction(self.alignj_action)

        format_menu.addSeparator()

        self._format_actions = [
            self.fonts,
            self.fontsize,
            self.bold_action,
            self.italic_action,
            self.underline_action,
        ]

        self.flag = 1
        # Initialize.
        self.update_format()
        self.update_title()
        self.show()

    def onButtonClick(self):
        sender = self.sender()
        print(sender.text() + ' 被按下了')
        if(self.flag == 1):
            al = self.editor1.toPlainText()
            al = al.upper()
            print(al)
            cnt = 0
            while (1):
                cnt += 1
                print('第%d次生成序列' % cnt)
                new = []
                new = LetterToCodon(LetterList, al, neihanzi1, neihanzi2, TerminationCodonList, n, mei1, suiji, _5UTR,
                                    _3UTR, wei, mei2)
                print(new)
                flag_mei1, flag_mei2 = checkEnzyme(new, mei1, mei2)
                if (flag_mei1 == 1 or flag_mei2 == 1):
                    continue
                flag_intron = checkIntron(new, mei1, mei2, neihanzi1, neihanzi2)
                if (flag_intron == 1):
                    continue
                break
            print('\n-------------------------------- checking over --------------------------------')
            text = new
            self.editor2.setText(text)
        else:
            codon = self.editor1.toPlainText()
            #print(codon)
            codon1 = dealCodon(codon, mei1, suiji, _5UTR, TerminationCodonList, neihanzi1, neihanzi2, _3UTR, wei, mei2)
            #print(codon1)
            le1 = ''
            le1 = CodonToLetter(LetterList, codon1)
            print(le1)
            self.editor2.setText(le1)
            #print(LetterList)
            #load()
            #text = segment(le1)
            #text2 = ''
            #for word in text:
                #text2 = text2+word+' '
            #print(segment(text))
            #text2 = text2.upper()
            #print(text2)
            #self.editor2.setText(text2)

    def btnstate(self, btn):
        if btn.text() == "letter to codon":
            if btn.isChecked() == True:
                print(btn.text() + " is selected")
                self.labl1.setText('please input the letter:')
                self.flag = 1
            else:
                print(btn.text() + " is deselected")

        if btn.text() == "codon to letter":
            if btn.isChecked() == True:
                print(btn.text() + " is selected")
                self.labl1.setText('please input the condon:')
                self.flag = 2
            else:
                print(btn.text() + " is deselected")

    def block_signals(self, objects, b):
        for o in objects:
            o.blockSignals(b)

    def update_format(self):
        self.block_signals(self._format_actions, True)

        self.fonts.setCurrentFont(self.editor1.currentFont())
        self.fontsize.setCurrentText(str(int(self.editor1.fontPointSize())))

        self.italic_action.setChecked(self.editor1.fontItalic())
        self.underline_action.setChecked(self.editor1.fontUnderline())
        self.bold_action.setChecked(self.editor1.fontWeight() == QFont.Bold)

        self.alignl_action.setChecked(self.editor1.alignment() == Qt.AlignLeft)
        self.alignc_action.setChecked(self.editor1.alignment() == Qt.AlignCenter)
        self.alignr_action.setChecked(self.editor1.alignment() == Qt.AlignRight)
        self.alignj_action.setChecked(self.editor1.alignment() == Qt.AlignJustify)

        self.block_signals(self._format_actions, False)

    def dialog_critical(self, s):
        dlg = QMessageBox(self)
        dlg.setText(s)
        dlg.setIcon(QMessageBox.Critical)
        dlg.show()

    def file_open(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open file", "", "HTML documents (*.html);Text documents (*.txt);All files (*.*)")

        try:
            with open(path, 'rU') as f:
                text = f.read()

        except Exception as e:
            self.dialog_critical(str(e))

        else:
            self.path = path
            self.editor1.setText(text)
            self.update_title()

    def file_save(self):
        if self.path is None:
            return self.file_saveas()

        text = self.editor1.toHtml() if splitext(self.path) in HTML_EXTENSIONS else self.editor1.toPlainText()

        try:
            with open(self.path, 'w') as f:
                f.write(text)

        except Exception as e:
            self.dialog_critical(str(e))

    def file_saveas(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save file", "", "HTML documents (*.html);Text documents (*.txt);All files (*.*)")

        if not path:
            return

        text = self.editor1.toHtml() if splitext(path) in HTML_EXTENSIONS else self.editor1.toPlainText()

        try:
            with open(path, 'w') as f:
                f.write(text)

        except Exception as e:
            self.dialog_critical(str(e))

        else:
            self.path = path
            self.update_title()

    def file_print(self):
        dlg = QPrintDialog()
        if dlg.exec_():
            self.editor1.print_(dlg.printer())

    def update_title(self):
        self.setWindowTitle("%s" % (os.path.basename(self.path) if self.path else "Biological Encryption"))

    def edit_toggle_wrap(self):
        self.editor1.setLineWrapMode( 1 if self.editor1.lineWrapMode() == 0 else 0 )


CodonList = []  # CodonList 密码子数组
TerminationCodonList = []  # TerminationCodonList 终止密码子数组
CodonList = read("mima.csv")
LenOfCodonList = len(CodonList)  # lenOfCodonList=64 64个密码子

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

LetterList = []  # 字母数组
LetterList = read("zimu.csv")
CodonList2 = CodonList.copy()

# 将频率最大的密码子分配给频率最小的几个字母
# 1 GAU 4.06 -- E(3.99)+Z(0.07)
# 2 GCU 3.47 -- E(3.37)+J(0.10)
# 3 GAA 3.00 -- I(2.89)+Q(0.11)
# 4 AUG 2.89 -- I(2.72)+X(0.17)
MinNumber = 0
MinNumber = allotMin(CodonList, LetterList)
for k in range(MinNumber, 26):
    cw = 0  # 当前总频率
    best = []  # 最好组合
    ps = []  # 当前组合
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
# out(LetterList)
t = sorted(LetterList, key=lambda LetterList: LetterList[0], reverse=False)
# drawing(t)


# 字母转密码子
getProbability(LetterList)
if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    app.exec_()
