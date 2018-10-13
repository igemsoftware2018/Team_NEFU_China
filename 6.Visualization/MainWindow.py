from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtPrintSupport import *
from PyQt5.QtWidgets import *
from PyQt5 import QtWidgets,QtGui
from wordsegment import load, segment
import os
import sys
import uuid
import re
import csv
import math
import random
import matplotlib.pyplot as plt

FONT_SIZES = [7, 8, 9, 10, 11, 12, 13, 14, 18, 24, 36, 48, 64, 72, 96, 144, ]
HTML_EXTENSIONS = ['.htm', '.html']
CodonList = []  # CodonList codons arrays
TerminationCodonList = []  # Termination arrays
LetterList = []  # Letters arrays
CodonList2 = []
MinNumber = 0
neihanzi1 = 'GUAUGU' \
            + 'ATGGGTAGAGTTAGAACCAAGACCGTCAAGCGTGCTTCTAAGGCTTTGATTGAACGTTACTATCCAAAGTTGACTTTG' \
            + 'GATTTCCAAACCAACAAGAGACTTTGTGATGAAATCGCCACTATCCAATCCAAGAGATTGAGAAACAAGATTGCTGGT' \
            + 'TACACCACCCATTTGATGAAGAGAATCCAAAAGGGTCCAGTTAGAGGTATCTCTTTCAAATTGCAAGAAGAAGAAAGA' \
            + 'GAAAGAAAGGACCAATACGTCCCAGAAGTCTCTGCTTTGGACTTGTCTCGTTCTAACGGTGTTTTGAACGTTGACAAC' \
            + 'CAAACTTCTGACTTGGTTAAATCTTTGGGTTTGAAGTTGCCATTATCTGTTATCAACGTTTCTGCCCAAAGAGACAGA' \
            + 'CGTTACAGAAAGAGAGTTTAA' \
            + 'UACUAACCAG'
neihanzi2 = 'GUAUGU' \
            + 'ATGGATTCTGGTATGTTCTAGCGCTTGCACCATCCCATTTAACTGTAAGAAGAATTGCACGGTCCCAATTGCTCGAGA' \
            + 'GATTTCTCTTTTACCTTTTTTTACTATTTTTCACTCTCCCATAACCTCCTATATTGACTGATCTGTAATAACCACGAT' \
            + 'ATTATTGGAATAAATAGGGGCTTGAAATTTGGAAAAAAAAAAAAAACTGAAATATTTTCGTGATAAGTGATAGTGATA' \
            + 'TTCTTCTTTTATTTGCTACTGTTACTAAGTCTCATGTACTAACATCGATTGCTTCATTCTTTTTGTTGCTATATTATA' \
            + 'TGTTTA' \
            + 'UACUAACCAG'

n = 2  # numbers of introns
mei1 = 'GAATTC'  # enzyme 1
suiji = ''  # random serious
_5UTR = 'TAGGTTGCTTCTTTTAGTGGTTTGCA'
_3UTR = 'TTTTCGTCTCTTATTATTAAACCTTTAAAAACGCTATCCTTGACTTTATCTGTACTTTGC'
wei = 'AATAAAAGCAGGCTCTGAGTGTTTAAATCTATTTTTCTTTCATTC'  # input('please enter wei:')
mei2 = 'GCTAGC'  # enzyme 2

cw = 0  # frequence
best = []  # optime
ps = []  # current result
val = 100
sum = 0
new = []


# read data from mima.csv,zimu.csv
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


# print results
def out(ls):
    t = sorted(ls, key=lambda ls: ls[0], reverse=False)
    for line in t:
        print(line[0], "%.2f" % line[1], end="|  ")
        for i in range(len(line[2]) - 1):
            if line[2][0] != 0:
                print(line[2][i][0], line[2][i][1], end="  ")

        print("")
        print("sum:%.2lf" % line[2][-1])


# plot
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


# letter to condon
def LetterToCodon(ls, al, neihanzi1, neihanzi2, end, n, mei1, suiji, _5UTR, _3UTR, wei, mei2):
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
    dis = []
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


def dealCodon(codon, mei1, suiji, _5UTR, end, neihanzi1, neihanzi2, _3UTR, wei, mei2):
    new1 = codon[len(mei1) + len(suiji) + len(_5UTR) + 3: -(len(_3UTR) + len(wei) + len(mei2))]
    # print('\nCondons：')
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

        elif new1[i:i + len(neihanzi2)] == neihanzi2:
            i = i + len(neihanzi2)
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


def checkEnzyme(new, mei1, mei2):
    line = new[len(mei1): -len(mei2)]
    print(line)
    # line = ''
    flag_mei1 = 0
    flag_mei2 = 0
    matchobj = re.match(r"(.*)GAATTC(.*)", line)
    if (matchobj == None):
        flag_mei1 = 0
    else:
        flag_mei2 = 1

    matchobj = re.match(r"(.*)GCTAGC(.*)", line)
    if (matchobj == None):
        flag_mei2 = 0
    else:
        flag_mei2 = 1
    return flag_mei1, flag_mei2


def checkIntron(new, mei1, mei2, neiahnzi1, neihanzi2):
    # line=" AAAGUAUGUCCCUACUAACCAGTTT"
    flag_long = 0
    # flag_short = 0
    line = new[len(mei1): -len(mei2)]
    matchobj = re.match(r"(.*)GUAUGU(.*)UACUAAC(.*)CAG(.*)GUAUGU(.*)UACUAAC(.*)CAG(.*)", new)
    if (matchobj == None):
        flag_long = 1
    else:
        temp_neihanzi1 = 'GUAUGU' + matchobj.group(2) + 'UACUAAC' + 'CAG'
        temp_neihanzi2 = 'GUAUGU' + matchobj.group(5) + 'UACUAAC' + 'CAG'
        if ((temp_neihanzi1 == neihanzi1 or temp_neihanzi1 == neihanzi2) and (
                temp_neihanzi2 == neihanzi2 or temp_neihanzi2 == neihanzi1)):
            flag_long = 0
            print('intron 1：')
            print(matchobj.group(2))
            print('intron 2：')
            print(matchobj.group(5))
        else:
            flag_long = 1


        if (flag_long == 0):
            temp_check = matchobj.group(1)
            temp_matchobj = re.match(r"(.*)GUAUGU(.*)", temp_check)
            if (temp_matchobj == None):
                flag_long = 0
            else:
                flag_long = 1
                print('A-A-B-C')

        if (flag_long == 0):
            temp_check = matchobj.group(3)
            temp_matchobj = re.match(r"(.*)CAG(.*)", temp_check)
            if (temp_matchobj == None):
                flag_long = 0
            else:
                flag_long = 1
                print('A-B-C-C')

        if (flag_long == 0):
            temp_check = matchobj.group(4)
            temp_matchobj = re.match(r"(.*)UACUAAC(.*)", temp_check)
            if (temp_matchobj == None):
                flag_long = 0
            else:
                flag_long = 1
                print('A-B-C-B-A-B-C')

        if (flag_long == 0):
            temp_check = matchobj.group(4)
            temp_matchobj = re.match(r"(.*)GUAUGU(.*)", temp_check)
            if (temp_matchobj == None):
                flag_long = 0
            else:
                flag_long = 1
                print('A-2A-B-C')

        if (flag_long == 0):
            print('There are no other introns.')
        else:
            print('There are other introns.')
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

        layout = QGridLayout()
        hlayout1 = QHBoxLayout()
        vlayout1 = QVBoxLayout()
        vlayout2 = QVBoxLayout()
        hlayout3 = QHBoxLayout()

        # Radio Button
        self.lettertocodon = QRadioButton("Letters to Codons")
        font = QFont('Times', 15)
        self.lettertocodon.setFont(font)
        self.lettertocodon.setChecked(True)
        self.lettertocodon.toggled.connect(lambda: self.btnstate(self.lettertocodon))

        self.codontoletter = QRadioButton("Codons to Letters")
        font = QFont('Times', 15)
        self.codontoletter.setFont(font)
        self.codontoletter.setChecked(False)
        self.codontoletter.toggled.connect(lambda: self.btnstate(self.codontoletter))

        self.labl1 = QLabel("Please input the letters:")
        font = QFont('Times', 15)
        self.labl1.setFont(font)

        self.editor1 = TextEdit()
        self.editor1.setAcceptRichText(False)
        self.editor1.setAutoFormatting(QTextEdit.AutoAll)
        font = QFont('Times', 12)
        self.editor1.setFont(font)
        self.editor1.setFontPointSize(12)


        self.editor2 = TextEdit()
        self.editor2.setAcceptRichText(False)
        self.editor2.setAutoFormatting(QTextEdit.AutoAll)
        font = QFont('Times', 12)
        self.editor2.setFont(font)
        self.editor2.setFontPointSize(12)


        '''
        self.button_translate = QPushButton('Translate')
        font = QFont('Times', 15)
        self.button_translate.setFont(font)
        self.button_translate.clicked.connect(self.onButtonClick)'''

        self.button_to = QPushButton("-->")
        font = QFont('Times', 20)
        font.setBold(True)
        self.button_to.setFont(font)
        self.button_to.setStatusTip("Click to translate")
        self.button_to.setToolTip("Click to translate")
        self.button_to.setStyleSheet("color: black;"
                                                            "background-color: #FF6347;"
                                                            "selection-color: red;"
                                                            "selection-background-color: blue;"
                                                            "border-width: 2px;"
                                                            "border-radius: 10px;"
                                                            "border-color: beige;"
                                                            "font: bold 14 px;"
                                                            "min-width: 5em;"
                                                            "padding: 6px;"
                                                            );
        self.button_to.clicked.connect(self.onButtonClick)
        '''
        self.codon = QtWidgets.QLabel()
        codons = QtGui.QPixmap('codons.png')
        self.codon.setPixmap(codons)

        self.letter = QtWidgets.QLabel()
        letters = QtGui.QPixmap('letters.png')
        self.codon.setPixmap(letters)
     '''

        self.path = None

        vlayout1.addWidget(self.lettertocodon)
        vlayout1.addWidget(self.codontoletter)
        #vlayout1.addWidget(self.codon)
        #hlayout1.addWidget(self.button_translate)
        vlayout2.addWidget(self.labl1)
        hlayout3.addWidget(self.editor1)
        hlayout3.addWidget(self.button_to)
        hlayout3.addWidget(self.editor2)

        vwg1 = QWidget()
        hwg1 = QWidget()
        vwg2 = QWidget()
        hwg3 = QWidget()

        vwg1.setLayout(vlayout1)
        #hwg1.setLayout(hlayout1)
        vwg2.setLayout(vlayout2)
        hwg3.setLayout(hlayout3)

        layout.addWidget(vwg1,0,0)
        layout.addWidget(hwg1,0,1)
        layout.addWidget(vwg2,1,0)
        layout.addWidget(hwg3,2,0)
        layout.setVerticalSpacing(0)

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

        self.flag = 1
        self.update_title()
        self.resize(1000,500)
        #self.setStyleSheet('QWidget{background-color:grey}' )
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("nefu.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.setWindowIcon(icon)
        self.show()


    def onButtonClick(self):
        sender = self.sender()
        print(sender.text() + ' has been ' + "pressed")
        if (self.flag == 1):
            al = self.editor1.toPlainText()
            al = al.upper()
            if (al == ""):
                QMessageBox.information(self, "Warning", "The entered data is empty.", QMessageBox.Ok)
                self.editor2.setText("")
            elif(re.search(r'\d', al)):
                QMessageBox.information(self, "Warning", "Please enter letters, you may enter some numbers.", QMessageBox.Ok)
                self.editor2.setText("")
            elif (len(al) > 50):
                QMessageBox.information(self, "Warning",
                                        "Sorry! The letters are too long, we can't translate them successfully!",
                                        QMessageBox.Ok)
                self.editor2.setText("")
            else:
                print(al)
                cnt = 0
                while (1):
                    cnt += 1
                    if(cnt==5):
                        break
                    new = []
                    new = LetterToCodon(LetterList, al, neihanzi1, neihanzi2, TerminationCodonList, n, mei1, suiji,
                                        _5UTR,
                                        _3UTR, wei, mei2)
                    print(new)
                    flag_mei1, flag_mei2 = checkEnzyme(new, mei1, mei2)
                    if (flag_mei1 == 1 or flag_mei2 == 1):
                        continue
                    flag_intron = checkIntron(new, mei1, mei2, neihanzi1, neihanzi2)
                    if (flag_intron == 1):
                        continue
                    break
                if(cnt==5):
                    text="666"
                    QMessageBox.information(self, "Warning", "Sorry! The letters are too long, we can't translate them successfully!", QMessageBox.Ok)
                    self.editor2.setText("")
                else:
                    text = new
                    self.editor2.setText(text)
        else:
            codon = self.editor1.toPlainText()
            if (codon == ""):
                QMessageBox.information(self, "Warning", "The entered data is empty.", QMessageBox.Ok)
                self.editor2.setText("")
            else:
                codon1 = dealCodon(codon, mei1, suiji, _5UTR, TerminationCodonList, neihanzi1, neihanzi2, _3UTR, wei,
                                   mei2)
                le1 = ''
                le1 = CodonToLetter(LetterList, codon1)
                print(le1)
                load()
                text = segment(le1)
                text2 = ''
                for word in text:
                    text2 = text2+word+' '
                text2 = text2.upper()
                self.editor2.setText(text2)

    def btnstate(self, btn):
        if btn.text() == "Letters to Codons":
            if btn.isChecked() == True:
                print(btn.text() + " is selected")
                self.labl1.setText('Please input the letters:')
                self.flag = 1
            else:
                print(btn.text() + " is deselected")

        if btn.text() == "Codons to Letters":
            if btn.isChecked() == True:
                print(btn.text() + " is selected")
                self.labl1.setText('Please input the condons:')
                self.flag = 2
            else:
                print(btn.text() + " is deselected")

    def block_signals(self, objects, b):
        for o in objects:
            o.blockSignals(b)

    def dialog_critical(self, s):
        dlg = QMessageBox(self)
        dlg.setText(s)
        dlg.setIcon(QMessageBox.Critical)
        dlg.show()

    def file_open(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open file", "",
                                              "HTML documents (*.html);Text documents (*.txt);All files (*.*)")

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
        path, _ = QFileDialog.getSaveFileName(self, "Save file", "",
                                              "HTML documents (*.html);Text documents (*.txt);All files (*.*)")

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
        self.editor1.setLineWrapMode(1 if self.editor1.lineWrapMode() == 0 else 0)


CodonList = []  # Codons arrays
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
# out(LetterList)
t = sorted(LetterList, key=lambda LetterList: LetterList[0], reverse=False)
# drawing(t)

getProbability(LetterList)
if __name__ == '__main__':
    app =QApplication(sys.argv)

    window = MainWindow()
    app.exec_()
