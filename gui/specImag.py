# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'specImag.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(979, 672)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.centralwidget)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.mplwiget = MatplotlibWidget(self.centralwidget)
        self.mplwiget.setMinimumSize(QtCore.QSize(700, 500))
        self.mplwiget.setObjectName(_fromUtf8("mplwiget"))
        self.horizontalLayout.addWidget(self.mplwiget)
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.label_4 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_4.setFont(font)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.verticalLayout_3.addWidget(self.label_4)
        self.listView = QtGui.QListView(self.centralwidget)
        self.listView.setMaximumSize(QtCore.QSize(16777215, 300))
        self.listView.setObjectName(_fromUtf8("listView"))
        self.verticalLayout_3.addWidget(self.listView)
        self.label_3 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_3.setFont(font)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.verticalLayout_3.addWidget(self.label_3)
        self.dataSelectBox = QtGui.QComboBox(self.centralwidget)
        self.dataSelectBox.setObjectName(_fromUtf8("dataSelectBox"))
        self.dataSelectBox.addItem(_fromUtf8(""))
        self.verticalLayout_3.addWidget(self.dataSelectBox)
        self.label_5 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_5.setFont(font)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.verticalLayout_3.addWidget(self.label_5)
        self.cropMinBox = QtGui.QSpinBox(self.centralwidget)
        self.cropMinBox.setObjectName(_fromUtf8("cropMinBox"))
        self.verticalLayout_3.addWidget(self.cropMinBox)
        self.label_6 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_6.setFont(font)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.verticalLayout_3.addWidget(self.label_6)
        self.cropMaxBox = QtGui.QSpinBox(self.centralwidget)
        self.cropMaxBox.setObjectName(_fromUtf8("cropMaxBox"))
        self.verticalLayout_3.addWidget(self.cropMaxBox)
        self.label_8 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_8.setFont(font)
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.verticalLayout_3.addWidget(self.label_8)
        self.scanWidthBox = QtGui.QSpinBox(self.centralwidget)
        self.scanWidthBox.setObjectName(_fromUtf8("scanWidthBox"))
        self.verticalLayout_3.addWidget(self.scanWidthBox)
        self.label_9 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_9.setFont(font)
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.verticalLayout_3.addWidget(self.label_9)
        self.scanHeightBox = QtGui.QSpinBox(self.centralwidget)
        self.scanHeightBox.setObjectName(_fromUtf8("scanHeightBox"))
        self.verticalLayout_3.addWidget(self.scanHeightBox)
        self.horizontalLayout.addLayout(self.verticalLayout_3)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.label_7 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_7.setFont(font)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.verticalLayout_2.addWidget(self.label_7)
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_2.addItem(spacerItem)
        self.label = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label.setFont(font)
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout_2.addWidget(self.label)
        self.filterWinBox = QtGui.QSpinBox(self.centralwidget)
        self.filterWinBox.setObjectName(_fromUtf8("filterWinBox"))
        self.verticalLayout_2.addWidget(self.filterWinBox)
        spacerItem1 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_2.addItem(spacerItem1)
        self.label_2 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_2.setFont(font)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout_2.addWidget(self.label_2)
        self.filterOrderBox = QtGui.QComboBox(self.centralwidget)
        self.filterOrderBox.setObjectName(_fromUtf8("filterOrderBox"))
        self.filterOrderBox.addItem(_fromUtf8(""))
        self.filterOrderBox.addItem(_fromUtf8(""))
        self.verticalLayout_2.addWidget(self.filterOrderBox)
        spacerItem2 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_2.addItem(spacerItem2)
        self.filterTestButton = QtGui.QPushButton(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.filterTestButton.setFont(font)
        self.filterTestButton.setObjectName(_fromUtf8("filterTestButton"))
        self.verticalLayout_2.addWidget(self.filterTestButton)
        self.horizontalLayout.addLayout(self.verticalLayout_2)
        self.horizontalLayout_2.addLayout(self.horizontalLayout)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 979, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.actionOpen = QtGui.QAction(MainWindow)
        self.actionOpen.setObjectName(_fromUtf8("actionOpen"))
        self.actionExit = QtGui.QAction(MainWindow)
        self.actionExit.setObjectName(_fromUtf8("actionExit"))
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addAction(self.actionExit)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "SpecImaging", None))
        self.label_4.setText(_translate("MainWindow", "List of avaliable Data", None))
        self.label_3.setText(_translate("MainWindow", "Selected Data", None))
        self.dataSelectBox.setItemText(0, _translate("MainWindow", "test", None))
        self.label_5.setText(_translate("MainWindow", "Crop Minimal", None))
        self.label_6.setText(_translate("MainWindow", "Crop Maximum", None))
        self.label_8.setText(_translate("MainWindow", "Width", None))
        self.label_9.setText(_translate("MainWindow", "Height", None))
        self.label_7.setText(_translate("MainWindow", "Filter Settings", None))
        self.label.setText(_translate("MainWindow", "Filter Window", None))
        self.label_2.setText(_translate("MainWindow", "Filter Order", None))
        self.filterOrderBox.setItemText(0, _translate("MainWindow", "2", None))
        self.filterOrderBox.setItemText(1, _translate("MainWindow", "3", None))
        self.filterTestButton.setText(_translate("MainWindow", "Test", None))
        self.menuFile.setTitle(_translate("MainWindow", "File", None))
        self.actionOpen.setText(_translate("MainWindow", "Open", None))
        self.actionExit.setText(_translate("MainWindow", "Exit", None))

from matplotlibwidget import MatplotlibWidget
