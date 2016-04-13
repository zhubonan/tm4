# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 17:16:27 2016

@author: Bonan
"""

from PyQt4 import QtGui
import sys
import specImag

class ExampleApp(QtGui.QMainWindow, specImag.Ui_MainWindow):
    def __init__(self, parent=None):
        super(ExampleApp, self).__init__(parent)
        self.setupUi(self)
            
def main():
    app = QtGui.QApplication(sys.argv)
    form = ExampleApp()
    form.show()
    app.exec_()
    

main()