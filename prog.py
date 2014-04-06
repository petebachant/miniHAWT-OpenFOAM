#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 21:34:46 2013

@author: pete
"""
from PyQt4 import QtGui
import foampy

if __name__ == '__main__':
    import sys
    app = QtGui.QApplication(sys.path)
    pbarwin = foampy.ProgressBar()
    pbarwin.show()
    app.exec_()