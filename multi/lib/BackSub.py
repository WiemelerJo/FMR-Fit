# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'BackSub.ui'
#
# Created by: PyQt5 UI code generator 5.15.1
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_BackgroundSubtraction(object):
    def setupUi(self, BackgroundSubtraction):
        BackgroundSubtraction.setObjectName("BackgroundSubtraction")
        BackgroundSubtraction.resize(1063, 720)
        self.gridLayout_2 = QtWidgets.QGridLayout(BackgroundSubtraction)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_2.setSpacing(0)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.Button_load_main = QtWidgets.QPushButton(BackgroundSubtraction)
        self.Button_load_main.setEnabled(False)
        self.Button_load_main.setObjectName("Button_load_main")
        self.horizontalLayout_2.addWidget(self.Button_load_main)
        self.Button_load_secondary = QtWidgets.QPushButton(BackgroundSubtraction)
        self.Button_load_secondary.setObjectName("Button_load_secondary")
        self.horizontalLayout_2.addWidget(self.Button_load_secondary)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.View_Orig = Single_Plot(BackgroundSubtraction)
        self.View_Orig.setObjectName("View_Orig")
        self.verticalLayout.addWidget(self.View_Orig)
        self.View_Mod = Single_Plot(BackgroundSubtraction)
        self.View_Mod.setObjectName("View_Mod")
        self.verticalLayout.addWidget(self.View_Mod)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.View_Colour = Colour_Plot(BackgroundSubtraction)
        self.View_Colour.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.View_Colour.sizePolicy().hasHeightForWidth())
        self.View_Colour.setSizePolicy(sizePolicy)
        self.View_Colour.setObjectName("View_Colour")
        self.verticalLayout_3.addWidget(self.View_Colour)
        self.horizontalLayout.addLayout(self.verticalLayout_3)
        self.gridLayout.addLayout(self.horizontalLayout, 0, 0, 1, 1)
        self.gridLayout_2.addLayout(self.gridLayout, 0, 0, 1, 1)

        self.retranslateUi(BackgroundSubtraction)
        QtCore.QMetaObject.connectSlotsByName(BackgroundSubtraction)

    def retranslateUi(self, BackgroundSubtraction):
        _translate = QtCore.QCoreApplication.translate
        BackgroundSubtraction.setWindowTitle(_translate("BackgroundSubtraction", "Subtraction Tool"))
        self.Button_load_main.setText(_translate("BackgroundSubtraction", "Load main data"))
        self.Button_load_secondary.setText(_translate("BackgroundSubtraction", "Load secondary data"))
from lib.CustomWidgets import Colour_Plot, Single_Plot
