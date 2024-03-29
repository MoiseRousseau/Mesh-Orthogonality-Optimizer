#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# Author : Moise Rousseau (2019), email at moise.rousseau@polymtl.ca

import qtsalome
from salome.gui import helper
import salome.smesh.smeshBuilder as smb

class UI:

  def __init__(self,context):
  
    self.selectMesh = False #seeds mesh
    self.mesh = None #seeds mesh
    self.sg = context.sg

    Dialog = qtsalome.QDialog()
    Dialog.setObjectName("Dialog")
    Dialog.resize(600, 180)
    Dialog.setSizeGripEnabled(False)
    
    #principal layout
    self.gridLayout_main = qtsalome.QGridLayout(Dialog)
    
    #mesh
    self.label_mesh = qtsalome.QLabel(Dialog)
    self.pb_origMeshFile = qtsalome.QPushButton(Dialog)
    self.pb_origMeshFile.setCheckable(True)
    self.le_origMeshFile = qtsalome.QLineEdit(Dialog)
    self.le_origMeshFile.setReadOnly(True)
    self.gridLayout_main.addWidget(self.label_mesh, 0, 0)
    self.gridLayout_main.addWidget(self.pb_origMeshFile, 0, 1)
    self.gridLayout_main.addWidget(self.le_origMeshFile, 0, 2)
    
    #function type
    self.label_function = qtsalome.QLabel(Dialog)
    self.gridLayout_main.addWidget(self.label_function , 1, 0)
    self.menu_function = qtsalome.QComboBox(Dialog)
    self.menu_function.addItems(["Power", "Inverse", "Logarithm", "Exponential"])
    self.gridLayout_main.addWidget(self.menu_function , 1, 2)
    
    #penalization
    self.label_penalization = qtsalome.QLabel(Dialog)
    self.le_penalization = qtsalome.QLineEdit(Dialog)
    self.gridLayout_main.addWidget(self.label_penalization , 2, 0)
    self.gridLayout_main.addWidget(self.le_penalization, 2, 2)
    #max it
    self.label_it = qtsalome.QLabel(Dialog)
    self.label_it.setObjectName("label_it")
    self.le_it = qtsalome.QLineEdit(Dialog)
    self.gridLayout_main.addWidget(self.label_it , 3, 0)
    self.gridLayout_main.addWidget(self.le_it, 3, 2)
    
    #ok and cancel button
    self.splitter = qtsalome.QSplitter(Dialog)
    self.splitter.setOrientation(qtsalome.Qt.Horizontal)
    self.pb_okCancel = qtsalome.QDialogButtonBox(self.splitter)
    self.pb_okCancel.setOrientation(qtsalome.Qt.Horizontal)
    self.pb_okCancel.setStandardButtons(qtsalome.QDialogButtonBox.Cancel|qtsalome.QDialogButtonBox.Ok)
    self.gridLayout_main.addWidget(self.splitter, 4, 0)
    self.pb_okCancel.accepted.connect(Dialog.accept)
    self.pb_okCancel.rejected.connect(Dialog.reject)
    
    self.translateUi(Dialog)
    self.pb_origMeshFile.clicked.connect(self.setMeshInput)
    
    self.Dialog = Dialog
    return
    
  def translateUi(self, Dialog):
    Dialog.setWindowTitle("OrthOpt")
    
    self.pb_origMeshFile.setText("Select")
    self.label_mesh.setText("Select mesh:")
    self.label_function.setText("Penalization function type:")
    self.label_penalization.setText("Penalization power:")
    self.le_penalization.setText("3")
    self.label_it.setText("Max iteration:")
    self.le_it.setText("30")
    return
  
  def select(self):
      #sg.getObjectBrowser().selectionChanged.disconnect(self.select)
      objId = helper.getSObjectSelected()[0].GetObject()
      if self.selectMesh:
        self._selectMeshInput(objId)
      elif self.selectSurface:
        self._selectSurfaceInput( objId)
      return
      
      
  def _selectMeshInput(self, objId):
    self.mesh = objId
    if isinstance(self.mesh,smb.meshProxy):
      name = smb.GetName(self.mesh)
    elif isinstance(self.mesh,SMESH._objref_SMESH_Group):
      name = smb.smesh.smeshBuilder.GetName(self.mesh)
    elif isinstance(self.mesh,smb.submeshProxy):
      name = smb.smesh.smeshBuilder.GetName(self.mesh)
    else:
      self.mesh = None
      name = ""
    self.le_origMeshFile.setText(name)
    return
  
  def setMeshInput(self):
    if self.selectMesh == True:
      self.selectMesh = False
      self.sg.getObjectBrowser().selectionChanged.disconnect(self.select)
    else:
      self.selectMesh = True
      self.sg.getObjectBrowser().selectionChanged.connect(self.select)
      self.select()
    return
  
  def show(self):
    self.Dialog.show()
    res = self.Dialog.exec_()
    return res
    
    
