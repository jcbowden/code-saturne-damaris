# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# Standard modules
# -------------------------------------------------------------------------------

import logging

# -------------------------------------------------------------------------------
# Third-party modules
# -------------------------------------------------------------------------------

from code_saturne.Base.QtCore import *
from code_saturne.Base.QtGui import *
from code_saturne.Base.QtWidgets import *

# -------------------------------------------------------------------------------
# Application modules import
# -------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.Pages.VolumicConditionsForm import Ui_VolumicConditionsForm

from code_saturne.model.LocalizationModel import LocalizationModel
from code_saturne.model.GroundwaterModel import GroundwaterModel

# -------------------------------------------------------------------------------
# Widgets import
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# log config
# -------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("VolumicConditionsView")
log.setLevel(GuiParam.DEBUG)


# -------------------------------------------------------------------------------
# Main class
# -------------------------------------------------------------------------------

class VolumicConditionsView(QWidget, Ui_VolumicConditionsForm):
    """ Display available volumic treatments for a given zone """

    def __init__(self, parent, case, zone_name):
        QWidget.__init__(self, parent)
        Ui_VolumicConditionsForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.parent = parent
        self.zone_name = zone_name
        self.case.undoStopGlobal()

        for zone in LocalizationModel('VolumicZone', self.case).getZones():
            if zone.getLabel() == zone_name:
                self.zone = zone

        if GroundwaterModel(self.case).getGroundwaterModel() != "groundwater":
            # self.tabWidget.removeTab(4)
            pass  # TODO
        else:
            pass  # TODO

        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            self.saturneInitializationWidget.hide()
            self.neptuneInitializationWidget.setup(self.case, self.zone_name)
        else:
            self.neptuneInitializationWidget.hide()
            self.saturneInitializationWidget.setup(self.case, self.zone_name)
        self.porosityPage.setup(self.case, self.zone_name)
        self.headLossesPage.setup(self.case, self.zone_name)
        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            self.saturneSourceTermsWidget.hide()
            self.neptuneSourceTermsWidget.setup(self.case, self.zone_name)
        else:
            self.neptuneSourceTermsWidget.hide()
            self.saturneSourceTermsWidget.setup(self.case, self.zone_name)
        self.groundwaterLawPage.setup(self.case, self.zone_name)

        if not (self.zone.isNatureActivated("groundwater_law")):
            self.tabWidget.removeTab(4)
        if not (self.zone.isNatureActivated("source_term")):
            self.tabWidget.removeTab(3)
        if not (self.zone.isNatureActivated("head_losses")):
            self.tabWidget.removeTab(2)
        if not (self.zone.isNatureActivated("porosity")):
            self.tabWidget.removeTab(1)
        if not (self.zone.isNatureActivated("initialization")):
            self.tabWidget.removeTab(0)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
