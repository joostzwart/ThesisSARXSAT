import sys
sys.path.append('../..')
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import ReaderData
import identify
import threading
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import random
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

StyleSheet = '''
QCheckBox {
    spacing: 15px;
    font: 12pt Lato;
    font-weight: 500;
    color: #fff;
}

QPushButton#loadoutput {
        border-width:4px;
        border-style: solid;
        background-color: white;
        border-color: #008CBA;
        border-radius: 8px;
        padding: 10px 15px;
        font: 12pt Lato;
        font-weight: 800;
        }

QPushButton#loadoutput::hover {
        border-width: 3px;
        border-style: solid;
        background-color: #008CBA;
        border-color: #008CBA;
        border-radius: 8px;
        padding: 10px 20px;
        font: 12pt Lato;
        font-weight: 800;
        color: white;
        }

QPushButton#loadinput {
        border-width: 1px;
        border-style: solid;
        background-color: white;
        border-color: #008CBA;
        border-radius: 8px;
        padding: 10px 20px;
        font: 12pt Lato;
        font-weight: 500;
        }

QPushButton#loadinput::hover {
        border-width: 5px;
        border-style: solid;
        background-color: #008CBA;
        border-color: #008CBA;
        border-radius: 8px;
        padding: 10px 20px;
        font: 12pt Lato;
        font-weight: 500;
        color: white;
        }

QCheckBox::indicator {
    width:  11px;
    height: 11px;
    background-color: rgb(255, 255, 255);
    border-style: solid;
    border-width: 5px;
    border-color: rgb(255, 255, 255);
    border-radius: 10.4px;
}

QCheckBox::indicator:checked {
    background-color: rgb(0, 136, 247);
}


QLabel {
    font: 12pt Lato;
    color: #fff;
}

QGroupBox {
    font: 15pt Lato;
    color: #fff;
    font-weight: 600; 
}

QRadioButton {
    spacing: 15px;
    font: 12pt Lato;
    color: #fff;   
}

QRadioButton::indicator {
    width:  11px;
    height: 11px;
    background-color: rgb(255, 255, 255);
    border-style: solid;
    border-width: 5px;
    border-color: rgb(255, 255, 255);
    border-radius: 10.4px;
}

QRadioButton::indicator:checked {
    background-color: rgb(0, 136, 247);
}
'''


class App(QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = 'SARX identification'
        self.left = 0
        self.top = 0
        self.boxes = {}
        self.labels = {}
        self.initUI()
        self.showMaximized()
        #self.pal = self.palette()
        #self.pal.setColor(self.backgroundRole(), QColor(205, 255, 255))
        #self.setPalette(self.pal)
        Image = QImage("pythonfiles/GUI/background.jpeg")
        #sImage = oImage.scaled(QSize(self.primaryScreen().size().width(), self.primaryScreen().size().height()))
        palette = QPalette()
        palette.setBrush(10, QBrush(Image))
        self.setPalette(palette)

    def initUI(self):
        self.centralwidget = QWidget()
        self.setCentralWidget(self.centralwidget)
        self.setWindowTitle(self.title)
        self.show()

        ## Create plot for dataset
        self.fig, self.ax = plt.subplots(figsize=(15,15))
        self.plotOutput = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.plotOutput, self)
        self.ax.set_title('Dataset')
        self.ax.set_xlabel('[k]')

        ## Create plot for switching sequence
        self.fig_switching, self.ax_switching = plt.subplots(figsize=(15,15))
        self.plot_switching = FigureCanvas(self.fig_switching)
        self.toolbar_switching = NavigationToolbar(self.plot_switching,self)
        self.ax_switching.set_title('Switching Sequence')

        ## Create informationscreen for logging
        self.informationScreen=QFrame()
        self.screenLabels=QVBoxLayout()
        self.informationScreen.setLayout(self.screenLabels)
        self.logLabel1=QLabel("---")
        self.logLabel2 = QLabel("---")
        self.logLabel3 = QLabel("---")
        self.screenLabels.addWidget(self.logLabel1)
        self.screenLabels.addWidget(self.logLabel2)
        self.screenLabels.addWidget(self.logLabel3)
        self.pal = self.informationScreen.palette()
        self.pal.setColor(self.informationScreen.backgroundRole(), QColor(40, 40,40))
        self.informationScreen.setPalette(self.pal)
        self.informationScreen.setAutoFillBackground(True)
        self.logLabel1.setWordWrap(True)
        self.logLabel2.setWordWrap(True)
        self.logLabel3.setWordWrap(True)

        ## Create two panels
        self.panels = QHBoxLayout(self.centralwidget)
        self.leftPanel=QGroupBox("Parameters")
        self.panels.addWidget(self.leftPanel,1)
        self.rightPanel = QGroupBox("Plots and information")
        self.panels.addWidget(self.rightPanel,2)

        ## Layout of right panel
        self.groupsRight = QVBoxLayout()
        self.rightPanel.setLayout(self.groupsRight)
        self.groupsRight.addWidget(self.plotOutput)
        self.groupsRight.addWidget(self.toolbar)
        self.groupsRight.addWidget(self.plot_switching)
        self.groupsRight.addWidget(self.toolbar_switching)
        self.groupsRight.addWidget(self.informationScreen)
        self.bottomRight = QFrame()
        self.groupsRight.addStretch(1)
        self.groupsRight.addWidget(self.bottomRight)

        ## Bottom right
        self.groupBottomRight=QHBoxLayout()
        self.bottomRight.setLayout(self.groupBottomRight)

        # Create a button in the window
        self.button = QPushButton('Start', self)
        self.button.setStyleSheet("background-color: green")
        self.button.clicked.connect(self.on_click)
        self.groupBottomRight.addWidget(self.button)

        # Create a progressbar
        self.progressbar=QProgressBar()
        self.progressbar.setMaximum(100)
        self.groupBottomRight.addWidget(self.progressbar)

        ## Layout of left panel
        self.modelsLeft=QVBoxLayout()
        self.leftPanel.setLayout(self.modelsLeft)
        self.panelgeneral=QFrame()
        self.panelsarx=QFrame()
        self.panelsarx.setFrameShape(QFrame.StyledPanel| QFrame.Raised)
        self.panelpwarx=QFrame()
        self.panelpwarx.setFrameShape(QFrame.StyledPanel)
        self.bottomLeft = QFrame()
        self.panelsimulate = QFrame()
        self.pal = self.panelpwarx.palette()
        self.pal.setColor(self.panelpwarx.backgroundRole(), QColor(40, 40,40, 150))
        self.panelpwarx.setPalette(self.pal)
        self.panelsimulate.setPalette(self.pal)
        self.modelsLeft.addWidget(self.panelgeneral, 2)
        self.modelsLeft.addWidget(self.panelsarx, 1)
        self.modelsLeft.addWidget(self.panelpwarx, 1)
        self.modelsLeft.addWidget(self.panelsimulate,3)
        self.modelsLeft.addWidget(self.bottomLeft,2)
        self.panelpwarx.setAutoFillBackground(True)
        self.pal = self.bottomLeft.palette()
        self.pal.setColor(self.bottomLeft.backgroundRole(), QColor(40, 40,40))
        self.bottomLeft.setPalette(self.pal)
        self.bottomLeft.setAutoFillBackground(True)
        self.panelsimulate.setAutoFillBackground(True)
        self.pal = self.panelgeneral.palette()
        self.pal.setColor(self.panelgeneral.backgroundRole(), QColor(70, 96, 220, 100))
        self.panelgeneral.setPalette(self.pal)
        self.panelgeneral.setAutoFillBackground(True)
        self.panelsarx.setPalette(self.pal)
        self.panelsarx.setAutoFillBackground(True)

        ## Button Groups
        self.modeltype=QButtonGroup()
        self.ownDataSet=QButtonGroup()

        ## Layout of general
        self.splitgeneral=QGridLayout()
        self.panelgeneral.setLayout(self.splitgeneral)
        self.splitgeneral.setColumnStretch(0,5)
        self.splitgeneral.setColumnStretch(1, 2)

        #ny
        self.ny = QSpinBox()
        self.ny.setMinimum(1)
        self.ny.setMaximum(100)
        self.splitgeneral.addWidget(self.ny,0,1)
        self.labelny=QLabel("Order of the output:")
        self.splitgeneral.addWidget(self.labelny, 0, 0)

        #nu
        self.nu = QSpinBox()
        self.nu.setMinimum(1)
        self.nu.setMaximum(100)
        self.splitgeneral.addWidget(self.nu,1,1)
        self.labelnu = QLabel("Order of the input:")
        self.splitgeneral.addWidget(self.labelnu, 1, 0)

        #delta
        self.labelDelta = QLabel("Bound on the error: ")
        self.splitgeneral.addWidget(self.labelDelta, 2, 0)
        self.delta = QDoubleSpinBox()
        self.delta.setDecimals(4)
        self.delta.setMinimum(0.001)
        self.delta.setLocale(QLocale('English'))
        self.delta.setSingleStep(0.001)
        self.splitgeneral.addWidget(self.delta, 2, 1)

        #split large datasets
        self.checkSplitDataset = QCheckBox("Split Large Datasets")
        self.splitgeneral.addWidget(self.checkSplitDataset,3,0)
        self.checkSplitDataset.stateChanged.connect(self.on_check_split_dataset)

        #blocks
        self.labelBoxes = QLabel("Number of blocks: ")
        self.splitgeneral.addWidget(self.labelBoxes,4,0)
        self.blocks = QSpinBox()
        self.blocks.setMinimum(1)
        self.blocks.setMaximum(100)
        self.splitgeneral.addWidget(self.blocks,4,1)

        ## Layout of pwarx
        self.splitpwarx=QGridLayout()
        self.panelpwarx.setLayout(self.splitpwarx)
        self.splitpwarx.setColumnStretch(0,5)
        self.splitpwarx.setColumnStretch(1, 2)
        self.checkPwarx = QRadioButton("PWARX")
        self.splitpwarx.addWidget(self.checkPwarx, 0, 0)
        self.modeltype.addButton(self.checkPwarx)
        self.checkPwarx.toggled.connect(self.on_check_pwarx)

        #switching limit
        self.Ls = QSpinBox()
        self.Ls.setMinimum(1)
        self.Ls.setMaximum(100)
        self.splitpwarx.addWidget(self.Ls, 2, 1)
        self.labelLs = QLabel("Maximum number of switches:")
        self.splitpwarx.addWidget(self.labelLs, 2, 0)

        ## Layout of simulate
        self.splitsimulate=QGridLayout()
        self.panelsimulate.setLayout(self.splitsimulate)

        #Labels
        self.checkGenerateModels = QRadioButton("Generate random system")
        self.splitsimulate.addWidget(self.checkGenerateModels,0,0)
        self.ownDataSet.addButton(self.checkGenerateModels)
        self.checkGenerateModels.setChecked(False)
        self.labelLengthDataset = QLabel("Length of the dataset: ")
        self.splitsimulate.addWidget(self.labelLengthDataset,1,0)
        self.labelNumberOfModels = QLabel("Number of models: ")
        self.splitsimulate.addWidget(self.labelNumberOfModels,2,0)
        self.labelNoise = QLabel("Bound on the noise: ")
        self.splitsimulate.addWidget(self.labelNoise,3,0)
        self.labelNumberOfInputs = QLabel("Number of inputs: ")
        self.splitsimulate.addWidget(self.labelNumberOfInputs,4,0)
        self.labelNumberOfOutputs = QLabel("Number of outputs: ")
        self.splitsimulate.addWidget(self.labelNumberOfOutputs,5,0)
        self.labelSeed = QLabel("Seed for randomness. 0 for a random seed")
        self.splitsimulate.addWidget(self.labelSeed,6,0)

        # length of dataset
        self.LengthDataset = QSpinBox()
        self.LengthDataset.setMinimum(10)
        self.LengthDataset.setMaximum(10000)
        self.splitsimulate.addWidget(self.LengthDataset,1,1)

        # number of models
        self.numberOfModels = QSpinBox()
        self.numberOfModels.setRange(1, 20)
        self.splitsimulate.addWidget(self.numberOfModels,2,1)

        # noise
        self.noise = QDoubleSpinBox()
        self.noise.setDecimals(4)
        self.noise.setMinimum(0.001)
        self.noise.setLocale(QLocale('English'))
        self.noise.setSingleStep(0.001)
        self.splitsimulate.addWidget(self.noise,3,1)

        # number of inputs
        self.numberOfInputs = QSpinBox()
        self.numberOfInputs.setRange(1, 10)
        self.splitsimulate.addWidget(self.numberOfInputs,4,1)

        # number of outputs
        self.numberOfOutputs = QSpinBox()
        self.numberOfOutputs.setRange(1, 10)
        self.splitsimulate.addWidget(self.numberOfOutputs,5,1)

        # seed
        self.seed = QSpinBox()
        self.seed.setRange(0, 10000)
        self.splitsimulate.addWidget(self.seed,6,1)
        self.checkGenerateModels.toggled.connect(self.on_check_generate_models)

        # hide objects from start
        self.blocks.hide()
        self.numberOfModels.hide()
        self.labelBoxes.hide()
        self.labelNumberOfModels.hide()
        self.noise.hide()
        self.labelNoise.hide()
        self.numberOfInputs.hide()
        self.labelNumberOfInputs.hide()
        self.numberOfOutputs.hide()
        self.labelNumberOfOutputs.hide()
        self.labelLengthDataset.hide()
        self.LengthDataset.hide()
        self.labelSeed.hide()
        self.seed.hide()
        self.Ls.hide()
        self.labelLs.hide()

        ## SARX
        self.splitsarx = QGridLayout()
        self.panelsarx.setLayout(self.splitsarx)
        self.splitsarx.setColumnStretch(0,5)
        self.splitsarx.setColumnStretch(1, 2)


        # dw
        self.dw = QSpinBox()
        self.dw.setMinimum(1)
        self.dw.setMaximum(100)
        self.splitsarx.addWidget(self.dw, 2, 1)
        self.labeldw = QLabel("Minimum dwell time:")
        self.splitsarx.addWidget(self.labeldw, 2, 0)

        # radiobutton
        self.radioSarx=QRadioButton("SARX")
        self.radioSarx.toggled.connect(self.on_check_sarx)
        self.splitsarx.addWidget(self.radioSarx,0,0)
        self.modeltype.addButton(self.radioSarx)
        self.radioSarx.setChecked(True)

        ## Bottom Left
        self.gridLayout=QGridLayout()
        self.gridLayout.setColumnStretch(0,5)
        self.gridLayout.setColumnStretch(1, 2)
        self.gridLayout.setColumnStretch(2, 5)
        self.bottomLeft.setLayout(self.gridLayout)
        self.radioLoadDataset=QRadioButton("Use own dataset")
        self.loadButtonInput = QPushButton('Load Input',objectName="loadoutput")
        self.gridLayout.addWidget(self.radioLoadDataset,0,0)
        self.ownDataSet.addButton(self.radioLoadDataset)
        self.radioLoadDataset.toggled.connect(self.on_check_own_dataset)
        self.gridLayout.addWidget(self.loadButtonInput,1,0)
        self.loadButtonInput.clicked.connect(self.on_click_load_input)

        self.loadButtonOutput = QPushButton('Load Output',objectName="loadoutput")
        self.gridLayout.addWidget(self.loadButtonOutput,1,2,)
        self.loadButtonOutput.clicked.connect(self.on_click_load_output)
        self.loadButtonOutput.hide()
        self.loadButtonInput.hide()



    @pyqtSlot()
    def on_click(self):
        Ls = int(self.Ls.text())
        dw = int(self.dw.text())
        nu = int(self.nu.text())
        ny = int(self.ny.text())
        delta=float(self.delta.text())
        split=int(self.checkSplitDataset.isChecked())
        blocks=int(self.blocks.text())
        inputs=int(self.numberOfInputs.text())
        outputs=int(self.numberOfOutputs.text())
        nt=float(self.noise.text())
        models=int(self.numberOfModels.text())
        T=int(self.LengthDataset.text())
        if self.checkGenerateModels.isChecked():
            seed=int(self.seed.text())
        else:
            seed=0
        pwarx=int(self.checkPwarx.isChecked())
        modelgeneration=2+int(self.checkGenerateModels.isChecked())
        #QMessageBox.question(self, 'Parameters', "Number of inputs: " + NOI + ", dwelltime = " + dw, QMessageBox.Ok,
        #                     QMessageBox.Ok)
        print(modelgeneration)
        if self.checkGenerateModels.isChecked():
            identify.main(self,True,nu=nu,ny=ny,dw=dw,seed=seed,pwarx=pwarx,delta=delta,split=split,blocks=blocks,inputs=inputs,outputs=outputs,
                      nt=nt,models=models,T=T,modelgeneration=modelgeneration,Ls=Ls)
        else:
            identify.main(self,True,nu=nu,ny=ny,dw=dw,seed=seed,pwarx=pwarx,delta=delta,split=split,blocks=blocks,modelgeneration=modelgeneration,inputfile=self.input_file,outputfile=self.output_file,Ls=Ls)



    @pyqtSlot()
    def on_check_split_dataset(self):
        if self.checkSplitDataset.isChecked():
            self.blocks.show()
            self.labelBoxes.show()
        else:
            self.blocks.hide()
            self.labelBoxes.hide()

    @pyqtSlot()
    def on_check_sarx(self):
        self.pal = self.panelsarx.palette()
        if self.radioSarx.isChecked():
            self.dw.show()
            self.labeldw.show()
            self.pal.setColor(self.panelsarx.backgroundRole(), QColor(70, 96, 220, 100))
            self.panelsarx.setPalette(self.pal)
        else:
            self.labeldw.hide()
            self.dw.hide()
            self.pal.setColor(self.panelsarx.backgroundRole(),  QColor(40, 40, 40, 150))
            self.panelsarx.setPalette(self.pal)

    @pyqtSlot()
    def on_check_pwarx(self):
        self.pal = self.panelpwarx.palette()
        if self.checkPwarx.isChecked():
            self.pal.setColor(self.panelpwarx.backgroundRole(), QColor(70, 96, 220, 100))
            self.panelpwarx.setPalette(self.pal)
            self.Ls.show()
            self.labelLs.show()
        else:
            self.pal.setColor(self.panelpwarx.backgroundRole(),  QColor(40, 40, 40, 150))
            self.panelpwarx.setPalette(self.pal)
            self.Ls.hide()
            self.labelLs.hide()

    @pyqtSlot()
    def on_check_generate_models(self):
        if self.checkGenerateModels.isChecked():
            self.numberOfModels.show()
            self.labelNumberOfModels.show()
            self.noise.show()
            self.labelNoise.show()
            self.numberOfInputs.show()
            self.labelNumberOfInputs.show()
            self.numberOfOutputs.show()
            self.labelNumberOfOutputs.show()
            self.labelLengthDataset.show()
            self.LengthDataset.show()
            self.labelSeed.show()
            self.seed.show()
        else:
            self.numberOfModels.hide()
            self.labelNumberOfModels.hide()
            self.noise.hide()
            self.labelNoise.hide()
            self.numberOfInputs.hide()
            self.labelNumberOfInputs.hide()
            self.numberOfOutputs.hide()
            self.labelNumberOfOutputs.hide()
            self.labelLengthDataset.hide()
            self.LengthDataset.hide()
            self.labelSeed.hide()
            self.seed.hide()

    @pyqtSlot()
    def on_check_own_dataset(self):
        if self.radioLoadDataset.isChecked():
            self.loadButtonInput.show()
            self.loadButtonOutput.show()
        else:
            self.loadButtonInput.hide()
            self.loadButtonOutput.hide()
            self.ax.clear()
            self.plotOutput.draw()


    @pyqtSlot()
    def on_click_load_input(self):
        self.input_file = QFileDialog.getOpenFileName(self, 'Open file', 'Data')[0]
        if hasattr(self, 'output_file') and hasattr(self, 'input_file'):
            (self.input, self.output, self.inputs, self.outputs, self.T, self.r)=ReaderData.Read(3, 3, 2, input=self.input_file, output=self.output_file)
            for i in range(len(self.input)):
                if i==0:
                    self.ax.plot(self.input[i], 'r-', label="input")
                else:
                    self.ax.plot(self.input[i], 'r-')
            for i in range(len(self.output)):
                if i==0:
                    self.ax.plot(self.output[i], 'r-', label="output")
                else:
                    self.ax.plot(self.output[i], 'r-')
            self.ax.legend(loc=2)
            self.plotOutput.draw()


    @pyqtSlot()
    def on_click_load_output(self):
        self.output_file = QFileDialog.getOpenFileName(self, 'Open file', 'Data')[0]
        if hasattr(self, 'output_file') and hasattr(self, 'input_file'):
            (self.input, self.output, self.inputs, self.outputs, self.T, self.r) = ReaderData.Read(3, 3, 2,
                        input=self.input_file,output=self.output_file)
            for i in range(len(self.input)):
                if i==0:
                    self.ax.plot(self.input[i], 'r-', label="input")
                else:
                    self.ax.plot(self.input[i], 'r-')
            for i in range(len(self.output)):
                if i==0:
                    self.ax.plot(self.output[i], 'b-', label="output")
                else:
                    self.ax.plot(self.output[i], 'b-')
            self.ax.legend(loc=2)
            self.plotOutput.draw()

    def openFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
        return fileName

if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyleSheet(StyleSheet)
    ex = App()
    sys.exit(app.exec_())