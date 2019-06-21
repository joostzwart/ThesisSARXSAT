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
    spacing: 5px;
    font: 14pt Lato;
}

QCheckBox::indicator {
    width:  20px;
    height: 20px;
}

QLabel {
    font: 14pt Lato;
}
'''


class App(QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = 'SARX identification'
        self.left = 0
        self.top = 0
        self.width = 1920
        self.height = 1080
        self.boxes = {}
        self.labels = {}
        self.initUI()

    def initUI(self):
        self.centralwidget = QWidget()
        self.setCentralWidget(self.centralwidget)
        self.setWindowTitle(self.title)
        self.resize(self.width,self.height)
        self.show()

        # Create boxes
        self.create_box("dw", 300, "Dwell time: ")
        self.create_box("nu", 400, "Order of the input: ")
        self.create_box("ny", 500, "Order of the output: ")

        # Create plot
        self.fig, self.ax = plt.subplots()
        self.plotOutput = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.plotOutput, self)
        ## Create two panels
        self.panels = QHBoxLayout(self.centralwidget)
        self.leftPanel=QGroupBox("leftPanel")
        self.panels.addWidget(self.leftPanel,1)
        self.rightPanel = QGroupBox("rightPanel")
        self.panels.addWidget(self.rightPanel,2)

        ## Layout of right panel
        self.groupsRight = QVBoxLayout()
        self.rightPanel.setLayout(self.groupsRight)
        self.groupsRight.addWidget(self.plotOutput)
        self.groupsRight.addWidget(self.toolbar)
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
        self.splitLeft=QVBoxLayout()
        self.leftPanel.setLayout(self.splitLeft)
        self.topLeft=QFrame()
        self.bottomLeft=QGroupBox()
        self.splitLeft.addWidget(self.topLeft,2)
        self.splitLeft.addWidget(self.bottomLeft,1)

        ## Top left
        self.columnsLeft = QHBoxLayout()
        self.topLeft.setLayout(self.columnsLeft)
        self.col1=QFrame()
        self.col2=QFrame()
        self.columnsLeft.addWidget(self.col1, 2)
        self.columnsLeft.addWidget(self.col2, 1)

        ## Bottom Left
        self.gridLayout=QGridLayout()
        self.gridLayout.setAlignment(Qt.AlignCenter)
        self.bottomLeft.setLayout(self.gridLayout)
        self.loadButtonInput = QPushButton('Load Input')
        self.gridLayout.addWidget(self.loadButtonInput,0,0,1,1)
        self.loadButtonInput.clicked.connect(self.on_click_load_input)
        self.loadButtonInput.setStyleSheet("background-color: blue")

        self.loadButtonOutput = QPushButton('Load Output')
        self.gridLayout.addWidget(self.loadButtonOutput,0,1,1,1)
        self.loadButtonOutput.clicked.connect(self.on_click_load_output)
        self.loadButtonOutput.setStyleSheet("background-color: blue")

        ## Column 1
        self.slotsCol1 = QVBoxLayout()
        self.col1.setLayout(self.slotsCol1)
        for k in self.labels:
            self.slotsCol1.addWidget(self.labels[k])
        self.labelSeed = QLabel("Seed for randomness. 0 for a random seed")
        self.slotsCol1.addWidget(self.labelSeed)
        self.checkSplitDataset = QCheckBox("Split Large Datasets")
        self.slotsCol1.addWidget(self.checkSplitDataset)
        self.labelBoxes=QLabel("Number of blocks: ")
        self.slotsCol1.addWidget(self.labelBoxes)
        self.labelDelta=QLabel("Bound on the error: ")
        self.slotsCol1.addWidget(self.labelDelta)
        self.checkGenerateModels=QCheckBox("Generate random system")
        self.slotsCol1.addWidget(self.checkGenerateModels)
        self.labelLengthDataset=QLabel("Length of the dataset: ")
        self.slotsCol1.addWidget(self.labelLengthDataset)
        self.labelNumberOfModels = QLabel("Number of models: ")
        self.slotsCol1.addWidget(self.labelNumberOfModels)
        self.labelNoise = QLabel("Bound on the noise: ")
        self.slotsCol1.addWidget(self.labelNoise)
        self.labelNumberOfInputs = QLabel("Number of inputs: ")
        self.slotsCol1.addWidget(self.labelNumberOfInputs)
        self.labelNumberOfOutputs = QLabel("Number of outputs: ")
        self.slotsCol1.addWidget(self.labelNumberOfOutputs)
        self.checkPwarx=QCheckBox("State depended switching sequence (PWARX)")
        self.slotsCol1.addWidget(self.checkPwarx)
        self.slotsCol1.addStretch(1)

        ## Column 2
        self.slotsCol2 = QVBoxLayout()
        self.col2.setLayout(self.slotsCol2)
        for k in self.boxes:
            self.slotsCol2.addWidget(self.boxes[k])

        self.seed=QSpinBox()
        self.slotsCol2.addWidget(self.seed)
        self.spacer = QSpacerItem(200, 25, QSizePolicy.Expanding)
        self.slotsCol2.addItem(self.spacer)

        #blocks
        self.blocks = QSpinBox()
        self.blocks.setMinimum(1)
        self.blocks.setMaximum(100)
        self.slotsCol2.addWidget(self.blocks)
        #delta
        self.delta=QDoubleSpinBox()
        self.delta.setDecimals(4)
        self.delta.setMinimum (0.001)
        self.delta.setLocale(QLocale('English'))
        self.delta.setSingleStep(0.001)
        self.slotsCol2.addWidget(self.delta)
        self.slotsCol2.addItem(self.spacer)
        #length of dataset
        self.LengthDataset = QSpinBox()
        self.LengthDataset.setMinimum(10)
        self.LengthDataset.setMaximum(10000)
        self.slotsCol2.addWidget(self.LengthDataset)
        #number of models
        self.numberOfModels = QSpinBox()
        self.numberOfModels.setRange(1,20)
        self.slotsCol2.addWidget(self.numberOfModels)
        #noise
        self.noise = QDoubleSpinBox()
        self.noise.setDecimals(4)
        self.noise.setMinimum(0.001)
        self.noise.setLocale(QLocale('English'))
        self.noise.setSingleStep(0.001)
        self.slotsCol2.addWidget(self.noise)
        #number of inputs
        self.numberOfInputs = QSpinBox()
        self.numberOfInputs.setRange(1,10)
        self.slotsCol2.addWidget(self.numberOfInputs)
        #number of outputs
        self.numberOfOutputs = QSpinBox()
        self.numberOfOutputs.setRange(1,10)
        self.slotsCol2.addWidget(self.numberOfOutputs)
        #misc
        self.slotsCol2.addStretch(1)
        self.checkSplitDataset.stateChanged.connect(self.on_check_split_dataset)
        self.checkGenerateModels.stateChanged.connect(self.on_check_generate_models)

        ## hide objects from start
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


    @pyqtSlot()
    def on_click(self):
        dw = int(self.boxes["dw"].text())
        nu = int(self.boxes["nu"].text())
        ny = int(self.boxes["ny"].text())
        delta=float(self.delta.text())
        split=int(self.checkSplitDataset.isChecked())
        blocks=int(self.blocks.text())
        inputs=int(self.numberOfInputs.text())
        outputs=int(self.numberOfOutputs.text())
        nt=float(self.noise.text())
        models=int(self.numberOfModels.text())
        T=int(self.LengthDataset.text())
        seed=int(self.seed.text())
        pwarx=int(self.checkPwarx.isChecked())
        modelgeneration=2+int(self.checkGenerateModels.isChecked())
        #QMessageBox.question(self, 'Parameters', "Number of inputs: " + NOI + ", dwelltime = " + dw, QMessageBox.Ok,
        #                     QMessageBox.Ok)
        if self.checkGenerateModels.isChecked():
            identify.main(self,nu=nu,ny=ny,dw=dw,seed=seed,pwarx=pwarx,delta=delta,split=split,blocks=blocks,inputs=inputs,outputs=outputs,
                      nt=nt,models=models,T=T,modelgeneration=modelgeneration)
        else:
            identify.main(self,nu=nu,ny=ny,dw=dw,seed=seed,pwarx=pwarx,delta=delta,split=split,blocks=blocks,modelgeneration=modelgeneration)



    @pyqtSlot()
    def on_check_split_dataset(self):
        if self.checkSplitDataset.isChecked():
            self.blocks.show()
            self.labelBoxes.show()
        else:
            self.blocks.hide()
            self.labelBoxes.hide()

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


    @pyqtSlot()
    def on_click_load_input(self):
        self.input_file=self.openFileNameDialog()
        if hasattr(self, 'output_file') and hasattr(self, 'input_file'):
            if not hasattr(self, 'input'):
                (self.input, self.output, self.inputs, self.outputs, self.T, self.r)=ReaderData.Read(3, 3, 2, input=self.input_file, output=self.output_file)
                for i in range(len(self.input)):
                    self.ax.plot(self.input[i], 'r-')
                self.ax.plot(self.output[0], 'r-')
                self.plotOutput.draw()



    @pyqtSlot()
    def on_click_load_output(self):
        self.output_file = QFileDialog.getOpenFileName(self, 'Open file', 'Data')[0]
        if hasattr(self, 'output_file') and hasattr(self, 'input_file'):
            if not hasattr(self, 'input'):
                (self.input, self.output, self.inputs, self.outputs, self.T, self.r) = ReaderData.Read(3, 3, 2,
                                                                                                       input=self.input_file,
                                                                                                       output=self.output_file)
                for i in range(len(self.input)):
                    self.ax.plot(self.input[i], 'r-')
                self.ax.plot(self.output[0], 'r-')
                self.plotOutput.draw()

    def create_box(self,name,y,text):
        self.boxes[name] = QSpinBox(self)
        self.boxes[name].setMinimum(1)
        self.create_label(name,text)

    def create_label(self,name,text):
        self.labels[name] = QLabel(self)
        self.labels[name].setStyleSheet("font: 14pt Lato")
        self.labels[name].setText(text)

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