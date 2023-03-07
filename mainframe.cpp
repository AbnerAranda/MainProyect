#include "mainframe.h"
#include "ui_mainframe.h"
#include <QLabel>
#include <QGridLayout>
#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QString>
#include <QFile>
#include <QFileDialog>
#include <QMessageBox>
#include <QTextStream>
#include <QMap>
#include <QDebug>
#include <iostream>
#include <QFile>
#include <QScrollArea>
#include <QList>
#include <QComboBox>
#include <string.h>

MainFrame::MainFrame(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::MainFrame)
{
    ui->setupUi(this);

    //number of variables
    numOfVar = sizeof(listOfVar)/sizeof(const char*);

    //main scroll area
    QScrollArea *scrollArea = new QScrollArea;
    scrollArea->setWidgetResizable(true);

    //widget that holds the scroll area
    scrollWidget = new QWidget;
    scrollLayout = new QVBoxLayout(scrollWidget);
    scrollArea->setWidget(scrollWidget);
    scrollWidget->setLayout(scrollLayout);
    scrollArea->setStyleSheet("QScrollBar:horizontal { width: 5px; }");

    //layouts inside the scroll area
    int lineCount = 0;
    int dropCount = 0;
    for(int i = 0; i < numOfVar; i++){
        entries[i] = new QWidget(scrollWidget);
        layouts[i] = new QHBoxLayout(entries[i]);
        label[i] = new QLabel(tr(listOfVar[i]));
        //if for both qline and qcombobox
        if(strcmp(listOfVar[i],"EXOMODE: ") == 0){
            dropInputs[i-lineCount] = new QComboBox;
            dropInputs[i-lineCount]->addItem("0");
            dropInputs[i-lineCount]->addItem("1");
            dropInputs[i-lineCount]->addItem("2");
            dropInputs[i-lineCount]->addItem("3");
            layouts[i]->addWidget(label[i]);
            layouts[i]->addWidget(dropInputs[i-lineCount]);
            scrollLayout->addWidget(entries[i]);
            dropCount++;
        }else if(strcmp(listOfVar[i],"DISPLAY: ") == 0){
            dropInputs[i-lineCount] = new QComboBox;
            dropInputs[i-lineCount]->addItem("0");
            dropInputs[i-lineCount]->addItem("1");
            dropInputs[i-lineCount]->addItem("2");
            layouts[i]->addWidget(label[i]);
            layouts[i]->addWidget(dropInputs[i-lineCount]);
            scrollLayout->addWidget(entries[i]);
            dropCount++;
        }else if(strcmp(listOfVar[i],"MODEL: ") == 0){
            dropInputs[i-lineCount] = new QComboBox;
            dropInputs[i-lineCount]->addItem("Acustic");
            dropInputs[i-lineCount]->addItem("Elastic");
            dropInputs[i-lineCount]->addItem("Frac");
            dropInputs[i-lineCount]->addItem("Elastoacust");
            dropInputs[i-lineCount]->addItem("Anisotropic");
            dropInputs[i-lineCount]->addItem("Aniso+Frac");
            layouts[i]->addWidget(label[i]);
            layouts[i]->addWidget(dropInputs[i-lineCount]);
            scrollLayout->addWidget(entries[i]);
            dropCount++;
        }else if(strcmp(listOfVar[i],"METHOD: ") == 0){
            dropInputs[i-lineCount] = new QComboBox;
            dropInputs[i-lineCount]->addItem("SGFD");
            dropInputs[i-lineCount]->addItem("SEM");
            dropInputs[i-lineCount]->addItem("DG");
            dropInputs[i-lineCount]->addItem("IGA");
            dropInputs[i-lineCount]->addItem("VSFD");
            dropInputs[i-lineCount]->addItem("EG");
            dropInputs[i-lineCount]->addItem("HEG");
            layouts[i]->addWidget(label[i]);
            layouts[i]->addWidget(dropInputs[i-lineCount]);
            scrollLayout->addWidget(entries[i]);
            dropCount++;
        }else if(strcmp(listOfVar[i],"MESHTYPE: ") == 0){
            dropInputs[i-lineCount] = new QComboBox;
            dropInputs[i-lineCount]->addItem("0");
            dropInputs[i-lineCount]->addItem("4");
            layouts[i]->addWidget(label[i]);
            layouts[i]->addWidget(dropInputs[i-lineCount]);
            scrollLayout->addWidget(entries[i]);
            dropCount++;
        }else if(strcmp(listOfVar[i],"SAVEMESH: ") == 0){
            dropInputs[i-lineCount] = new QComboBox;
            dropInputs[i-lineCount]->addItem("no");
            dropInputs[i-lineCount]->addItem("yes");
            layouts[i]->addWidget(label[i]);
            layouts[i]->addWidget(dropInputs[i-lineCount]);
            scrollLayout->addWidget(entries[i]);
            dropCount++;
        }else if(strcmp(listOfVar[i],"TSMETHOD: ") == 0){
            dropInputs[i-lineCount] = new QComboBox;
            dropInputs[i-lineCount]->addItem("FDM");
            dropInputs[i-lineCount]->addItem("RK4");
            dropInputs[i-lineCount]->addItem("LW4");
            layouts[i]->addWidget(label[i]);
            layouts[i]->addWidget(dropInputs[i-lineCount]);
            scrollLayout->addWidget(entries[i]);
            dropCount++;
        }else if(strcmp(listOfVar[i],"FREESURF: ") == 0){
            dropInputs[i-lineCount] = new QComboBox;
            dropInputs[i-lineCount]->addItem("no");
            dropInputs[i-lineCount]->addItem("yes");
            layouts[i]->addWidget(label[i]);
            layouts[i]->addWidget(dropInputs[i-lineCount]);
            scrollLayout->addWidget(entries[i]);
            dropCount++;
        }else if(strcmp(listOfVar[i],"BC: ") == 0){
            dropInputs[i-lineCount] = new QComboBox;
            dropInputs[i-lineCount]->addItem("Neumann");
            dropInputs[i-lineCount]->addItem("Taper");
            dropInputs[i-lineCount]->addItem("Periodic");
            layouts[i]->addWidget(label[i]);
            layouts[i]->addWidget(dropInputs[i-lineCount]);
            scrollLayout->addWidget(entries[i]);
            dropCount++;
        }else if(strcmp(listOfVar[i],"SRCTYPE: ") == 0){
            dropInputs[i-lineCount] = new QComboBox;
            dropInputs[i-lineCount]->addItem("pint src");
            dropInputs[i-lineCount]->addItem("compressional plane wave");
            dropInputs[i-lineCount]->addItem("shear plane wave");
            dropInputs[i-lineCount]->addItem("3");
            layouts[i]->addWidget(label[i]);
            layouts[i]->addWidget(dropInputs[i-lineCount]);
            scrollLayout->addWidget(entries[i]);
            dropCount++;
        }else if(strcmp(listOfVar[i],"SAVESRC: ") == 0){
            dropInputs[i-lineCount] = new QComboBox;
            dropInputs[i-lineCount]->addItem("no");
            dropInputs[i-lineCount]->addItem("yes");
            layouts[i]->addWidget(label[i]);
            layouts[i]->addWidget(dropInputs[i-lineCount]);
            scrollLayout->addWidget(entries[i]);
            dropCount++;
        }else if(strcmp(listOfVar[i],"SNAPSHFMT: ") == 0){
            dropInputs[i-lineCount] = new QComboBox;
            dropInputs[i-lineCount]->addItem("None");
            dropInputs[i-lineCount]->addItem("Binary");
            dropInputs[i-lineCount]->addItem("NetCDF");
            dropInputs[i-lineCount]->addItem("Slice");
            dropInputs[i-lineCount]->addItem("VTK");
            dropInputs[i-lineCount]->addItem("Exodus II ");
            layouts[i]->addWidget(label[i]);
            layouts[i]->addWidget(dropInputs[i-lineCount]);
            scrollLayout->addWidget(entries[i]);
            dropCount++;
        }else if(strcmp(listOfVar[i],"DXUSEIS: ") == 0){
            dropInputs[i-lineCount] = new QComboBox;
            dropInputs[i-lineCount]->addItem("displacement seismograms");
            dropInputs[i-lineCount]->addItem("X-derivative seismograms");
            layouts[i]->addWidget(label[i]);
            layouts[i]->addWidget(dropInputs[i-lineCount]);
            scrollLayout->addWidget(entries[i]);
            dropCount++;
        }else if(strcmp(listOfVar[i],"BLOCKTYPE_1: ") == 0){
            dropInputs[i-lineCount] = new QComboBox;
            dropInputs[i-lineCount]->addItem("layer");
            dropInputs[i-lineCount]->addItem("blk");
            layouts[i]->addWidget(label[i]);
            layouts[i]->addWidget(dropInputs[i-lineCount]);
            scrollLayout->addWidget(entries[i]);
            dropCount++;
        }else{
            lineInputs[i-dropCount] = new QLineEdit;
            layouts[i]->addWidget(label[i]);
            layouts[i]->addWidget(lineInputs[i-dropCount]);
            scrollLayout->addWidget(entries[i]);
            lineCount++;
        }
    }

    //buttons
    selectButton = new QPushButton(tr("Select File"));
    saveButton = new QPushButton(tr("Save"));
    newButton = new QPushButton(tr("New"));

    connect(selectButton, &QPushButton::clicked, this, &MainFrame::openFile);
    connect(saveButton, &QPushButton::clicked, this, &MainFrame::saveFile);
    connect(newButton, &QPushButton::clicked, this, &MainFrame::newFile);

    //Layout
    QHBoxLayout *selectLayout = new QHBoxLayout;
    selectLayout->addWidget(selectButton);

    QHBoxLayout *saveNewLayout = new QHBoxLayout;
    saveNewLayout->addWidget(saveButton);
    saveNewLayout->addWidget(newButton);

    //main layout
    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addLayout(selectLayout);
    mainLayout->addWidget(scrollArea);
    mainLayout->addLayout(saveNewLayout);

    setLayout(mainLayout);
    setWindowTitle("Input Window");
}

void MainFrame::saveFile(){
    QString fileName = QFileDialog::getOpenFileName(this, "Open a file", "/Users/abner/");
    QFile data(fileName);
    data.open(QFile::WriteOnly|QFile::Truncate);
    print(&data);
}

void MainFrame::newFile(){

    QString fileName = QFileDialog::getSaveFileName(this,tr("Save Text"), "",tr("Text (*.txt);;All Files (*)"));

    if(fileName.isEmpty())
                return;
            else{
                    QFile file(fileName);
                    if (!file.open(QIODevice::WriteOnly)) {
                        QMessageBox::information(this, tr("Unable to open file"),
                        file.errorString());
                        return;
                    }
                    print(&file);
                }
}

void MainFrame::openFile(){
    QString fileName = QFileDialog::getOpenFileName(this, "Open a file", "/Users/abner/");
    QFile data(fileName);
    data.open(QFile::ReadOnly);

    QTextStream in(&data);
    QString line[numOfVar];
    QString inp[numOfVar];
    QString outputs[numOfVar];
    int lineCount = 0;
    int dropCount = 0;
    for(int i=0; i<numOfVar; i++)
    {
        line[i]  = in.readLine();
        inp[i].append(line[i]);
        QStringList partString = inp[i].split(": ");
        outputs[i] = partString.at(1);
        if((strcmp(listOfVar[i],"EXOMODE: ") == 0)||(strcmp(listOfVar[i],"DISPLAY: ") == 0)||(strcmp(listOfVar[i],"MODEL: ") == 0)
                ||(strcmp(listOfVar[i],"METHOD: ") == 0)||(strcmp(listOfVar[i],"MESHTYPE: ") == 0)||(strcmp(listOfVar[i],"SAVEMESH: ") == 0)
                ||(strcmp(listOfVar[i],"TSMETHOD: ") == 0)||(strcmp(listOfVar[i],"FREESURF: ") == 0)||(strcmp(listOfVar[i],"BC: ") == 0)
                ||(strcmp(listOfVar[i],"SRCTYPE: ") == 0)||(strcmp(listOfVar[i],"SAVESRC: ") == 0)||(strcmp(listOfVar[i],"SNAPSHFMT: ") == 0)
                ||(strcmp(listOfVar[i],"DXUSEIS: ") == 0)||(strcmp(listOfVar[i],"BLOCKTYPE_1: ") == 0)){
            int index = outputs[i].toInt();
            //dropInputs[i-lineCount]->setCurrentText(outputs[i]);
            dropInputs[i-lineCount]->setCurrentIndex(index);
            dropCount++;
        }else{
            lineInputs[i-dropCount]->setText(outputs[i]);
            lineCount++;
        }
    }
}

void MainFrame::print(QFile *qfile){
    QString printX;
    int lineCount = 0;
    int dropCount = 0;
    for(int i=0; i<numOfVar; i++){
        if((strcmp(listOfVar[i],"EXOMODE: ") == 0)||(strcmp(listOfVar[i],"DISPLAY: ") == 0)||(strcmp(listOfVar[i],"MODEL: ") == 0)
                ||(strcmp(listOfVar[i],"METHOD: ") == 0)||(strcmp(listOfVar[i],"MESHTYPE: ") == 0)||(strcmp(listOfVar[i],"SAVEMESH: ") == 0)
                ||(strcmp(listOfVar[i],"TSMETHOD: ") == 0)||(strcmp(listOfVar[i],"FREESURF: ") == 0)||(strcmp(listOfVar[i],"BC: ") == 0)
                ||(strcmp(listOfVar[i],"SRCTYPE: ") == 0)||(strcmp(listOfVar[i],"SAVESRC: ") == 0)||(strcmp(listOfVar[i],"SNAPSHFMT: ") == 0)
                ||(strcmp(listOfVar[i],"DXUSEIS: ") == 0)||(strcmp(listOfVar[i],"BLOCKTYPE_1: ") == 0)){ 
            int index = dropInputs[i-lineCount]->currentIndex();
            printX = printX + label[i]->text() + QString::number(index) + "\n";
            dropCount++;
        }else{
            printX = printX + label[i]->text() + lineInputs[i-dropCount]->displayText() + "\n";
            lineCount++;
        }
    }
    QTextStream out(qfile);
    out << printX;
}

MainFrame::~MainFrame()
{
    delete ui;
}

