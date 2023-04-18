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

    //main scroll area
    QScrollArea *scrollArea = new QScrollArea(this);
    scrollArea->setWidgetResizable(true);

    //widget that holds the scroll area
    scrollWidget = new QWidget(this);
    scrollLayout = new QVBoxLayout(scrollWidget);
    scrollArea->setWidget(scrollWidget);
    scrollWidget->setLayout(scrollLayout);
    scrollArea->setStyleSheet("QScrollBar:horizontal { width: 5px; }");

    //layouts inside the scroll area
    for(int i = EXOMODE; i <= BLOCKTYPE_1; i++){
        entriesi.insert(i,new QWidget(scrollWidget));
        layoutsi.insert(i,new QHBoxLayout(entriesi[i]));
        labeli.insert(i,new QLabel(listOfVar[i] + " = "));
        switch(i){
        case EXOMODE:
            printDropItems(exomodeEntries, exomodeSize, i);
            break;
        case DISPLAY:
            printDropItems(displayEntries, displaySize, i);
            break;
        case MODEL:
            printDropItems(modelEntries, modelSize, i);
            break;
        case METHOD:
            printDropItems(methodEntries, methodSize, i);
            break;
        case MESHTYPE:
            printDropItems(meshtypeEntries, meshtypeSize, i);
            break;
        case SAVEMESH: case FREESURF: case SAVESRC:
            printDropItems(boolEntries, boolSize, i);
            break;
        case DGTYPE:
            printDropItems(dgtypeEntries, dgtypeSize, i);
            break;
        case BASISTYPE:
            printDropItems(basistypeEntries, basistypeSize, i);
            break;
        case ENRICHTYPE:
            printDropItems(basistypeEntries, basistypeSize, i);
            break;
        case TSMETHOD:
            printDropItems(tsmethodEntries, tsmethodSize, i);
            break;
        case BC:
            printDropItems(bcEntries, bcSize, i);
            break;
        case SRCTYPE:
            printDropItems(srctypeEntries, srctypeSize, i);
            break;
        case SNAPSHFMT:
            printDropItems(snapshfmtEntries, snapshfmtSize, i);
            break;
        case DXUSEIS:
            printDropItems(dxuseisEntries, dxuseisSize, i);
            break;
        case BLOCKTYPE_1:
            printDropItems(blocktype_1Entries, blocktype_1Size, i);
            break;
        default:
            lineInputsi.insert(i-dropCountLayout, new QLineEdit);
            layoutsi[i]->addWidget(labeli[i]);
            layoutsi[i]->addWidget(lineInputsi[i-dropCountLayout]);
            scrollLayout->addWidget(entriesi[i]);
            lineCountLayout++;
            break;
        }
    }

    //buttons
    selectButton = new QPushButton(tr("Select File"), this);
    refreshButton = new QPushButton(tr("Refresh"), this);
    saveButton = new QPushButton(tr("Save"), this);
    newButton = new QPushButton(tr("New"), this);

    connect(selectButton, &QPushButton::clicked, this, &MainFrame::openFile);
    connect(refreshButton, &QPushButton::clicked, this, &MainFrame::refresh);
    connect(saveButton, &QPushButton::clicked, this, &MainFrame::saveFile);
    connect(newButton, &QPushButton::clicked, this, &MainFrame::newFile);

    //Layout
    QHBoxLayout *selectLayout = new QHBoxLayout;
    selectLayout->addWidget(selectButton);
    selectLayout->addWidget(refreshButton);

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
    for(int i = EXOMODE; i <= BLOCKTYPE_1; i++)
    {
        line[i]  = in.readLine();
        inp[i].append(line[i]);
        QStringList partString = inp[i].split(" = ");
        outputs[i] = partString.at(1);
        int index;
        switch(i){
        case EXOMODE: case DISPLAY: case MODEL: case METHOD: case MESHTYPE: case SAVEMESH: case DGTYPE: case BASISTYPE: case ENRICHTYPE:
        case TSMETHOD: case FREESURF: case BC: case SRCTYPE: case SAVESRC: case SNAPSHFMT: case DXUSEIS: case BLOCKTYPE_1:
            index = outputs[i].toInt();
            dropInputsi[i-lineCountOpen]->setCurrentIndex(index);
            dropCountOpen++;
            break;
        default:
            lineInputsi[i-dropCountOpen]->setText(outputs[i]);
            lineCountOpen++;
            break;
        }
    }
}

void MainFrame::refresh(){
    QString nReps = lineInputsi[31]->text();
    int numReps = nReps.toInt();
    if(numOfVar > NUMVAR){
        do{
            entriesi.removeLast();
            layoutsi.removeLast();
            labeli.removeLast();
            lineInputsi.removeLast();
            lineCountLayout--;
            numOfVar--;
        }while(numOfVar > NUMVAR);
    }

    if(numReps > 0){
        for(int l = 0; l < numReps; l++){
            entriesi.insert(NUMVAR + l,new QWidget(scrollWidget));
            layoutsi.insert(NUMVAR + l,new QHBoxLayout(entriesi[NUMVAR + l]));
            labeli.insert(NUMVAR + l,new QLabel(listOfVar[X_1] + " = "));
            lineInputsi.insert(NUMVAR + l -dropCountLayout, new QLineEdit);
            layoutsi[NUMVAR + l]->addWidget(labeli[NUMVAR + l]);
            layoutsi[NUMVAR + l]->addWidget(lineInputsi[NUMVAR + l -dropCountLayout]);
            scrollLayout->addWidget(entriesi[NUMVAR + l]);
            lineCountLayout++;
            numOfVar++;
        }
    }
}

void MainFrame::print(QFile *qfile){
    QString printX;
    for(int i = EXOMODE; i <= BLOCKTYPE_1; i++){
        int index;
        QString input;
        QMessageBox errorMessage;
        bool isInt;
        bool isCord;
        bool isTensor;
        switch(i){
        case EXOMODE: case DISPLAY: case MODEL: case METHOD: case MESHTYPE: case SAVEMESH: case DGTYPE: case BASISTYPE: case ENRICHTYPE:
        case TSMETHOD: case FREESURF: case BC: case SRCTYPE: case SAVESRC: case SNAPSHFMT: case DXUSEIS: case BLOCKTYPE_1:
            index = dropInputsi[i-lineCountPrint]->currentIndex();
            printX = printX + labeli[i]->text() + QString::number(index) + "\n";
            dropCountPrint++;
            break;
        case LOGFILE: case OUTDIR:
            input = lineInputsi[i-dropCountPrint]->displayText();
            isInt = checkInteger(input);
            if(!isInt){
                printX = printX + labeli[i]->text() + input + "\n";
                lineCountPrint++;
            }else{
                errorMessage.setText("The variable "+labeli[i]->text()+" only recives string");
                errorMessage.exec();
                return;
            }
            break;
        case DTSAVE: case DGPENALTY: case TMAX: case CFL: case TAPALPHA: case TAPBETA: case TWIDTH: case PKFREQ: case SRCAMP:
        case X_1 ... Z_1: case VMAX_1 ... RHO_1://0.2
            input = lineInputsi[i-dropCountPrint]->displayText();
            isInt = checkInteger(input);
            if(!isInt){
                printX = printX + labeli[i]->text() + input + "\n";
                lineCountPrint++;
            }else{
                errorMessage.setText("The variable "+labeli[i]->text()+" only recives real numbers");
                errorMessage.exec();
                return;
            }
            break;
        case SWIDTH: case SRCLOC: case SRCVECTOR: case REC0: case RECD: case SEISMOi: case NELEM_1:
            input = lineInputsi[i-dropCountPrint]->displayText();
            isCord = checkXYZCordinates(input);
            if(isCord){
                printX = printX + labeli[i]->text() + input + "\n";
                lineCountPrint++;
            }else{
                errorMessage.setText("The variable "+labeli[i]->text()+" only recives 3d coordinates (x,y,z)");
                errorMessage.exec();
                return;
            }
            break;
        case SRCTENSOR:
            input = lineInputsi[i-dropCountPrint]->displayText();
            isTensor = checkTensor(input);
            if(isTensor){
                printX = printX + labeli[i]->text() + input + "\n";
                lineCountPrint++;
            }else{
                errorMessage.setText("The variable "+labeli[i]->text()+" only recives 6 components (xx,yy,zz,xy,xz,yz)");
                errorMessage.exec();
                return;
            }
            break;
        case NTSAVE: case DIMENSIONS: case FEMORD: case ENRICHORD: case TAPERLEN: case SRCFUNC: case SHOTS: case SNAPSHRES: case NREC:
        case NSEISMO: case MESHFAC: case NBLOCKS:
            input = lineInputsi[i-dropCountPrint]->displayText();
            isInt = checkInteger(input);
            if(isInt){
                printX = printX + labeli[i]->text() + input + "\n";
                lineCountPrint++;
            }else{
                errorMessage.setText("The variable "+labeli[i]->text()+" only recives integers");
                errorMessage.exec();
                return;
            }
            break;
        default:
            input = lineInputsi[i-dropCountPrint]->displayText();
            printX = printX + labeli[i]->text() + input + "\n";
            lineCountPrint++;
            break;
        }
    }
    QTextStream out(qfile);
    out << printX;
}

void MainFrame::printDropItems(QList <QString> entry, int arraySize, int index){
    dropInputsi.insert(index-lineCountLayout, new QComboBox);
    for(int j = 0; j < arraySize; j++){
        dropInputsi[index-lineCountLayout]->addItem(QString(entry[j]));
    }
    layoutsi[index]->addWidget(labeli[index]);
    layoutsi[index]->addWidget(dropInputsi[index-lineCountLayout]);
    scrollLayout->addWidget(entriesi[index]);
    dropCountLayout++;
}

bool MainFrame::checkInteger(QString input){
    int itr = 0;
    if(input.size() == 0){
        return false;
    }
    if(input[0] == '-'){
        itr = 1;
    }
    for(int i=itr; i<input.size(); i++){
        if(input[i]<'0' || input[i]>'9'){
            return false;
        }
    }
    return true;
}

bool MainFrame::checkPositive(QString input){
    bool isNeg = false;
    if(input.size() == 0){
        return false;
    }
    if(input[0] == '-'){
        isNeg = true;
    }
    return isNeg;
}

bool MainFrame::checkFloat(QString input){
    bool isFloat = false;
    if(input.contains(".")){
        isFloat = true;
    }
    return isFloat;
}

bool MainFrame::checkXYZCordinates(QString input){
    bool isXYZ = false;
    int space = 0;
    for(int i=0; i<input.size(); i++){
        if(input[i] == ' '){
            space++;
        }
    }

    if(space == 2){
        isXYZ = true;
    }
    return isXYZ;
}

bool MainFrame::checkTensor(QString input){
    bool isTensor = false;
    int space = 0;
    for(int i=0; i<input.size(); i++){
        if(input[i] == ' '){
            space++;
        }
    }
    if(space == 5){
        isTensor = true;
    }
    return isTensor;
}

MainFrame::~MainFrame()
{
    delete ui;
}

