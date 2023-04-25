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
#include <iostream>

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
    for(int i = EXOMODE; i < numOfVar; i++){
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

//function that reads file
void MainFrame::openFile(){
    QString fileName = QFileDialog::getOpenFileName(this, "Open a file", "/Users/abner/");
    QFile data(fileName);
    data.open(QFile::ReadOnly);

    QTextStream in(&data);
    QString line[numOfVar];
    QString inp[numOfVar];
    QString outputs[numOfVar];
    lineCountOpen = 0;
    dropCountOpen = 0;

    for(int i = EXOMODE; i < NUMVAR; i++)
    {
        line[i]  = in.readLine();
        inp[i].append(line[i]);
        QStringList partString = inp[i].split(" = ");
        outputs[i] = partString.at(1);
        int index;
        switch(i){
        case EXOMODE: case DISPLAY: case MODEL: case METHOD: case MESHTYPE: case SAVEMESH: case DGTYPE: case BASISTYPE: case ENRICHTYPE:
        case TSMETHOD: case FREESURF: case BC: case SRCTYPE: case SAVESRC: case SNAPSHFMT: case DXUSEIS:
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
    QString nReps = lineInputsi[31]->text();
    numReps = nReps.toInt();
    //reading the blocks
    if(numReps > 0){
        //refresh();
        //reading the blocks
        for(int l = 1; l <= numReps; l++){
            int s = 1;
            for(int v = 1; v <= 9; v++){
                line[NBLOCKS + 9*l - 9 + s]  = in.readLine();
                inp[NBLOCKS + 9*l - 9 + s].append(line[NBLOCKS + 9*l - 9 + s]);
                QStringList partString = inp[NBLOCKS + 9*l - 9 + s].split(" = ");
                outputs[NBLOCKS + 9*l - 9 + s] = partString.at(1);
                int index;
                switch (v) {
                case 9:
                    index = outputs[NBLOCKS + 9*l - 9 + s].toInt();
                    dropInputsi[NBLOCKS + 9*l - 9 + s-lineCountOpen]->setCurrentIndex(index);
                    dropCountOpen++;
                    break;
                default:
                    lineInputsi[NBLOCKS + 9*l - 9 + s-dropCountOpen]->setText(outputs[NBLOCKS + 9*l - 9 + s]);
                    lineCountOpen++;
                    break;
                }
                s++;
            }
        }
    }else{
        return;
    }
}

//New function that actualizes the interface
void MainFrame::refresh(){
    QString nReps = lineInputsi[31]->text();
    numReps = nReps.toInt();
    if(numOfVar > NUMVAR){
        do{
            QWidget *temp = entriesi.last();
            delete(temp);
            entriesi.removeLast();
            layoutsi.removeLast();
            labeli.removeLast();
            lineInputsi.removeLast();
            lineCountLayout--;
            numOfVar--;
        }while(numOfVar > NUMVAR);
    }
//Creates the new labels and inputs entrances
    if(numReps > 0){
        for(int l = 1; l <= numReps; l++){
            int s = 1;
            for(int v = NUMVAR; v <= (NUMVAR + 8); v++){
                entriesi.insert(NBLOCKS + 9*l - 9 + s,new QWidget(scrollWidget));
                layoutsi.insert(NBLOCKS + 9*l - 9 + s,new QHBoxLayout(entriesi[NBLOCKS + 9*l - 9 + s]));
                switch(v){
                case (NUMVAR + 8):
                    dropInputsi.insert(NBLOCKS + 9*l - 9 + s-lineCountLayout, new QComboBox);
                    for(int j = 0; j < blocktype_1Size; j++){
                        dropInputsi[NBLOCKS + 9*l - 9 + s-lineCountLayout]->addItem(QString(blocktype_1Entries[j]));
                    }
                    labeli.insert(NBLOCKS + 9*l - 9 + s,new QLabel(listOfVar[v] + QString::number(l) + " = "));
                    layoutsi[NBLOCKS + 9*l - 9 + s]->addWidget(labeli[NBLOCKS + 9*l - 9 + s]);
                    layoutsi[NBLOCKS + 9*l - 9 + s]->addWidget(dropInputsi[NBLOCKS + 9*l - 9 + s-lineCountLayout]);
                    scrollLayout->addWidget(entriesi[NBLOCKS + 9*l - 9 + s]);
                    dropCountLayout++;
                    break;
                default:
                    labeli.insert(NBLOCKS + 9*l - 9 + s,new QLabel(listOfVar[v] + QString::number(l) + " = "));
                    lineInputsi.insert(NBLOCKS + 9*l - 9 + s -dropCountLayout, new QLineEdit);
                    layoutsi[NBLOCKS + 9*l - 9 + s]->addWidget(labeli[NBLOCKS + 9*l - 9 + s]);
                    layoutsi[NBLOCKS + 9*l - 9 + s]->addWidget(lineInputsi[NBLOCKS + 9*l - 9 + s -dropCountLayout]);
                    scrollLayout->addWidget(entriesi[NBLOCKS + 9*l - 9 + s]);
                    lineCountLayout++;
                    break;
                }
                numOfVar++;
                s++;
            }
        }
    }
}

void MainFrame::print(QFile *qfile){
    QString printX;
    lineCountPrint = 0;
    dropCountPrint = 0;
    for(int i = EXOMODE; i < NUMVAR; i++){
        int index;
        QString input;
        QMessageBox errorMessage;
        bool isInt;
        bool isFlt;
        bool isCord;
        bool isTensor;
        //Switch prints variables outside the block
        switch(i){
        case EXOMODE: case DISPLAY: case MODEL: case METHOD: case MESHTYPE: case SAVEMESH: case DGTYPE: case BASISTYPE: case ENRICHTYPE:
        case TSMETHOD: case FREESURF: case BC: case SRCTYPE: case SAVESRC: case SNAPSHFMT: case DXUSEIS:
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
            input = lineInputsi[i-dropCountPrint]->displayText();
            isFlt = checkFloat(input);
            if(isFlt){
                printX = printX + labeli[i]->text() + input + "\n";
                lineCountPrint++;
            }else{
                errorMessage.setText("The variable "+labeli[i]->text()+" only recives real numbers");
                errorMessage.exec();
                return;
            }
            break;
        case SWIDTH: case SRCLOC: case SRCVECTOR: case REC0: case RECD: case SEISMOi:
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
    //Printing the blocks of variables
    if(numReps > 0){
        for(int l = 1; l <= numReps; l++){
            int s = 1;
            for(int v = 1; v <= 9; v++){
                int index;
                QString input;
                bool isCord;
                bool isFlt;
                QMessageBox errorMessage;
                switch (v) {
                case 4:
                    input = lineInputsi[NBLOCKS + 9*l - 9 + s-dropCountPrint]->displayText();
                    isCord = checkXYZCordinates(input);
                    if(isCord){
                        printX = printX + labeli[NBLOCKS + 9*l - 9 + s]->text() + input + "\n";
                        lineCountPrint++;
                    }else{
                        errorMessage.setText("The variable "+labeli[NBLOCKS + 9*l - 9 + s]->text()+" only recives 3d coordinates (x,y,z)");
                        errorMessage.exec();
                        return;
                    }
                    break;
                case 9:
                    index = dropInputsi[NBLOCKS + 9*l - 9 + s-lineCountPrint]->currentIndex();
                    printX = printX + labeli[NBLOCKS + 9*l - 9 + s]->text() + QString::number(index) + "\n";
                    dropCountPrint++;
                    break;
                default:
                    input = lineInputsi[NBLOCKS + 9*l - 9 + s-dropCountPrint]->displayText();
                    isFlt = checkFloat(input);
                    if(isFlt){
                        printX = printX + labeli[NBLOCKS + 9*l - 9 + s]->text() + input + "\n";
                        lineCountPrint++;
                    }else{
                        errorMessage.setText("The variable "+labeli[NBLOCKS + 9*l - 9 + s]->text()+" only recives real numbers");
                        errorMessage.exec();
                        return;
                    }
                    break;
                }
                s++;
            }
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
        }   //std::ranges::all_of(s.begin(), s.end(),[](char c){ return isdigit(c) != 0; })
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

