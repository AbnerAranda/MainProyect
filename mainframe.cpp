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
        //switch statement for qline and qcombobox
        switch(i){
        case 0:
            printDropItems(exomodeEntries, exomodeSize, i, lineCount);
            dropCount++;
            break;
        case 1:
            printDropItems(displayEntries, displaySize, i, lineCount);
            dropCount++;
            break;
        case 6:
            printDropItems(modelEntries, modelSize, i, lineCount);
            dropCount++;
            break;
        case 8:
            printDropItems(methodEntries, methodSize, i, lineCount);
            dropCount++;
            break;
        case 10:
            printDropItems(meshtypeEntries, meshtypeSize, i, lineCount);
            dropCount++;
            break;
        case 11: case 22: case 28:
            printDropItems(boolEntries, boolSize, i, lineCount);
            dropCount++;
            break;
        case 19:
            printDropItems(tsmethodEntries, tsmethodSize, i, lineCount);
            dropCount++;
            break;
        case 23:
            printDropItems(bcEntries, bcSize, i, lineCount);
            dropCount++;
            break;
        case 27:
            printDropItems(srctypeEntries, srctypeSize, i, lineCount);
            dropCount++;
            break;
        case 39:
            printDropItems(snapshfmtEntries, snapshfmtSize, i, lineCount);
            dropCount++;
            break;
        case 40:
            printDropItems(dxuseisEntries, dxuseisSize, i, lineCount);
            dropCount++;
            break;
        case 56:
            printDropItems(blocktype_1Entries, blocktype_1Size, i, lineCount);
            dropCount++;
            break;
        default:
            lineInputs[i-dropCount] = new QLineEdit;
            layouts[i]->addWidget(label[i]);
            layouts[i]->addWidget(lineInputs[i-dropCount]);
            scrollLayout->addWidget(entries[i]);
            lineCount++;
            break;
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
        int index;
        switch(i){
        case 0: case 1: case 6: case 8: case 10: case 11: case 19:
        case 22: case 23: case 27:case 28:case 39:case 40: case 56:
            index = outputs[i].toInt();
            dropInputs[i-lineCount]->setCurrentIndex(index);
            dropCount++;
            break;
        default:
            lineInputs[i-dropCount]->setText(outputs[i]);
            lineCount++;
            break;
        }
    }
}

void MainFrame::print(QFile *qfile){
    QString printX;
    int lineCount = 0;
    int dropCount = 0;

    for(int i=0; i<numOfVar; i++){
        int index;
        QString input;
        QMessageBox errorMessage;
        bool isInt;
        switch(i){
        case 0: case 1: case 6: case 8: case 10: case 11: case 19:
        case 22: case 23: case 27:case 28:case 39:case 40: case 56:
            index = dropInputs[i-lineCount]->currentIndex();
            printX = printX + label[i]->text() + QString::number(index) + "\n";
            dropCount++;
            break;
        case 4: case 7: case 9: case 13: case 17: case 20: case 21:
        case 24: case 25: case 26: case 29 ... 34: case 38: case 41:
        case 44: case 46: case 47: case 52 ... 55:
            input = lineInputs[i-dropCount]->displayText();
            isInt = checkInteger(input);
            if(isInt){
                printX = printX + label[i]->text() + input + "\n";
                lineCount++;
            }else{
                errorMessage.setText("The variable "+label[i]->text()+" only recives integers");
                errorMessage.exec();
                return;
            }
            break;
        default:
            input = lineInputs[i-dropCount]->displayText();
            printX = printX + label[i]->text() + input + "\n";
            lineCount++;
            break;
        }
    }
    QTextStream out(qfile);
    out << printX;
}

void MainFrame::printDropItems(const char ** entry, int arraySize, int index, int lnCount){
    dropInputs[index-lnCount] = new QComboBox;
    for(int j = 0; j < arraySize; j++){
        dropInputs[index-lnCount]->addItem(QString(entry[j]));
    }
    layouts[index]->addWidget(label[index]);
    layouts[index]->addWidget(dropInputs[index-lnCount]);
    scrollLayout->addWidget(entries[index]);
}

bool MainFrame::checkInteger(QString input){
    //bool isNeg = false;
    int itr = 0;
    if(input.size() == 0){
        return false;
    }
    if(input[0] == '-'){
        //isNeg = true;
        itr = 1;
    }
    for(int i=itr; i<input.size(); i++){
        if(input[i]<'0' || input[i]>'9'){
            return false;
        }
    }
    return true;
}

MainFrame::~MainFrame()
{
    delete ui;
}

