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

MainFrame::MainFrame(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::MainFrame)
{
    ui->setupUi(this);

    //Label list
    labelList[0] = new QLabel(tr("EXOMODE: "));
    labelList[1] = new QLabel(tr("DISPLAY: "));
    labelList[2] = new QLabel(tr("LOGFILE: "));
    labelList[3] = new QLabel(tr("OUTDIR: "));
    labelList[4] = new QLabel(tr("NTSAVE: "));
    labelList[5] = new QLabel(tr("DTSAVE: "));
    labelList[6] = new QLabel(tr("MODEL: "));
    labelList[7] = new QLabel(tr("DIMENSIONS: "));
    labelList[8] = new QLabel(tr("METHOD: "));
    labelList[9] = new QLabel(tr("FEMORD: "));
    labelList[10] = new QLabel(tr("MESHTYPE: "));
    labelList[11] = new QLabel(tr("SAVEMESH: "));
    labelList[12] = new QLabel(tr("DGTYPE: "));
    labelList[13] = new QLabel(tr("DGPENALTY: "));
    labelList[14] = new QLabel(tr("DGLUMPPING: "));
    numOfVar = sizeof(labelList)/sizeof(QLabel*);

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
    for(int i = 0; i < numOfVar; i++){
        entries[i] = new QWidget(scrollWidget);
        layouts[i] = new QHBoxLayout(entries[i]);
        label[i] = labelList[i];
        inputs[i] = new QLineEdit;
        layouts[i]->addWidget(label[i]);
        layouts[i]->addWidget(inputs[i]);
        scrollLayout->addWidget(entries[i]);
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
    for(int i=0; i<numOfVar; i++)
    {
        line[i]  = in.readLine();
        inp[i].append(line[i]);
        QStringList partString = inp[i].split(": ");
        outputs[i] = partString.at(1);
        inputs[i]->setText(outputs[i]);
    }
}

void MainFrame::print(QFile *qfile){
    QString printX;
    for(int i=0; i<numOfVar; i++){
        printX = printX + label[i]->text() + inputs[i]->displayText() + "\n";
    }
    QTextStream out(qfile);
    out << printX;
}

MainFrame::~MainFrame()
{
    delete ui;
}

