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

    //labels and inputs
    QLabel *exomodeLabel = new QLabel(tr("EXOMODE: "));
    exomodeInput = new QLineEdit;
    QLabel *displayLabel = new QLabel(tr("DISPLAY: "));
    displayInput = new QLineEdit;
    QLabel *logfileLabel = new QLabel(tr("LOGFILE: "));
    logfileInput = new QLineEdit;
    /*QLabel *outdirLabel = new QLabel(tr("OUTDIR: "));
    outdirInput = new QLineEdit;
    QLabel *ntsaveLabel = new QLabel(tr("NTSAVE: "));
    ntsaveInput = new QLineEdit;
    QLabel *dtsaveLabel = new QLabel(tr("DTSAVE: "));
    dtsaveInput = new QLineEdit;
    QLabel *modelLabel = new QLabel(tr("MODEL: "));
    modelInput = new QLineEdit;*/

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

    QScrollArea *scroll = new QScrollArea;

    QVBoxLayout *scrollLayout = new QVBoxLayout;

    QHBoxLayout *exomodeLayout = new QHBoxLayout;
    exomodeLayout->addWidget(exomodeLabel);
    exomodeLayout->addWidget(exomodeInput);

    QHBoxLayout *displayLayout = new QHBoxLayout;
    displayLayout->addWidget(displayLabel);
    displayLayout->addWidget(displayInput);

    QHBoxLayout *logfileLayout = new QHBoxLayout;
    displayLayout->addWidget(logfileLabel);
    displayLayout->addWidget(logfileInput);

    /*QHBoxLayout *outdirLayout = new QHBoxLayout;
    displayLayout->addWidget(outdirLabel);
    displayLayout->addWidget(outdirInput);

    QHBoxLayout *ntsaveLayout = new QHBoxLayout;
    displayLayout->addWidget(ntsaveLabel);
    displayLayout->addWidget(ntsaveInput);

    QHBoxLayout *dtsaveLayout = new QHBoxLayout;
    displayLayout->addWidget(dtsaveLabel);
    displayLayout->addWidget(dtsaveInput);

    QHBoxLayout *modelLayout = new QHBoxLayout;
    displayLayout->addWidget(modelLabel);
    displayLayout->addWidget(modelInput);*/

    scrollLayout->addLayout(exomodeLayout);
    scrollLayout->addLayout(displayLayout);
    scrollLayout->addLayout(logfileLayout);
    //scrollLayout->addLayout(outdirLayout);
    //scrollLayout->addLayout(ntsaveLayout);
    //scrollLayout->addLayout(dtsaveLayout);
    //scrollLayout->addLayout(modelLayout);


    QWidget *widget = new QWidget;
    widget->setLayout(scrollLayout);
    scroll->setWidget(widget);
    scroll->setStyleSheet("QScrollBar:horizontal { width: 5px; }");
    //scroll->setStyleSheet("QScrollBar:vertical { width: 5px; }");

    QHBoxLayout *saveNewLayout = new QHBoxLayout;
    saveNewLayout->addWidget(saveButton);
    saveNewLayout->addWidget(newButton);

    //main layout
    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addLayout(selectLayout);
    mainLayout->addWidget(scroll);
    mainLayout->addLayout(saveNewLayout);

    setLayout(mainLayout);
    setWindowTitle("Input Window");
}

void MainFrame::saveFile(){
    QString fileName = QFileDialog::getOpenFileName(this, "Open a file", "/Users/abner/");
    QFile data(fileName);
    data.open(QFile::WriteOnly|QFile::Truncate);

    QString printEXOMODE;
    QString printDISPLAY;

    printEXOMODE = "EXOMODE: " + exomodeInput->displayText() + "\n";
    printDISPLAY = "DISPLAY: " + displayInput->displayText() + "\n";
    QString printText = printEXOMODE + printDISPLAY;

    QTextStream out(&data);
    out << printText;
}

void MainFrame::newFile(){
    QString printEXOMODE;
    QString printDISPLAY;

    printEXOMODE = "EXOMODE: " + exomodeInput->displayText() + "\n";
    printDISPLAY = "DISPLAY: " + displayInput->displayText() + "\n";
    QString printText = printEXOMODE + printDISPLAY;

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
                    QTextStream out(&file);
                    out << printText;
                }
}

void MainFrame::openFile(){
    QString fileName = QFileDialog::getOpenFileName(this, "Open a file", "/Users/abner/");
    QFile data(fileName);
    data.open(QFile::ReadOnly);

    QTextStream in(&data);
    QString line[2];
    QString inputs[2];
    QString outputs[2];
    for(int i=0; i<2; i++)
    {
        line[i]  = in.readLine();
        inputs[i].append(line[i]);
        QStringList partString = inputs[i].split(": ");
        outputs[i] = partString.at(1);
    }
    exomodeInput->setText(outputs[0]);
    displayInput->setText(outputs[1]);
}

MainFrame::~MainFrame()
{
    delete ui;
}

