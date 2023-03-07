#ifndef MAINFRAME_H
#define MAINFRAME_H

#include <QWidget>
#include <QLineEdit>
#include <QTextEdit>
#include <QPushButton>
#include <QMap>
#include <QLabel>
#include <QFile>
#include <QVBoxLayout>
#include <QList>
#include <QComboBox>

QT_BEGIN_NAMESPACE
namespace Ui { class MainFrame; }
QT_END_NAMESPACE

class MainFrame : public QWidget
{
    Q_OBJECT

public:
    MainFrame(QWidget *parent = nullptr);
    ~MainFrame();


public slots:
   void openFile();
   void saveFile();
   void newFile();
   void print(QFile*);

private:
    Ui::MainFrame *ui;
    QVBoxLayout *scrollLayout;
    QWidget *scrollWidget;
    int numOfVar;
    //"EXOMODE: ",
    const char *listOfVar[57] = {"EXOMODE: ","DISPLAY: ","LOGFILE: ","OUTDIR: ","NTSAVE: ","DTSAVE: ","MODEL: ","DIMENSIONS: ","METHOD: ",
                                 "FEMORD: ","MESHTYPE: ","SAVEMESH: ","DGTYPE: ","DGPENALTY: ","DGLUMPPING: ","BASISTYPE: ", "DGPENVEL: ",
                                 "ENRICHORD: ","ENRICHTYPE: ","TSMETHOD: ","TMAX: ","CFL: ","FREESURF: ","BC: ","TAPERLEN: ","TAPALPHA: ",
                                 "TAPBETA: ","SRCTYPE: ","SAVESRC: ","SRCFUNC: ","TWIDTH: ","SWIDTH: ","PKFREQ: ","SRCAMP: ","SHOTS: ","SRCLOC: ",
                                 "SRCVECTOR: ","SRCTENSOR: ","SNAPSHRES: ","SNAPSHFMT: ","DXUSEIS: ","NREC: ","REC0: ","RECD: ","NSEISMO: ",
                                 "SEISMO[i]: ","MESHFAC: ","NBLOCKS: ","X_1: ","Y_1: ","Z_1: ","NELEM_1: ","VMAX_1: ","VP_1: ","VS_1: ","RHO_1: ",
                                 "BLOCKTYPE_1: "};
    //Entries
    QWidget *entries[57];
    //QHBoxLayouts
    QHBoxLayout *layouts[57];
    //Inputs
    QLineEdit *lineInputs[57];
    QComboBox *dropInputs[57];
    //Labels
    QLabel *label[57];
    //Buttons
    QPushButton *selectButton;
    QPushButton *saveButton;
    QPushButton *newButton;
};
#endif // MAINFRAME_H
