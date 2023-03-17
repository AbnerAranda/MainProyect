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
   void printDropItems(const char **, int, int, int);
   bool checkInteger(QString);
   //bool checkPositive(QString);

private:
    Ui::MainFrame *ui;
    QVBoxLayout *scrollLayout;
    QWidget *scrollWidget;
    int numOfVar;
    //list of variables
    const char *listOfVar[57] = {"EXOMODE = "/*0*/,"DISPLAY = "/*1*/,"LOGFILE = "/*2*/,"OUTDIR = "/*3*/,"NTSAVE = "/*4*/,"DTSAVE = "/*5*/,
                                 "MODEL = "/*6*/,"DIMENSIONS = "/*7*/,"METHOD = "/*8*/,"FEMORD = "/*9*/,"MESHTYPE = "/*10*/,"SAVEMESH = "/*11*/,
                                 "DGTYPE = "/*12*/,"DGPENALTY = "/*13*/,"DGLUMPPING = "/*14*/,"BASISTYPE = "/*15*/, "DGPENVEL = "/*16*/,
                                 "ENRICHORD = "/*17*/,"ENRICHTYPE = "/*18*/,"TSMETHOD = "/*19*/,"TMAX = "/*20*/,"CFL = "/*21*/,"FREESURF = "/*22*/,
                                 "BC = "/*23*/,"TAPERLEN = "/*24*/,"TAPALPHA = "/*25*/,"TAPBETA = "/*26*/,"SRCTYPE = "/*27*/,"SAVESRC = "/*28*/,
                                 "SRCFUNC = "/*29*/,"TWIDTH = "/*30*/,"SWIDTH = "/*31*/,"PKFREQ = "/*32*/,"SRCAMP = "/*33*/,"SHOTS = "/*34*/,
                                 "SRCLOC = "/*35*/,"SRCVECTOR = "/*36*/,"SRCTENSOR = "/*37*/,"SNAPSHRES = "/*38*/,"SNAPSHFMT = "/*39*/,"DXUSEIS = "/*40*/,
                                 "NREC = "/*41*/,"REC0 = "/*42*/,"RECD = "/*43*/,"NSEISMO = "/*44*/,"SEISMO[i] =  "/*45*/,"MESHFAC = "/*46*/,
                                 "NBLOCKS = "/*47*/,"X_1 = "/*48*/,"Y_1 = "/*49*/,"Z_1 = "/*50*/,"NELEM_1 = "/*51*/,"VMAX_1 = "/*52*/,"VP_1 = "/*53*/,
                                 "VS_1 = "/*54*/,"RHO_1 = "/*55*/,"BLOCKTYPE_1 = "/*56*/};
    //Entries arrays
    const char *exomodeEntries[4] = {"0","1","2","3"};
    int exomodeSize = sizeof(exomodeEntries)/sizeof(const char*);
    const char *displayEntries[3] = {"0","1","2"};
    int displaySize = sizeof(displayEntries)/sizeof(const char*);
    const char *modelEntries[6] = {"Acustic","Elastic","Frac","Elastoacust","Anisotropic","Aniso+Frac"};
    int modelSize = sizeof(modelEntries)/sizeof(const char*);
    const char *methodEntries[7] = {"SGFD","SEM","DG","IGA","VSFD","EG","HEG"};
    int methodSize = sizeof(methodEntries)/sizeof(const char*);
    const char *meshtypeEntries[2] = {"0","4"};
    int meshtypeSize = sizeof(meshtypeEntries)/sizeof(const char*);
    const char *boolEntries[2] = {"no","yes"};
    int boolSize = sizeof(boolEntries)/sizeof(const char*);
    const char *tsmethodEntries[3] = {"FDM","RK4","LW4"};
    int tsmethodSize = sizeof(tsmethodEntries)/sizeof(const char*);
    const char *bcEntries[3] = {"Neumann","Taper","Periodic"};
    int bcSize = sizeof(bcEntries)/sizeof(const char*);
    const char *srctypeEntries[4] = {"pint src","compressional plane wave","shear plane wave","3"};
    int srctypeSize = sizeof(srctypeEntries)/sizeof(const char*);
    const char *snapshfmtEntries[6] = {"None","Binary","NetCDF","Slice","VTK","Exodus II "};
    int snapshfmtSize = sizeof(snapshfmtEntries)/sizeof(const char*);
    const char *dxuseisEntries[2] = {"displacement seismograms","X-derivative seismograms"};
    int dxuseisSize = sizeof(dxuseisEntries)/sizeof(const char*);
    const char *blocktype_1Entries[2] = {"layer","blk"};
    int blocktype_1Size = sizeof(blocktype_1Entries)/sizeof(const char*);
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
