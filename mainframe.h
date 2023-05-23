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
   void refresh();
   void print(QFile*);
   void printDropItems(QList <QString>, int, int);
   bool checkInteger(QString);
   bool checkPositive(QString);
   bool checkFloat(QString);
   bool checkXYZCordinates(QString);
   bool checkTensor(QString);
private:
    Ui::MainFrame *ui;
    QVBoxLayout *scrollLayout;
    QWidget *scrollWidget;
    //Enum for every variable
    enum enumVar {EXOMODE, DISPLAY, LOGFILE, OUTDIR, NTSAVE, DTSAVE, MODEL, DIMENSIONS, METHOD, FEMORD, MESHTYPE, SAVEMESH, DGTYPE, DGPENALTY,
                 DGLUMPPING, BASISTYPE, DGPENVEL, ENRICHORD, ENRICHTYPE, TSMETHOD, TMAX, CFL, FREESURF, BC, TAPERLEN, TAPALPHA, TAPBETA, SRCTYPE,
                 SAVESRC, SRCFUNC, TWIDTH, SWIDTH, PKFREQ, SRCAMP, SHOTS, SRCLOC, SRCVECTOR, SRCTENSOR, SNAPSHRES, SNAPSHFMT, DXUSEIS, NREC, REC0,
                 RECD, NSEISMO, SEISMOi, MESHFAC, NBLOCKS, NUMVAR};
    //X_, Y_1, Z_1, NELEM_1, VMAX_1, VP_1, VS_1, RHO_1, BLOCKTYPE_1,
    //number of variables
    int numOfVar = NUMVAR;
    int numReps;
    //line counters
    int lineCountLayout = 0;
    int dropCountLayout = 0;
    //line open counter
    int lineCountOpen = 0;
    int dropCountOpen = 0;
    //print line counter
    int lineCountPrint = 0;
    int dropCountPrint = 0;
    //Number of reptitions
    QString nReps;
    //list of variables
    QList <QString> listOfVar {"EXOMODE"/*0*/,"DISPLAY"/*1*/,"LOGFILE"/*2*/,"OUTDIR"/*3*/,"NTSAVE"/*4*/,"DTSAVE"/*5*/,
                               "MODEL"/*6*/,"DIMENSIONS"/*7*/,"METHOD"/*8*/,"FEMORD"/*9*/,"MESHTYPE"/*10*/,"SAVEMESH"/*11*/,
                               "DGTYPE"/*12*/,"DGPENALTY"/*13*/,"DGLUMPPING"/*14*/,"BASISTYPE"/*15*/, "DGPENVEL"/*16*/,
                               "ENRICHORD"/*17*/,"ENRICHTYPE"/*18*/,"TSMETHOD"/*19*/,"TMAX"/*20*/,"CFL"/*21*/,"FREESURF"/*22*/,
                               "BC"/*23*/,"TAPERLEN"/*24*/,"TAPALPHA"/*25*/,"TAPBETA"/*26*/,"SRCTYPE"/*27*/,"SAVESRC"/*28*/,
                               "SRCFUNC"/*29*/,"TWIDTH"/*30*/,"SWIDTH"/*31*/,"PKFREQ"/*32*/,"SRCAMP"/*33*/,"SHOTS"/*34*/,
                               "SRCLOC"/*35*/,"SRCVECTOR"/*36*/,"SRCTENSOR"/*37*/,"SNAPSHRES"/*38*/,"SNAPSHFMT"/*39*/,"DXUSEIS"/*40*/,
                               "NREC"/*41*/,"REC0"/*42*/,"RECD"/*43*/,"NSEISMO"/*44*/,"SEISMO[i]"/*45*/,"MESHFAC"/*46*/,
                               "NBLOCKS"/*47*/,"X_"/*48*/,"Y_"/*49*/,"Z_"/*50*/,"NELEM_"/*51*/,"VMAX_"/*52*/,"VP_"/*53*/,
                               "VS_"/*54*/,"RHO_"/*55*/,"BLOCKTYPE_"/*56*/};
    //Entries arrays
    //EXOMODE
    QList <QString> exomodeEntries {"0","1","2","3"};
    int exomodeSize = exomodeEntries.size();
    //DISPLAY
    QList <QString> displayEntries {"0","1","2"};
    int displaySize = displayEntries.size();
    //MODEL
    QList <QString> modelEntries {"Acustic","Elastic","Frac","Elastoacust","Anisotropic","Aniso+Frac"};
    int modelSize = modelEntries.size();
    //METHOD
    QList <QString> methodEntries {"SGFD","SEM","DG","IGA","VSFD","EG","HEG"};
    int methodSize = methodEntries.size();
    //MESHTYPE
    QList <QString> meshtypeEntries {"Rec","THEX"};
    int meshtypeSize = meshtypeEntries.size();
    //DGTYPE
    QList <QString> dgtypeEntries {"NIPG","SIPG","IIPG","OBB"};
    int dgtypeSize = dgtypeEntries.size();
    //BASISTYPE
    QList <QString> basistypeEntries {"Nodal-Equisp","GLL","Gauss","Legendre"};
    int basistypeSize = basistypeEntries.size();
    //BOOL
    QList <QString> boolEntries {"no","yes"};
    int boolSize = boolEntries.size();
    //TSMETHOD
    QList <QString> tsmethodEntries {"FDM","RK4","LW4"};
    int tsmethodSize = tsmethodEntries.size();
    //BC
    QList <QString> bcEntries {"Neumann","Taper","Periodic"};
    int bcSize = bcEntries.size();
    //SRCTYPE
    QList <QString> srctypeEntries {"pint src","compressional plane wave","shear plane wave","3"};
    int srctypeSize = srctypeEntries.size();
    //SNAPSHFMT
    QList <QString> snapshfmtEntries {"None","Binary","NetCDF","Slice","VTK","Exodus II "};
    int snapshfmtSize = snapshfmtEntries.size();
    //DXUSEIS
    QList <QString> dxuseisEntries {"displacement seismograms","X-derivative seismograms"};
    int dxuseisSize = dxuseisEntries.size();
    //BLOCKTYPE_1
    QList <QString> blocktype_1Entries {"layer","blk"};
    int blocktype_1Size = blocktype_1Entries.size();
    //Entries
    QList <QWidget*> entriesi;
    //QHBoxLayouts
    QList <QHBoxLayout*> layoutsi;
    //Inputs
    QList <QLineEdit*> lineInputsi;
    QList <QComboBox*> dropInputsi;
    //Labels
    QList <QLabel*> labeli;
    //Buttons
    QPushButton *selectButton;
    QPushButton *refreshButton;
    QPushButton *saveButton;
    QPushButton *newButton;
};
#endif // MAINFRAME_H
