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
    //Entries
    QWidget *entries[15];
    //QHBoxLayouts
    QHBoxLayout *layouts[15];
    //Inputs
    QLineEdit *inputs[15];
    //list labels
    QLabel *labelList[15];
    //Labels
    QLabel *label[15];
    //Buttons
    int numOfVar;
    QPushButton *selectButton;
    QPushButton *saveButton;
    QPushButton *newButton;
};
#endif // MAINFRAME_H
