#ifndef MAINFRAME_H
#define MAINFRAME_H

#include <QWidget>
#include <QLineEdit>
#include <QTextEdit>
#include <QPushButton>
#include <QMap>

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

private:
    Ui::MainFrame *ui;
    QLineEdit *exomodeInput;
    QLineEdit *displayInput;
    QLineEdit *logfileInput;
    QLineEdit *outdirInput;
    QLineEdit *ntsaveInput;
    QLineEdit *dtsaveInput;
    QLineEdit *modelInput;
    QPushButton *selectButton;
    QPushButton *saveButton;
    QPushButton *newButton;
};
#endif // MAINFRAME_H
