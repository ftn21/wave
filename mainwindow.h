#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <complex>
#include <QtGlobal>

#include <QtCharts/QChartView>
#include <QtCharts/QBarSeries>
#include <QtCharts/QBarSet>
#include <QtCharts/QLegend>
#include <QtCharts/QBarCategoryAxis>
#include <QtCharts/QHorizontalStackedBarSeries>
#include <QtCharts/QLineSeries>
#include <QtCharts/QCategoryAxis>

#include <iostream>
#include <string>
#include <fstream>
#include <cstdint>

using std::cin;
using std::cout;
using std::endl;
using std::fstream;
using std::string;

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    void read_wav(QString pathname, QList<double> *data, double *freq_dis);
    void clr_status_text(QString text);

private slots:
    void on_open_btn_clicked();

    void on_fourier_btn_clicked();

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
