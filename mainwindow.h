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
    void read_wav(QString pathname, QList<double> *data, double *time_k, double *freq_d);
    void clr_status_text(QString text);

    QtCharts::QLineSeries *filtered_series;
    QtCharts::QChart *spectr_chart;
    QtCharts::QChart *wave_chart;
    QtCharts::QLineSeries *corr_series;
    QtCharts::QLineSeries *amp_series;
    QtCharts::QLineSeries *data_series;
    QtCharts::QLineSeries *SPECTR_SERIES;

private slots:
    void on_open_btn_clicked();

    void on_fourier_btn_clicked();

    void on_filter_btn_clicked();

    void on_filtered_checkBox_stateChanged(int arg1);

    void on_corr_btn_clicked();

    void on_amps_checkBox_stateChanged(int arg1);

    void on_corr_checkBox_stateChanged(int arg1);

    void on_wave_checkBox_stateChanged(int arg1);

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H
