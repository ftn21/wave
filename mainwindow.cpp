#include "mainwindow.h"
#include "ui_mainwindow.h"

//"consts"
const double pi = 3.1416;
const int BLOCK = 1024;
QStringList html_colours = { "#DA70D6" , "#FF00FF" , "#FF00FF" , "#BA55D3" , "#9370DB" , "#8A2BE2" ,
                             "#9400D3" , "#9932CC" , "#8B008B" , "#800080" , "#4B0082" , "#6A5ACD" ,
                             "#483D8B",  "#3CB371" , "#2E8B57" , "#228B22" , "#008000" , "#006400" ,
                             "#9ACD32" , "#6B8E23" , "#808000" , "#556B2F" , "#66CDAA" , "#8FBC8F" ,
                             "#20B2AA" , "#008B8B" , "#008080" , "#5F9EA0" , "#4682B4" , "#0000FF" ,
                             "#0000CD" , "#00008B" , "#000080" , "#191970" , "#DAA520" , "#B8860B" ,
                             "#CD853F" , "#D2691E" };
//global
QString NAME = "";
double FD = 0;
double TIME_K = 0;
QList<double> DATA;
QList<double> T;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ui->filtered_checkBox->hide();
    ui->amps_checkBox->hide();
    ui->corr_checkBox->hide();
    ui->spectr_checkBox->hide();
    ui->wave_checkBox->hide();
}

MainWindow::~MainWindow()
{
    delete ui;
}

struct amp_freq {
    QVector<double> amp;
    QVector<double> freq;
    QVector<double> im;
    QVector<double> re;
};

//one more global
amp_freq ampfreq;

typedef struct  WAV_HEADER
{
    /* RIFF Chunk Descriptor */
    uint8_t         RIFF[4];        // RIFF Header Magic header
    uint32_t        ChunkSize;      // RIFF Chunk Size
    uint8_t         WAVE[4];        // WAVE Header
    /* "fmt" sub-chunk */
    uint8_t         fmt[4];         // FMT header
    uint32_t        Subchunk1Size;  // Size of the fmt chunk
    uint16_t        AudioFormat;    // Audio format 1=PCM,6=mulaw,7=alaw,     257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM
    uint16_t        NumOfChan;      // Number of channels 1=Mono 2=Sterio
    uint32_t        SamplesPerSec;  // Sampling Frequency in Hz
    uint32_t        bytesPerSec;    // bytes per second
    uint16_t        blockAlign;     // 2=16-bit mono, 4=16-bit stereo
    uint16_t        bitsPerSample;  // Number of bits per sample
    /* "data" sub-chunk */
    uint8_t         Subchunk2ID[4]; // "data"  string
    uint32_t        Subchunk2Size;  // Sampled data length
} wav_hdr;

//one more global
wav_hdr HEADER;

int getFileSize(FILE* inFile)
{
    int fileSize = 0;
    fseek(inFile, 0, SEEK_END);

    fileSize = ftell(inFile);

    fseek(inFile, 0, SEEK_SET);
    return fileSize;
}

void MainWindow::clr_status_text(QString text) {
    QString clr = html_colours.at(qrand() % (html_colours.length() + 1));
    ui->statusBar->setStyleSheet("color: " + clr);
    ui->statusBar->showMessage(text);
}

void MainWindow::read_wav(QString pathname, QList<double> *data, double *time_k, double *freq_d)
{
    wav_hdr wavHeader;
    int headerSize = sizeof(wav_hdr), filelength = 0;

    string fp = pathname.toStdString();
    const char* filePath = fp.c_str();

    FILE* wavFile = fopen(filePath, "r");
    if (wavFile == nullptr)
    {
        fprintf(stderr, "Unable to open wave file: %s\n", filePath);
    }

    //Read the header
    size_t bytesRead = fread(&wavHeader, 1, headerSize, wavFile);
    //one more for later saving
    fread(&HEADER, 1, headerSize, wavFile);
    if (bytesRead > 0)
    {
        //Read the data
        uint16_t bytesPerSample = wavHeader.bitsPerSample / 8;      //Number of bytes per sample
        uint64_t numSamples = wavHeader.ChunkSize / bytesPerSample; //How many samples are in the wav file?
        static const uint16_t BUFFER_SIZE = 2048;
        int16_t *buffer = new int16_t[BUFFER_SIZE];
        while ((bytesRead = fread(buffer, sizeof buffer[0], BUFFER_SIZE / (sizeof buffer[0]), wavFile)) > 0)
        {
            //info
            qDebug() <<  "Read " + QString::number(bytesRead) + " bytes";

            //fill data array
            for (int i = 0; i < int(bytesRead); i++) {
                data->append(buffer[i]/32768.0);
            }
        }

        delete [] buffer;
        buffer = nullptr;
        filelength = getFileSize(wavFile);

        ui->textBrowser->append("<B>" + pathname + "</B>");
        ui->textBrowser->append("");
        ui->textBrowser->append("File is:  " + QString::number(filelength) + " bytes");
        ui->textBrowser->append( "Data size:  " + QString::number(wavHeader.ChunkSize) );

        // 1/samplingRate and samplingrate
        *time_k = 1.0 / wavHeader.SamplesPerSec;
        *freq_d = wavHeader.SamplesPerSec;

        // Display the sampling Rate from the header
        ui->textBrowser->append( "Sampling Rate:  " + QString::number(wavHeader.SamplesPerSec) );
        ui->textBrowser->append( "Number of bits used:  " + QString::number(wavHeader.bitsPerSample ) );
        ui->textBrowser->append( "Number of channels:  " + QString::number(wavHeader.NumOfChan ) );
        ui->textBrowser->append( "Number of bytes per second :  " + QString::number(wavHeader.bytesPerSec ) );
        ui->textBrowser->append( "Data length:  " + QString::number(wavHeader.Subchunk2Size ) );
        ui->textBrowser->append( "Audio Format:  " + QString::number(wavHeader.AudioFormat ) );
        // Audio format 1=PCM,6=mulaw,7=alaw, 257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM
        //ui->textBrowser->append( "Block align:  " + QString::number(wavHeader.blockAlign ) );
    }
    fclose(wavFile);
}

void MainWindow::save_wav(QList<double> data, double time_k, QString pathname) {
        wav_hdr wavHeader;
        int headerSize = sizeof(wav_hdr), filelength = 0;

        char c = 'R';
        char *pc = &c;
        wavHeader.RIFF[0] = *pc;
        c = 'I';
        wavHeader.RIFF[1] = *pc;
        c = 'F';
        wavHeader.RIFF[2] = *pc;
        c = 'F';
        wavHeader.RIFF[3] = *pc;

        c = 'W';
        wavHeader.WAVE[0] = *pc;
        c = 'A';
        wavHeader.WAVE[1] = *pc;
        c = 'V';
        wavHeader.WAVE[2] = *pc;
        c = 'E';
        wavHeader.WAVE[3] = *pc;

        c = 'f';
        wavHeader.fmt[0] = *pc;
        c = 'm';
        wavHeader.fmt[1] = *pc;
        c = 't';
        wavHeader.fmt[2] = *pc;
        c = ' ';
        wavHeader.fmt[3] = *pc;

        wavHeader.Subchunk1Size = 16;
        wavHeader.AudioFormat = 1;
        wavHeader.NumOfChan = 1;
        wavHeader.SamplesPerSec = round(1 / time_k);

        wavHeader.blockAlign = 2;
        wavHeader.bytesPerSec= wavHeader.SamplesPerSec * wavHeader.blockAlign;

        /*
        wavHeader.bytesPerSec = 2;
        wavHeader.blockAlign = wavHeader.SamplesPerSec * wavHeader.blockAlign;
        */

        wavHeader.bitsPerSample = wavHeader.blockAlign*8;

        c = 'd';
        wavHeader.Subchunk2ID[0] = *pc;
        c = 'a';
        wavHeader.Subchunk2ID[1] = *pc;
        c = 't';
        wavHeader.Subchunk2ID[2] = *pc;
        c = 'a';
        wavHeader.Subchunk2ID[3] = *pc;

        wavHeader.Subchunk2Size = wavHeader.blockAlign*data.length();
        wavHeader.ChunkSize = wavHeader.Subchunk2Size + 36;

        cout << "File is                    :" << filelength << " bytes." << endl;
        cout << "RIFF header                :" << wavHeader.RIFF[0] << wavHeader.RIFF[1] << wavHeader.RIFF[2] << wavHeader.RIFF[3] << endl;
        cout << "WAVE header                :" << wavHeader.WAVE[0] << wavHeader.WAVE[1] << wavHeader.WAVE[2] << wavHeader.WAVE[3] << endl;
        cout << "FMT                        :" << wavHeader.fmt[0] << wavHeader.fmt[1] << wavHeader.fmt[2] << wavHeader.fmt[3] << endl;
        cout << "Data size                  :" << wavHeader.ChunkSize << endl;
        cout << "Sampling Rate              :" << wavHeader.SamplesPerSec << endl;
        cout << "Number of bits used        :" << wavHeader.bitsPerSample << endl;
        cout << "Number of channels         :" << wavHeader.NumOfChan << endl;
        cout << "Number of bytes per second :" << wavHeader.bytesPerSec << endl;
        cout << "Data length                :" << wavHeader.Subchunk2Size << endl;
        cout << "Audio Format               :" << wavHeader.AudioFormat << endl;
        cout << "Block align                :" << wavHeader.blockAlign << endl;
        cout << "Data string                :" << wavHeader.Subchunk2ID[0] << wavHeader.Subchunk2ID[1] << wavHeader.Subchunk2ID[2] << wavHeader.Subchunk2ID[3] << endl;

        string fp = pathname.toStdString();
        const char* filePath = fp.c_str();

        FILE* wavFile = fopen(filePath, "w");
        fwrite(&wavHeader, 1, headerSize, wavFile);

        short int w = 0;
        for (int i = 0; i < data.length(); i++) {
            w = round(data[i]*32768);
            fwrite(&w, 2 , 1, wavFile);
        }
        cout << sizeof(w);

        fclose(wavFile);
}

void do_fourier(QList<double> data, double freq_d, amp_freq *AF) {
    double n = data.length();
    QVector<double> X_re(n+1), X_im(n+1), X_amp(n+1), f(n+1);

    for (int k = 1; k < n+1; k++) {
        for (int i = 0; i < n; i++) {
            X_re[k] = X_re[k] + data.at(i)*cos( (2*pi*k*i)/n );
            X_im[k] = X_im[k] - data.at(i)*sin( (2*pi*k*i)/n );
        }
        f[k] = double(k) * freq_d / n;
        X_amp[k] = sqrt( X_im[k]*X_im[k] + X_re[k]*X_re[k] ) / n;
    }

    AF->amp.append(X_amp);
    AF->freq.append(f);
    AF->re.append(X_re);
    AF->im.append(X_im);
}

void fill_series(QLineSeries *series, QVector<double> X, QVector<double> Y) {
    for (int i = 0; i < X.length(); i++) {
        series->append(X[i], Y[i]);
    }
}

void do_filter(amp_freq *AF, double limit) {
    double n = AF->amp.length();
    for (int i = 0; i < n; i++) {
        if (AF->amp.at(i) < limit) {
            AF->amp.replace(i, 0);
        }
    }
}


void MainWindow::on_open_btn_clicked()
{
    //open dialog
    QString pathname = QFileDialog::getOpenFileName(0, "Открыть файл", "", "*.wav *.WAV");
    QFile file(pathname);
    QFileInfo fileInfo(file.fileName());
    NAME = (fileInfo.fileName());

    //read_wav
    read_wav(pathname, &DATA, &TIME_K, &FD);

    //draw
    data_series = new QLineSeries();
    data_series->setName(NAME);

    //fill series
    for (int i = 0; i < DATA.length(); i++) {
        T.append(i*TIME_K);
        data_series->append(T[i], DATA[i]);
    }

    //info
    clr_status_text(NAME + " succsesfully loaded");

    //set chart
    wave_chart = new QChart();
    wave_chart->legend()->setVisible(true);
    wave_chart->legend()->setAlignment(Qt::AlignBottom);
    wave_chart->addSeries(data_series);
    wave_chart->createDefaultAxes();
    wave_chart->setTheme(QChart::ChartThemeBrownSand);
    wave_chart->setAnimationOptions(QChart::NoAnimation);

    //adding gui elements
    ui->wave_view->setChart(wave_chart);
    ui->wave_checkBox->show();
}

void MainWindow::on_fourier_btn_clicked()
{
    //if a file not big
    if (DATA.length() < 1600 /*2*BLOCK*/) {
        do_fourier(DATA, FD, &ampfreq);
    }
    else {
        //if a file BIG
        int count_block = DATA.length() / BLOCK;

        //all blocks
        for (int i = 0; i < count_block; i++) {

            //forming a block
            QList<double> block;
            for (int j = 0; j < BLOCK; j++) {
                block.append(DATA[i*BLOCK + j]);
            }

            //fourier for each block
            do_fourier(block, FD, &ampfreq);

            qDebug() << "block " + QString::number(i) + " processed";
            //clr_status_text("block " + QString::number(i) + " processed");
        }
    }
    clr_status_text("fourier done");

    //draw
    SPECTR_SERIES = new QLineSeries();
    fill_series(SPECTR_SERIES, ampfreq.freq, ampfreq.amp);
    SPECTR_SERIES->setName(NAME);

    //set chart
    spectre_chart = new QChart();
    spectre_chart->legend()->setVisible(true);
    spectre_chart->legend()->setAlignment(Qt::AlignBottom);
    spectre_chart->addSeries(SPECTR_SERIES);
    spectre_chart->createDefaultAxes();
    spectre_chart->setTheme(QChart::ChartThemeBrownSand);
    spectre_chart->setAnimationOptions(QChart::NoAnimation);

    //adding gui elements
    ui->spectr_view->setChart(spectre_chart);
    ui->fourier_btn->setEnabled(false);
    ui->spectr_checkBox->show();
}

void MainWindow::on_filter_btn_clicked()
{
    //filter
    do_filter(&ampfreq, 0.01);
    clr_status_text("data filtered");

    //draw
    filtered_series = new QLineSeries();
    filtered_series->setName("filtered " + NAME);
    fill_series(filtered_series, ampfreq.freq, ampfreq.amp);

    //adding gui elements
    spectre_chart->addSeries(filtered_series);
    ui->spectr_view->repaint();
    ui->filtered_checkBox->setEnabled(true);
    ui->filtered_checkBox->show();
    ui->filter_btn->setEnabled(false);
}

void MainWindow::on_filtered_checkBox_stateChanged(int arg1)
{
    if (arg1 == 2) {
        clr_status_text("filtered series are showed");
        spectre_chart->addSeries(filtered_series);
    }
    else {
        clr_status_text("filtered series are hidden");
        spectre_chart->removeSeries(filtered_series);
    }
}

void MainWindow::on_corr_btn_clicked()
{
    corr_series = new QLineSeries();
    amp_series = new QLineSeries();

    QList<double> data_sample;
    double time_k_sample, fd_sample;
    double mx, my, dx, dy, k;

    read_wav("/home/mitya/Documents/MAI/5sem/StatDin/clave_sample.wav", &data_sample, &time_k_sample, &fd_sample);

    double divk = time_k_sample / TIME_K;
    double count = (DATA.length() + 1) / divk - data_sample.length();
    double n = data_sample.length();

    //мат ожидание
    for (int i = 0; i < count; i++) {
        mx = 0;
        my = 0;
        for (int j = 0; j < n; j++) {
            mx += DATA.at(divk*i + j);
            my += data_sample.at(j);
        }
        mx /= n;
        my /= n;
        amp_series->append(time_k_sample*i, DATA.at(divk*i));
    }

    //дисперсия
    for (int i = 0; i < count; i++) {
        dx = 0;
        dy = 0;
        for (int j = 0; j < n; j++) {
            dx += ( DATA.at(divk*i + j) - mx ) * ( DATA.at(divk*i + j) - mx );
            dy += ( data_sample.at(j) - my ) * ( data_sample.at(j) - my );
        }
        dx /= n;
        dy /= n;
    }

    //корреляция
    QList<double> corr;
    for (int i = 0; i < count; i++) {
        k = 0;
        for (int j = 0; j < n; j++) {
            k += ( DATA.at(divk*i + j) - mx ) * ( data_sample.at(j) - my );
        }
        k /= sqrt(dx*dy)*n;
        corr.append(k);
        corr_series->append(time_k_sample*i, k);
    }

    //setting chart & adding gui elements
    amp_series->setName("amps");
    corr_series->setName("correlation");
    corr_chart = new QChart();
    corr_chart->addSeries(amp_series);
    corr_chart->addSeries(corr_series);
    corr_chart->legend()->setVisible(true);
    corr_chart->legend()->setAlignment(Qt::AlignBottom);
    corr_chart->createDefaultAxes();
    corr_chart->setTheme(QChart::ChartThemeBrownSand);
    corr_chart->setAnimationOptions(QChart::NoAnimation);

    ui->corr_view->setChart(corr_chart);

    //btns&checkBoxes
    ui->corr_btn->setEnabled(false);
    ui->corr_checkBox->show();
    ui->amps_checkBox->show();
}

void MainWindow::on_amps_checkBox_stateChanged(int arg1)
{
    if (arg1 == 2) {
        clr_status_text("amps series are showed");
        corr_chart->addSeries(amp_series);
    }
    else {
        clr_status_text("amp series are hidden");
        corr_chart->removeSeries(amp_series);
    }
}

void MainWindow::on_corr_checkBox_stateChanged(int arg1)
{
    if (arg1 == 2) {
        clr_status_text("correlation series are showed");
        corr_chart->addSeries(corr_series);
    }
    else {
        clr_status_text("correlation series are hidden");
        corr_chart->removeSeries(corr_series);
    }
}

void MainWindow::on_wave_checkBox_stateChanged(int arg1)
{
    if (arg1 == 2) {
        clr_status_text("wave series are showed");
        wave_chart->addSeries(data_series);
    }
    else {
        clr_status_text("wave series are hidden");
        wave_chart->removeSeries(data_series);
    }
}


void MainWindow::on_proc_btn_clicked()
{
    QList<double> data_sample;
    double time_k_sample, fd_sample;

    read_wav("/home/mitya/Documents/MAI/5sem/StatDin/clave_sample.wav", &data_sample, &time_k_sample, &fd_sample);

    int count = ceil(0.0377/time_k_sample);

    QVector<double> data1(count), data2(count);

    double v1, v2;
    for (int i = 0; i < count; i++) {
        v1 = data_sample.at(i) * sin( 2*pi*(i*time_k_sample) * 1297 );
        v2 = data_sample.at(i) * sin( 2*pi*(i*time_k_sample) * 3891 );

        data1[i] = v1;
        data2[i] = v2;
    }

    int st = 2*round(1/(3891*time_k_sample));

    int w1 = 4*round(1/(1297*time_k_sample));
    int w2 = 4*round(1/(3891*time_k_sample));

    int n = (count - 2*w1) / st;

    QVector<double> fr1_A(n), fr2_A(n), time_list(n);

    for (int i = 0; i < n; i++) {
        v1 = 0;
        for (int j = -w1; j < w1; j++) {
            v1 += data1.at(w1 + i*st + j);
        }
        v1 /= 2*w1 + 1;
        fr1_A[i] = v1;

        v2 = 0;
        for (int j = -w2; j < w2; j++) {
            v2 += data2.at(w2 + i*st + j);
        }
        v2 /= 2*w2 + 1;
        fr2_A[i] = v2;

        time_list[i] = (i*st*time_k_sample);
    }

    //set series
    QLineSeries *fr1_series = new QLineSeries();
    QLineSeries *fr2_series = new QLineSeries();
    fr1_series->setName("fr1");
    fr2_series->setName("fr2");

    fill_series(fr1_series, time_list, fr1_A);
    fill_series(fr2_series, time_list, fr2_A);

    //set chart
    proc_chart = new QChart();
    proc_chart->legend()->setVisible(true);
    proc_chart->legend()->setAlignment(Qt::AlignBottom);
    proc_chart->addSeries(fr1_series);
    proc_chart->addSeries(fr2_series);
    proc_chart->createDefaultAxes();
    proc_chart->setTheme(QChart::ChartThemeBrownSand);
    proc_chart->setAnimationOptions(QChart::NoAnimation);

    //adding gui elements
    ui->proc_view->setChart(proc_chart);
    ui->proc_btn->setEnabled(false);
    //ui->spectr_checkBox->show();

}

void MainWindow::on_spectr_checkBox_stateChanged(int arg1)
{
    if (arg1 == 2) {
        clr_status_text("spectre series are showed");
        spectre_chart->addSeries(SPECTR_SERIES);
    }
    else {
        clr_status_text("spectre series are hidden");
        spectre_chart->removeSeries(SPECTR_SERIES);
    }
}

void MainWindow::on_prepare_btn_clicked()
{

}

void MainWindow::on_save_btn_clicked()
{
    QString pathname = QFileDialog::getSaveFileName(0, "save file");
    save_wav(DATA, TIME_K, pathname);
}
