#include "mainwindow.h"
#include "ui_mainwindow.h"

//consts
const double pi = 3.1416;
//global
QString NAME = "";
double FD = 0;
double TIME_K = 0;
QList<double> DATA;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
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

int getFileSize(FILE* inFile)
{
    int fileSize = 0;
    fseek(inFile, 0, SEEK_END);

    fileSize = ftell(inFile);

    fseek(inFile, 0, SEEK_SET);
    return fileSize;
}

void MainWindow::read_wav(QString pathname, QList<double> *data, double *time_k)
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
    ui->textBrowser->append( "Header Read " + QString::number(bytesRead) + " bytes." );
    if (bytesRead > 0)
    {
        //Read the data
        uint16_t bytesPerSample = wavHeader.bitsPerSample / 8;      //Number     of bytes per sample
        uint64_t numSamples = wavHeader.ChunkSize / bytesPerSample; //How many samples are in the wav file?
        static const uint16_t BUFFER_SIZE = 4096;
        int16_t *buffer = new int16_t[BUFFER_SIZE];
        while ((bytesRead = fread(buffer, sizeof buffer[0], BUFFER_SIZE / (sizeof buffer[0]), wavFile)) > 0)
        {
            qDebug() <<  "Read " + QString::number(bytesRead) + " bytes";
            for (int i = 0; i < BUFFER_SIZE; i++) {
                data->append(buffer[i]);
            }
        }

        delete [] buffer;
        buffer = nullptr;
        filelength = getFileSize(wavFile);

        //ui->textBrowser->append()
        "File is:" + QString::number(filelength) + " bytes.";
        ui->textBrowser->append( "Data size:" + QString::number(wavHeader.ChunkSize) );

        *time_k = 1.0 / wavHeader.SamplesPerSec;

        // Display the sampling Rate from the header
        ui->textBrowser->append( "Sampling Rate:" + QString::number(wavHeader.SamplesPerSec) );
        ui->textBrowser->append( "Number of bits used:" + QString::number(wavHeader.bitsPerSample ) );
        ui->textBrowser->append( "Number of channels:" + QString::number(wavHeader.NumOfChan ) );
        ui->textBrowser->append( "Number of bytes per second :" + QString::number(wavHeader.bytesPerSec ) );
        ui->textBrowser->append( "Data length:" + QString::number(wavHeader.Subchunk2Size ) );
        ui->textBrowser->append( "Audio Format:" + QString::number(wavHeader.AudioFormat ) );
        // Audio format 1=PCM,6=mulaw,7=alaw, 257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM
        ui->textBrowser->append( "Block align:" + QString::number(wavHeader.blockAlign ) );
    }
    fclose(wavFile);
}

void do_fourier(QList<double> data, amp_freq *AF) {
    double n = data.length();
    QVector<double> X_re(n+1), X_im(n+1), X_amp(n+1), f(n+1);

    for(int k = 1; k < n+1; k++) {
        for (int i = 0; i < n; n++) {
            X_re[k] += data.at(i)*cos( (2*pi*k*i)/n );
            X_im[k] -= data.at(i)*sin( (2*pi*k*i)/n );
        }
        f[k] = double(k) * FD / n; //FD нужно тоже подавать, тк для сэмпла и лонга они разные
        X_amp[k] = sqrt( X_im[k]*X_im[k] + X_re[k]*X_re[k] ) / n;
    }

    AF->amp.append(X_amp);
    AF->freq.append(f);
    AF->re.append(X_re);
    AF->im.append(X_im);
}


void MainWindow::on_open_btn_clicked()
{
    //open dialog
    QString pathname = QFileDialog::getOpenFileName(0, "Открыть файл", "", "*.wav *.WAV");
    QFile file(pathname);
    QFileInfo fileInfo(file.fileName());
    NAME = (fileInfo.fileName());

    //read_wav
    read_wav(pathname, &DATA, &TIME_K);

    //draw
    QLineSeries *data_series = new QLineSeries();

}
