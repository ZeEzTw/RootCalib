#include "task2.h"
#include "sortEnergy.h"
#pragma once
#include <iostream>
#include <fstream>
#include <cmath> // Sau alte biblioteci necesare pentru calcule matematice
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"

float MaxDistance = 10;
float MinDistance = 1.9f;
long int iterations = 0;
long int goodGaus = 0;
long int badGaus = 0;
void fitGaussian(TF1 *gaus, TH1D *hist, double peakX, double rangeWidth)
{
    hist->Fit(gaus, "RQ+", "", peakX - rangeWidth, peakX + rangeWidth);
}

float extractBackgroundAreaForGaussian(TH1D *hist, int peakBin, TF1 *gaussian)
{
    /* float backgroundArea = 0;

     // Obținem limitele stânga și dreapta ale vârfului
     double leftLimit, rightLimit;
     findStartOfPeak(hist, peakBin, gaussian, leftLimit, rightLimit);

     // Calculăm aria totală sub gaussiană între aceste limite
     double totalArea = gaussian->Integral(leftLimit, rightLimit);
     */
}
float calculateAreaOfGaussian(TF1 *gaussian)
{
    float amplitude = gaussian->GetParameter(0);
    float mean = gaussian->GetParameter(1);
    float sd = gaussian->GetParameter(2);

    float area = amplitude * sd * sqrt(2 * TMath::Pi());
    return area;
}

float areaPeak(TF1 *gaussian, TH1D *hist, int peakBin) {
    double mean = gaussian->GetParameter(1);
    double sigma = gaussian->GetParameter(2);
    double leftLimit = mean - 2 * sigma;
    double rightLimit = mean + 2 * sigma;

    // Calculează aria sub gausiană în intervalul [leftLimit, rightLimit]
    double peakArea = gaussian->Integral(leftLimit, rightLimit);

    // Estimează aria fundalului în același interval
    // Se asumă că fundalul este modelat de termenii [3] și [4]*x din funcția gausiană
    TF1 background("background", "[0] + [1]*x", leftLimit, rightLimit);
    background.SetParameters(gaussian->GetParameter(3), gaussian->GetParameter(4));
    double backgroundArea = background.Integral(leftLimit, rightLimit);

    // Calcularea ariei vârfului eliminând contribuția fundalului
    double area = peakArea - backgroundArea;

    return area;
}
float calculateResolution(TF1 *gaussian)
{
    // Parametrul sigma (lățimea standard) al gaussienei
    float sigma = gaussian->GetParameter(2);
    // Parametrul mean (media) al gaussienei
    float mean = gaussian->GetParameter(1);
    // Calcularea FWHM
    float FWHM = 2.3548 * sigma;
    // Calcularea rezoluției
    float resolution = FWHM / mean;
    return resolution;
}

void findStartOfPeak(TH1D *hist, int maxBin, TF1 *gaussian, double &leftLimitPosition, double &rightLimitPosition)
{
    double mean = gaussian->GetParameter(1);
    double sigma = gaussian->GetParameter(2);

    leftLimitPosition = mean - 3 * sigma;
    rightLimitPosition = mean + 3 * sigma;
    double left = abs(leftLimitPosition - maxBin);
    double right = abs(rightLimitPosition - maxBin);
    if(left > right)
    {
        right = left;
    }
    else
    {
        left = right;
    }
    if (left > MaxDistance || left < MinDistance)
    {
        badGaus++;

        goodGaus--;
        leftLimitPosition = maxBin - MinDistance;
    }
    if (right > MaxDistance || right < MinDistance)
    {
        rightLimitPosition = maxBin + MinDistance;
    }
    goodGaus++;
}

void eliminatePeak(TH1D *hist, int maxBin, TF1 *gaussian)
{
    double mean = gaussian->GetParameter(1);
    double sigma = gaussian->GetParameter(2);
    int leftLimit = static_cast<int>(mean - 2 * sigma);
    int rightLimit = static_cast<int>(mean + 2 * sigma);

    // Asigurăm că limitele sunt în interiorul histogramului
    if (leftLimit < 1)
        leftLimit = 1;
    if (rightLimit > hist->GetNbinsX())
        rightLimit = hist->GetNbinsX();
    double leftLimitPosition, rightLimitPosition;
    findStartOfPeak(hist, maxBin, gaussian, leftLimitPosition, rightLimitPosition);
    // Iterăm și setăm conținutul la 0 pentru fiecare bin din interval
    for (int i = leftLimitPosition; i <= rightLimitPosition; i++)
    {
        // cout<<"bin "<<i<<" a devenit 0"<<endl;
        hist->SetBinContent(i, 0);
    }
}
float findPeak(TH1D *hist, int numBins, TH1D *mainHist, int peak, TF1 *gaus[])
{
    float maxPeakY = 0;
    int maxBin = 0;
    double peakWithoutBackground = 0;
    TF1 *gaussianTemp = nullptr;
    double leftLimit, rightLimit;
    leftLimit = MinDistance;
    rightLimit = MinDistance;
    for (int bin = 1; bin <= numBins; ++bin)
    {
        float binContent = hist->GetBinContent(bin);
        if (binContent == 0)
            continue;
        iterations++;
        // Calcularea înălțimii peak-ului fără background
        double peakWithoutBackgroundTemp = binContent - (mainHist->GetBinContent(mainHist->FindBin(bin - leftLimit)) + mainHist->GetBinContent(mainHist->FindBin(bin + rightLimit))) / 2;
        if (peakWithoutBackgroundTemp > peakWithoutBackground)
        {
            peakWithoutBackground = peakWithoutBackgroundTemp;
            maxPeakY = binContent;
            maxBin = bin;
        }
        delete gaussianTemp;
    }

    float maxPeakX = mainHist->GetXaxis()->GetBinCenter(maxBin);
    gaus[peak] = new TF1(Form("gausFit_%d", peak), "[0]*exp(-0.5*((x-[1])/[2])**2) + [3] + ([4]*x)", maxPeakX - 10, maxPeakX + 10);
    gaus[peak]->SetParameters(maxPeakY, maxPeakX, 0.1, 0.0, 0.0); // Setarea parametrilor funcției gaussian
    hist->Fit(gaus[peak], "RQ+");
    eliminatePeak(hist, maxBin, gaus[peak]);
    //<<"peak"<<gaus[peak]->GetParameter(1)<<endl;
    return gaus[peak]->GetParameter(1);
}

void printInFileJson(ofstream &file, int column, int peak, float peak_position, float area, float resolution)
{
    file << "{" << endl;
    file << "  \"Column\": " << column << "," << endl;
    file << "  \"Peak\": " << peak << "," << endl;
    file << "  \"Position\": " << peak_position << "," << endl;
    file << "  \"Area\": " << area << "," << endl;
    file << "  \"Resolution\": " << resolution << endl;
    file << "}" << endl;
}

void sortTxt(double *&energyArray, int &size)
{
    ifstream inputFile("energy.txt");
    ofstream outputFile("output.txt");
    getEnergyArray(inputFile, energyArray, size);
    sortEnergy(energyArray, size);
    printToFile(outputFile, energyArray, size);
}
/*void calculateRegres(double x, int y, double& intercept, double& slope) {
    cout<<"x: "<<x<<" y: "<<y<<endl;
    TGraph *graph = new TGraph(size);
        graph->SetPoint(x, y);

    // Crearea unui obiect TF1 pentru regresia liniară
    TF1 *linearFit = new TF1("linearFit", "pol1"); // Polinom de grad 1 pentru regresia liniară

    // Aplicarea regresiei liniare pe grafic
    graph->Fit(linearFit, "Q"); // "Q" pentru fit rapid, fără afișare grafică

    // Extrage coeficienții dreptei de regresie liniară
    intercept = linearFit->GetParameter(0);
    slope = linearFit->GetParameter(1);

    // Eliberare memorie
    delete graph;
    delete linearFit;

}


// Funcție pentru verificarea corectitudinii regresiei liniare
bool checkRegres(double intercept, double slope, int* peaks, double* energyArray, int size) {
    bool isCorrectlyCorrelated = true;

    for (int i = 0; i < size; ++i) {
        double yCalc = slope * energyArray[i] + intercept;

        // Verifică corectitudinea asociatiei, poți ajusta toleranța aici
        if (fabs(peaks[i] - yCalc) > 1.0) {
            isCorrectlyCorrelated = false;
            break;
        }
    }

    return isCorrectlyCorrelated;
}

// Funcție pentru calibrare
void calibration(int* peaks, int size, double* energyArray) {
    double intercept = 0;
    double slope = 0;

    // Iterăm prin vârfuri și energii pentru a calibra
    for (int i = 0; i < size; ++i) {
        int x = peaks[i];
        cout<<"x: "<<x<<endl;
        // Încercăm mai multe valori de energie (x) pentru a găsi o regresie liniară potrivită
        for (int j = 0; j < size; ++j) {
            double y = energyArray[j];

            // Calculăm regresia liniară
            calculateRegres(x, y, intercept, slope);//problema memorie aici

            // Verificăm corectitudinea regresiei
            //if (checkRegres(intercept, slope, peaks, energyArray, size)) {
                //std::cout << "Calibrated energy for peak " << y << " is " << x << std::endl;
                //break; // Ieșim din bucla internă dacă am găsit o corelație bună
            }
        }
    }


//}*/
double calculateError(double measured[], double known[], int size)
{
    double error = 0.0;
    for (int i = 0; i < size; ++i)
    {
        error += fabs(measured[i] - known[i]);
    }
    return error;
}

// Verifică dacă energia prezisă se află în apropiere de una dintre energiile cunoscute
// Verifică dacă energia prezisă se află în apropiere de una dintre energiile cunoscute
bool checkPredictedEnergies(double predictedEnergy, double knownEnergies[], int size)
{
    double minError = std::numeric_limits<double>::max();

    for (int i = 0; i < size; ++i)
    {
        if (fabs(predictedEnergy - knownEnergies[i]) < minError)
        {
            minError = fabs(predictedEnergy - knownEnergies[i]);
            // cout<<"energyKnown"<<knownEnergies[i]<<endl;
            // cout<<"predictedEnergy"<<predictedEnergy<<endl;
            // cout<<"minError"<<minError<<endl;
        }
    }

    return minError < 10.0;
}
void outputDiferences(double knownEnergies[], int size, int peaks[], int peakCount)
{
    for (int i = 0; i < peakCount; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            cout << "Diferența dintre " << knownEnergies[j] << endl; //<< " și " << peaks[i] << " este " << fabs(knownEnergies[j] - peaks[i]) << std::endl;
        }
    }
    cout << "---------------------------------" << endl;
}
void calibratePeaks(float peaks[], int peakCount, double knownEnergies[], int size, float &M, float &B)
{
    double bestM = 0.0;
    double bestB = 0.0;
    int bestCorrelation = 0;
    // outputDiferences(knownEnergies, size, peaks, peakCount);
    //  Iterăm prin fiecare valoare posibilă pentru panta (m)
    for (double m = 0.1; m <= 10.0; m += 0.01)
    {
        for (double n = 0.0; n <= 5.0; n += 0.1)
        {
            int correlations = 0;

            // Încercăm diferite valori pentru fiecare peak
            for (int i = 0; i < peakCount; ++i)
            {
                // Calculăm energia prezisă pentru panta curentă (m) și intercept-ul (n) curent
                double predictedEnergy = m * peaks[i] + n;

                // Verificăm dacă energia prezisă se potrivește cu una dintre energiile cunoscute
                if (checkPredictedEnergies(predictedEnergy, knownEnergies, size))
                {
                    ++correlations;
                }
            }

            // Actualizăm cea mai bună corelație găsită până acum
            if (correlations > bestCorrelation)
            {
                bestCorrelation = correlations;
                bestM = m; // Actualizăm panta pentru cea mai bună corelație
                bestB = n; // Actualizăm intercept-ul pentru cea mai bună corelație
            }

            // Dacă am găsit deja toate corelațiile posibile, întrerupem bucla
            if (bestCorrelation == peakCount)
            {
                break;
            }
        }

        // Dacă am găsit deja toate corelațiile posibile, întrerupem bucla
        if (bestCorrelation == peakCount)
        {
            break;
        }
    }
    M = bestM;
    B = bestB; 

    // Afișăm rezultatele
    //std::cout << "Best m: " << bestM << ", Best b: " << bestB << ", Best Correlation: " << bestCorrelation << std::endl;
}
bool checkConditions(float peakPosition, float Xmax, float Xmin, float FWHMmax, TF1 *gaus)
{
    bool condition1 = true;
    bool condition2 = true;
    double FWHM = 2.3548 * gaus->GetParameter(2);
    if (peakPosition > Xmax || peakPosition < Xmin)
    {
        condition1 = false;
    }
    if (FWHM > FWHMmax)
    {
        condition2 = false;
    }
    //std::cout << "For peak " << peakPosition << ": condition1 " << condition1 << " condition2 " << condition2 << std::endl;
    return condition1 && condition2;
}
void applyXCalibration(float m, float b, TH1D *originalHist, TH1D *&calibratedHist) {
    int numBinsOriginal = originalHist->GetNbinsX();
    int numBinsCalibrated = calibratedHist->GetNbinsX();

    // Iterează prin fiecare bin al histogramei calibrate
    for (int bin_calibrated = 1; bin_calibrated <= numBinsCalibrated; ++bin_calibrated) {
        // Calculează centrul binului recalibrat
        double binCenter_calibrated = calibratedHist->GetXaxis()->GetBinCenter(bin_calibrated);

        // Calculează poziția corespunzătoare pe axa x originală folosind formula de recalibrare inversă
        double binCenter_original = (binCenter_calibrated - b) / m;

        // Găsește binul corespunzător în histrograma originală folosind interpolarea liniară
        int bin_original = originalHist->GetXaxis()->FindBin(binCenter_original);

        // Adaugă conținutul binului original la binul recalibrat, utilizând metoda de interpolare
        double content_original = originalHist->GetBinContent(bin_original);
        double content_residual = originalHist->GetBinContent(bin_original + 1);
        double fraction = (binCenter_original - originalHist->GetXaxis()->GetBinLowEdge(bin_original)) / originalHist->GetXaxis()->GetBinWidth(bin_original);
        double interpolated_content = content_original + fraction * (content_residual - content_original);
        calibratedHist->SetBinContent(bin_calibrated, interpolated_content);
    }
}

void task1(int number_of_peaks, const char *file_path, double *energyArray, int size, float Xmax, float Xmin, float FWHMmax) {
    std::ofstream file("data.json");
    TFile *fr = new TFile(file_path, "READ");
    TFile *outputFileHistograms = new TFile("histogramsWithPeaks.root", "RECREATE");
    TFile *outputFileCalibrated = new TFile("calibratedHistograms.root", "RECREATE");

    TH2F *h1;
    fr->GetObject("mDelila_raw", h1);
    int number_of_columns = h1->GetNbinsX();

    for (int column = 1; column <= number_of_columns; column++) {
        TH1D *hist1D = h1->ProjectionY(Form("hist1D_col%d", column), column, column);
        if (hist1D->GetMean() < 5) {
            delete hist1D;
            continue;
        }

        float peaks[number_of_peaks];
        TF1 *gaus[number_of_peaks];
        TH1D *tempHist = (TH1D *)hist1D->Clone();
        for (int peak = 0; peak < number_of_peaks; peak++) {
            float peak_position = findPeak(tempHist, hist1D->GetNbinsX(), hist1D, peak, gaus);
            if (checkConditions(peak_position, Xmax, Xmin, FWHMmax, gaus[peak])) {
                peaks[peak] = peak_position;
                hist1D->Fit(gaus[peak], "RQ+");

                // Printează în fișierul JSON
                printInFileJson(file, column, peak, peak_position, areaPeak(gaus[peak], hist1D, peak_position), calculateResolution(gaus[peak]));
            } else {
                if (gaus[peak] != nullptr) {
                    hist1D->GetListOfFunctions()->Remove(gaus[peak]);
                    delete gaus[peak];
                    gaus[peak] = nullptr;
                }
            }
        }
        float M, N;
        calibratePeaks(peaks, number_of_peaks, energyArray, size, M, N);
        cout<<endl<<"M"<<M<<endl;
        cout<<"N"<<N<<endl;
        // Creează histograma calibrată pentru coloana curentă
        TH1D *calibratedHist = new TH1D(Form("calibrated_col%d", column), Form("calibrated_col%d", column), hist1D->GetNbinsX(), hist1D->GetXaxis()->GetXmin(), hist1D->GetXaxis()->GetXmax());
        applyXCalibration(M, N, hist1D, calibratedHist);

        // Resetarea statisticilor și salvarea în fișierul corespunzător
        calibratedHist->ResetStats();
        outputFileHistograms->cd();
        hist1D->Write();
        outputFileCalibrated->cd();
        calibratedHist->Write();

        delete hist1D;
        delete calibratedHist;
    }

    outputFileHistograms->Close();
    outputFileCalibrated->Close();
    fr->Close();
    file.close();

    delete outputFileHistograms;
    delete outputFileCalibrated;
    delete fr;
}
int main(int argc, char **argv)
{
    if (argc < 6)
    {
        std::cerr << "Usage: " << argv[0] << " <number_of_peaks> <file_path> <Xmin> <Xmax> <FWHMmax>" << std::endl;
        return 1;
    }
    int number_of_peaks = std::atoi(argv[1]);
    const char *file_path = argv[2];
    float Xmin = std::atof(argv[3]);
    float Xmax = std::atof(argv[4]);
    float FWHMmax = std::atof(argv[5]);

    double *energyArray = nullptr;
    int size;
    sortTxt(energyArray, size);

    double calibTest[10] = {69.8, 134.8, 199.8, 264.8, 329.8, 394.8, 459.8, 524.8, 589.8, 654.8};
    task1(number_of_peaks, file_path, energyArray, size, Xmax, Xmin, FWHMmax);

    delete[] energyArray;
    std::cout << std::endl
              << "Total number of iterations: " << iterations << std::endl;
    std::cout << "Good Gaussians: " << goodGaus << std::endl;
    std::cout << "Bad Gaussians: " << badGaus << std::endl;
    return 0;
}

/*

void eliminatePeak(TH1D *hist, int maxBin, TF1* gaussian)
{
    //cout<<"Eliminare peak"<<endl;
    double mean = gaussian->GetParameter(1);
    double sigma = gaussian->GetParameter(2);

    // Calculăm limitele inițiale
    int leftLimit = static_cast<int>(mean - 2 * sigma);
    int rightLimit = static_cast<int>(mean + 2 * sigma);

    // Verificăm dacă limitele depășesc o diferență față de bin-ul vârfului
    int bin = maxBin;
    const int maxAllowedDifference = 10;
    const int minAllowedDifference = 4;

    // Verificare pentru stânga
    int leftDifference = abs(leftLimit - bin);
    if (leftDifference > maxAllowedDifference)
    {
        leftLimit = bin - maxAllowedDifference;
    }
    else if (leftDifference < minAllowedDifference)
    {
        leftLimit = bin - minAllowedDifference;
    }

    // Verificare pentru dreapta
    int rightDifference = abs(bin - rightLimit);
    if (rightDifference > maxAllowedDifference)
    {
        rightLimit = bin + maxAllowedDifference;
    }
    else if (rightDifference < minAllowedDifference)
    {
        rightLimit = bin + minAllowedDifference;
    }

    // Asigurăm că limitele sunt în interiorul histogramului
    if (leftLimit < 1) leftLimit = 1;
    if (rightLimit > hist->GetNbinsX()) rightLimit = hist->GetNbinsX();
    // Iterăm și setăm conținutul la 0 pentru fiecare bin din interval
    for (int i = leftLimit; i <= rightLimit; i++)
    {
        //cout<<"bin "<<i<<" a devenit 0"<<endl;
        hist->SetBinContent(i, 0); // i deoarece hist->SetBinContent folosește indexare de la 1
    }
}

void findStartOfPeak(TH1D *hist, int maxBin, TF1 *gaussian, double &leftLimitPosition, double &rightLimitPosition)
{
        double mean = gaussian->GetParameter(1);
        double sigma = gaussian->GetParameter(2);

        leftLimitPosition = mean - 2 * sigma;
        rightLimitPosition = mean + 2 * sigma;
}
/*int findPeak(TH1D *hist, int numBins, TH1D *mainHist, int peak, TF1 *gaus[]) {
    float maxPeakY = 0;
    int maxBin = 0;
    double peakWithoutBackground = 0;
    TF1* guasianTemp = nullptr;
    TH1D *tempHist = (TH1D *)hist->Clone();
    float bin = tempHist->GetMaximumBin();
    double binContent = tempHist->GetBinContent(bin);

    while (binContent > 10) {
        guasianTemp = new TF1("gausTemp", "gaus", bin - 4, bin + 4);
        fitGaussian(guasianTemp, tempHist, bin, 3);
        cout<<"bin"<<bin<<endl;
        double leftLimit, rightLimit;
        findStartOfPeak(mainHist, bin, guasianTemp, leftLimit, rightLimit);

        int leftBin = mainHist->FindBin(leftLimit);
        int rightBin = mainHist->FindBin(rightLimit);

        // Calculul peakWithoutBackgroundTemp folosind tempHist pentru conținutul bine-lor
        double peakWithoutBackgroundTemp = binContent - (tempHist->GetBinContent(leftBin) + tempHist->GetBinContent(rightBin)) / 2;

        if (peakWithoutBackgroundTemp > peakWithoutBackground) {
            peakWithoutBackground = peakWithoutBackgroundTemp;
            maxPeakY = binContent;
            maxBin = bin;
        }
        cout<<"beanSters";
        eliminatePeak(tempHist, bin, guasianTemp); // Eliminare vârf testat din histograma clonată
        delete guasianTemp;
        bin = tempHist->GetMaximumBin(); // Recalcularea bin-ului cu cel mai mare conținut pentru următoarea iterație
        binContent = tempHist->GetBinContent(bin); // Recalcularea conținutului bin-ului
    }

    cout << "maxBin: " << maxBin << endl;

    float maxPeakX = mainHist->GetXaxis()->GetBinCenter(maxBin);
    gaus[peak] = new TF1(Form("gaus%d_peak%d", peak, maxBin), "gaus", maxPeakX - 2.5, maxPeakX + 2.5);
    gaus[peak]->SetParameters(maxPeakY, maxPeakX);
    mainHist->Fit(gaus[peak], "RQ+");

    double leftLimit, rightLimit;
    findStartOfPeak(mainHist, maxBin, gaus[peak], leftLimit, rightLimit);
    eliminatePeak(hist, maxBin, gaus[peak]);

    return maxBin;
}



int findPeak(TH1D *hist, int numBins, TH1D *mainHist, int peak, TF1 *gaus[]) {
    float maxPeakY = 0;
    int maxBin = 0;
    double peakWithoutBackground = 0;
    TF1* guasianTemp; // Declară variabila TF1

    for (int bin = 1; bin <= numBins; ++bin) { // Începe de la 1 pentru a evita bin-ul 0
        float binContent = hist->GetBinContent(bin);
        if (binContent == 0) continue;

        guasianTemp = new TF1("gausTemp", "gaus", bin - 4, bin + 4);
        fitGaussian(guasianTemp, hist, bin, 3);

        double leftLimit, rightLimit;
        findStartOfPeak(mainHist, bin, guasianTemp, leftLimit, rightLimit);
        if(hist->GetBinContent(leftLimit) == 0 || hist->GetBinContent(rightLimit) == 0) continue;
        double peakWithoutBackgroundTemp = binContent - (hist->GetBinContent(leftLimit) + hist->GetBinContent(rightLimit)) / 2;

        if (peakWithoutBackgroundTemp > peakWithoutBackground) {
            //cout<<"leftLimit"<<leftLimit<<endl;
            //cout<<"rightLimit"<<rightLimit<<endl;
            //cout<<"mainHist->LeftLimit"<<hist->GetBinContent(leftLimit)<<endl;
            //cout<<"mainHist->LeftLimit"<<hist->GetBinContent(rightLimit)<<endl;
            peakWithoutBackground = peakWithoutBackgroundTemp;
            maxPeakY = binContent;
            maxBin = bin;
        }
        delete guasianTemp;
    }
    cout<<"maxBin"<<maxBin<<endl;

    float maxPeakX = mainHist->GetXaxis()->GetBinCenter(maxBin);
    gaus[peak] = new TF1(Form("gaus%d_peak%d", peak, maxBin), "gaus", maxPeakX - 2.5, maxPeakX + 2.5); // Nume unic pentru fiecare fit Gaussian
    gaus[peak]->SetParameters(maxPeakY, maxPeakX);
    mainHist->Fit(gaus[peak], "RQ+"); // Ajustare și desenare fit Gaussian pe histogramă

    double leftLimit, rightLimit;
    findStartOfPeak(mainHist, maxBin, gaus[peak], leftLimit, rightLimit);
    eliminatePeak(hist, maxBin, gaus[peak]);

    return maxBin;
}

*/