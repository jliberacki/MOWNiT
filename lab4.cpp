#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <ctime>
#include <stdlib.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

template <typename T>
class FourierTransformation {
public:
    virtual ~FourierTransformation() = default;
    static std::vector<std::complex<T>> DFT_method(std::vector<std::complex<T>> x);
    static std::vector<std::complex<T>> FFT_method(std::vector<std::complex<T>> x);
};

template <typename T>
std::vector<std::complex<T>> FourierTransformation<T>::DFT_method(std::vector<std::complex<T>> x) {
    std::vector<std::complex<T>> res;
    //Dla każdego punktu obliczamy sume elementow Xn zgodnie z formuła Eulera
    for(int i = 0; i < x.size(); i++) {
        std::complex<T> sum = std::complex<T>(0, 0);
        for(int j = 0; j < x.size(); j++) {
            std::complex<T> formula = std::complex<T>(cos(2*M_PI*i*j / x.size()), -sin(2*M_PI*i*j / x.size()));
            sum += x[j] * formula;
        }
        res.push_back(sum);
    }
    return res;
}

template <typename T>
std::vector<std::complex<T>> FourierTransformation<T>::FFT_method(std::vector<std::complex<T>> x) {

    //algorytm oparty na https://jakevdp.github.io/blog/2013/08/28/understanding-the-fft/

    if(x.size() % 2 > 0) {
        throw "Input size must be 2^n";
    }

    //wykonujemy dopoki sie to opłaca, mogłaby to byc wieksza liczba, w przykladzie uzyte jest 32
    if(x.size() <= 2) {
        return DFT_method(x);
    }

    std::vector<std::complex<T>> res(x.size());

    //Dane wejściowe dzielimy na te o elementach parzystych i nieparzystych
    std::vector<std::complex<T>> x_even;
    std::vector<std::complex<T>> x_odd;

    for(int i = 0; i < x.size(); i++){
        if(i % 2 == 0) x_even.push_back(x[i]);
        else x_odd.push_back(x[i]);
    }

    //Rekurencyjnie wywołujemy algorytm FFT na podzielonych wczesniej danych
    x_even = FFT_method(x_even);
    x_odd = FFT_method(x_odd);

    //Dzieki temu iterujemy polowe mniej razy
    for(int i = 0; i < x.size()/2; i++) {
        //Wyznaczamy mnożnik zgodnie ze wzorem
        std::complex<T> factor = std::complex<T>(cos(2*M_PI*i / x.size()), -sin(2*M_PI*i / x.size()));
        //Łaczymy wyniki w calosc korzystajac z policzonego wczesniej mnoznika
        res[i] = x_even[i]+factor*x_odd[i];
        res[i+x.size() / 2] = x_even[i]-factor*x_odd[i];
    }

    return res;
}

int main() {
    try {
        // http://www.twojapogoda.pl/prognoza-16dni-polska/lodzkie-piatek/?page=3
        std::vector<std::complex<double>> x { 
            std::complex<double> (16.0, 0.0),
            std::complex<double> (9.0, 0.0),
            std::complex<double> (10.0, 0.0),
            std::complex<double> (8.0, 0.0),
            std::complex<double> (6.0, 0.0),
            std::complex<double> (2.0, 0.0),
            std::complex<double> (2.0, 0.0),
            std::complex<double> (3.0, 0.0),
            std::complex<double> (2.0, 0.0),
            std::complex<double> (1.0, 0.0),
            std::complex<double> (4.0, 0.0),
            std::complex<double> (3.0, 0.0),
            std::complex<double> (3.0, 0.0),
            std::complex<double> (2.0, 0.0),
            std::complex<double> (2.0, 0.0),
            std::complex<double> (1.0, 0.0) };

        std::vector<std::complex<double>> vect(1024);
        //std::vector<std::complex<double>> vect(2048);
        //std::vector<std::complex<double>> vect(4096);

        //pomiar czasu DFT
        clock_t begin = clock();
        FourierTransformation<double>::DFT_method(vect);
        clock_t end = clock();
        double time_DFT = double(end - begin) / CLOCKS_PER_SEC;

        //pomiar czasu FFT
        begin = clock();
        FourierTransformation<double>::FFT_method(vect);
        end = clock();
        double time_FFT = double(end - begin) / CLOCKS_PER_SEC;

        std::cout<< "DFT - time: " << time_DFT << ", FFT - time: " << time_FFT << std::endl;

        std::vector<std::complex<double>> results_DFT = FourierTransformation<double>::DFT_method(x);
        std::vector<std::complex<double>> results_FFT = FourierTransformation<double>::FFT_method(x);

        //DFT wyniki
        std::cout<< "DFT - results: " << std::endl;
        for(int i = 0; i < results_DFT.size(); ++i) {
            std::cout  << results_DFT[i] << std::endl;
        }

        //FFT wyniki
        std::cout<< "FFT - results: " << std::endl;
        for(int i = 0; i < results_FFT.size(); ++i) {
            std::cout  << results_FFT[i] << std::endl;
        }
    } catch (const char* msg) {
        std::cerr << msg << std::endl;
    }
    
    system ("PAUSE");
}

