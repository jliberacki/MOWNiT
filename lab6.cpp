#include <functional>
#include <random>
#include <math.h>
#include <iostream>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double function1(double x) {
    return x*x*x;
}

double function2(double x) {
    return x*x + 2 / (x*x*x);
}

double function3(double x) {
    return 3*x*x*x*x*x + 5*x*x*x + 2*x*x;
}

double rectangles(double xp, double xk, double n, const std::function<double(double)>& function) {
    double sum = 0.0;

    for(int i = 1; i <= n; i++) {
        sum += function(xp + i * (xk-xp) / n);
    }
    
    return (xk-xp) * sum / n;
}

double trapezes(double xp, double xk, double n, const std::function<double(double)>& function) {
  double sum = 0.0;

  for(int i = 0; i < n; i++) {
    sum += (function(xp + i * (xk-xp) / n) + function(xp + (i + 1.0) * (xk-xp) / n)) / 2.0;
  }
  return (xk-xp) / n * sum;
}

// based on https://eduinf.waw.pl/inf/alg/004_int/0004.php
double simpsons(double xp, double xk, double n, const std::function<double(double)>& function) {
    double sum = function(xp) + function(xk);
    double position = xp + (xk-xp) / n;
    double t = position - (xk-xp) / n / 2.0;

    for(int i = 1; i < n; i++) {
        sum += 2 * function(position) + 4 * function(t);
        position += (xk-xp) / n;
        t += (xk-xp) / n;
    }

    sum += 4 * function(t);
    return (xk-xp) * sum / (6 * n);
}

double measure_pi(int n) {

    //https://stackoverflow.com/questions/2704521/generate-random-double-numbers-in-c
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::default_random_engine re;

    int missed = 0;
    for(int i = 0; i < n; i++) {
        double x = unif(re);
        double y = unif(re);
        if(sqrt(x * x + y * y) < 1.0) missed++;
    }
    return 4.0 * (double)missed / (double)n;
}

double monte_carlo(double xp, double xk, int n, const std::function<double(double)> &function) {
    double sum = 0.0;

    //https://stackoverflow.com/questions/2704521/generate-random-double-numbers-in-c
    std::uniform_real_distribution<double> unif(xp, xk);
    std::default_random_engine re;

    for(int i = 0; i < n; i++) {
        sum += function(unif(re));
    }
    return (xk-xp) * sum / n;
}

int main() {

    // double result = 20.0;
    // int xp = 1;
    // int xk = 3;

    // double result = 37.0/12.0;
    // int xp = 1;
    // int xk = 2;

    double result = 225.0 / 4.0;
    int xp = -1;
    int xk = 2;

    for(int i = 1; i < 100; i += 5) {
        //std::cout << i << " " << fabs(result - rectangles(xp,xk,i,function3)) << std::endl;
        //std::cout << i << " " << fabs(result - trapezes(xp,xk,i,function3)) << std::endl;
        //std::cout << i << " " << fabs(result - simpsons(xp,xk,i,function3)) << std::endl;
        std::cout << i << " " << fabs(result - monte_carlo(xp,xk,i,function3)) << std::endl;
    }

    
    // for(int i = 100; i < 10000; i += 100) {
    //     std::cout << i << " " << fabs(M_PI - measure_pi(i)) << std::endl;
    // }

    system ("PAUSE");
}