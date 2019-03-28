#include <iostream>
#include <vector>
#include <array>
#include <math.h>
#include <algorithm>
#include <stdlib.h>

typedef std::vector<std::array<float, 2>> PointsList2D;
typedef std::vector<float> PointsList1D;
/// <summary>
/// Template function used for clamping specified datasets
/// </summary>
template <typename T>
T GetIndexClamped(const std::vector<T> points, int index)
{
	if (index < 0)
		return points[0];
	else if (index >= int(points.size()))
		return points.back();
	else
		return points[index];
}

/// <summary>
/// Defines structure for interpolation classes
/// </summary>
class IInterpolation
{
	// public construction and destruction methods
	public:
		virtual ~IInterpolation() = default;

	// public interface methods
	public:
        virtual void Interpolate1D(int pointsToInterpolate) =  0;
		virtual void Interpolate2D(int pointsToInterpolate) =  0;

};

/// <summary>
/// Defines Hermite Cubic Interpolation 
/// </summary>
class CubicInterpolation : public IInterpolation
{
	// public construction methods
	public:
		 CubicInterpolation(const PointsList2D& points2D, const PointsList1D& points1D) : pointsList2D(points2D), pointsList1D(points1D) {}
    // IInterpolation methods
    public:
        void Interpolate2D(int pointsToInterpolate) override
		{
			std::vector<int> index(pointsToInterpolate);
			std::vector<float> t;
			std::vector<float> tx;

			int i = 0, points_size = pointsList2D.size() - 1;
			std::generate(index.begin(), index.end(), [&i, &pointsToInterpolate, &points_size, &t, &tx]()
			{
				float percent = ((float)i) / (float(pointsToInterpolate - 1));
				tx.push_back((points_size)* percent);
				t.push_back(tx[i] - floor(tx[i]));
				return int(tx[i++]);
			});

			for (int i = 0; i < pointsToInterpolate; ++i)
			{
				std::array<PolynomialCoeffs, 2> coeffs;
				std::array<float, 2> A = GetIndexClamped(pointsList2D, index[i] - 1);
				std::array<float, 2> B = GetIndexClamped(pointsList2D, index[i]);
				std::array<float, 2> C = GetIndexClamped(pointsList2D, index[i] + 1);
				std::array<float, 2> D = GetIndexClamped(pointsList2D, index[i] + 2);

				for (int i = 0; i < 2; i++)
				{
					coeffs[i].A = A[i];
					coeffs[i].B = B[i];
					coeffs[i].C = C[i];
					coeffs[i].D = D[i];
				}

				float x = CubicHermite(coeffs[0], t[i]);
				float y = CubicHermite(coeffs[1], t[i]);

				std::cout << "Value at " << tx[i] << " = " << x << "  " << y << std::endl;
			}
		};

        void Interpolate1D(int pointsToInterpolate) override 
		{
            std::vector<int> index(pointsToInterpolate);
			std::vector<float> t;
			std::vector<float> tx;

			int i = 0, points_size = pointsList1D.size() - 1;
			std::generate(index.begin(), index.end(), [&i, &pointsToInterpolate, &points_size, &t, &tx]()
			{
				float percent = ((float)i) / (float(pointsToInterpolate - 1));
				tx.push_back((points_size)* percent);
				t.push_back(tx[i] - floor(tx[i]));
				return int(tx[i++]);
			});

			for (int i = 0; i < pointsToInterpolate; ++i)
			{
				PolynomialCoeffs coeffs;
				float A = GetIndexClamped(pointsList1D, index[i] - 1);
				float B = GetIndexClamped(pointsList1D, index[i]);
				float C = GetIndexClamped(pointsList1D, index[i] + 1);
				float D = GetIndexClamped(pointsList1D, index[i] + 2);

					coeffs.A = A;
					coeffs.B = B;
					coeffs.C = C;
					coeffs.D = D;

				float x = CubicHermite(coeffs, t[i]);

				std::cout << "Value at " << tx[i] << " = " << x << std::endl;
			}
        }
        
    // private methods
	private:
		struct PolynomialCoeffs
		{
			float A, B, C, D;
		};

        float CubicHermite(PolynomialCoeffs coeffs, float t) const
		{
			float a = -coeffs.A / 2.0f + (3.0f*coeffs.B) / 2.0f - (3.0f*coeffs.C) / 2.0f + coeffs.D / 2.0f;
			float b = coeffs.A - (5.0f*coeffs.B) / 2.0f + 2.0f*coeffs.C - coeffs.D / 2.0f;
			float c = -coeffs.A / 2.0f + coeffs.C / 2.0f;
			float d = coeffs.B;

			return a * pow(t, 3) + b * pow(t, 2) + c * t + d;
		}
        
	// private members
	private:
		const PointsList2D& pointsList2D;
		const PointsList1D& pointsList1D;

};

class LagrangeInterpolation : public IInterpolation {
public:
    LagrangeInterpolation(const PointsList1D& points1D) : pointsList1D(points1D) {}

    void Interpolate1D(int pointsToInterpolate) override
	{
		std::vector<int> index(pointsToInterpolate);
		std::vector<float> t;

		int i = 0, points_size = pointsList1D.size() - 1;
		std::generate(index.begin(), index.end(), [&i, &pointsToInterpolate, &points_size, &t]()
		{
			float percent = ((float)i) / (float(pointsToInterpolate - 1));
			t.push_back((points_size)* percent);
			return int(t[i++]);
		});

		for(int i = 0; i < pointsToInterpolate; i++)
		{
			float x = 0.0;
			for(int j = 0; j < pointsList1D.size(); j++)
			{
				float basePolynomial = 1.0;
				for(int k = 0; k < pointsList1D.size(); k++)
				{
					if(j != k)
						basePolynomial *= (t[i] - (float(k))) / (float(j - k));
				}
				x += basePolynomial * pointsList1D[j];
    		}
			std::cout << "Value at " << t[i] << " = " << x << std::endl;
		}
	}

	void Interpolate2D(int pointsToInterpolate) override
        {
            
        };

private:
    const PointsList1D& pointsList1D;
    float CalculateLagrange1D(float x);
};


const PointsList1D points1D =
{
    0.0f,
    1.6f,
    2.3f,
    3.5f,
    4.3f,
    5.9f,
    6.8f
};

const PointsList2D points2D =
{
    { 0.0f, 1.1f },
    { 1.6f, 8.3f },
    { 2.3f, 6.5f },
    { 3.5f, 4.7f },
    { 4.3f, 3.1f },
    { 5.9f, 7.5f },
    { 6.8f, 0.0f }
};

float mean(std::vector<float> x)
{
    float res = 0.0;
    for(int i = 0; i < x.size(); i++)
    {
        res += x[i];
    }
    res /= x.size();
    return res;
}

float rmse(std::vector<float> x, std::vector<float> y)
{	
	float res=0;
	int size = std::min(x.size(), y.size());
	for(int i=0; i<size; i++) {
		res+=pow((x[i]-y[i]), 2.0);
	}
	return sqrt(res/(float (size)));
}

std::pair<float, float> linearRegression(std::vector<float> x, std::vector<float> y)
{
    float a = 0.0;
	float b = 0.0;
    float numerator = 0.0;
    float denominator = 0.0;
    float xMean = mean(x);
    float yMean = mean(y);
	int size = std::min(x.size(), y.size());

    for(int i = 0; i < size; i++)
    {
        numerator += (x[i] - xMean) * (y[i] - yMean);
        denominator += (x[i] - xMean) * (x[i] - xMean);
    }
    a = numerator / denominator;
    b = yMean - a * xMean;
    return std::make_pair(a, b);
}

int main() {
	auto interpolationCubic = CubicInterpolation(points2D, points1D);
		interpolationCubic.Interpolate2D(13);
		std::cout << std::endl;
		interpolationCubic.Interpolate1D(13);
		std::cout << std::endl;

	auto interpolationLagrange = LagrangeInterpolation(points1D);
		interpolationLagrange.Interpolate1D(13);
		std::cout << std::endl;

	std::vector<float> x{ 1.0, 2.0, 4.0, 3.0, 5.0 }; 
	std::vector<float> y{ 1.0, 3.0, 3.0, 2.0, 5.0 };

	std::pair<float, float> result = linearRegression(x,y);

	std::vector<float> func_results(x.size());
	for(int i=0; i<x.size(); i++) {
		func_results[i]=(x[i]*result.first + result.second);
	} 

	std::cout<< "Result - y=" << result.first << "x+" << result.second << std::endl;
	std::cout<< "Error=" << rmse(func_results, y) << std::endl;

	system ("PAUSE");
}