#include <complex>
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>

using namespace std::complex_literals;

#include "picopng.cpp"
#include "SFML/Graphics.hpp"

using uint8_t = unsigned char;
using flt_t = float;
using cplx_t = std::complex<flt_t>;

#define POW2(x) ((x) * (x))

const int imidiateSteps = 1;
const int w = 1000, h = 600;
const cplx_t li(0.0f, 1.0f); //why am i unable to find this literal?


void NormalizeWF(cplx_t* wf)
{
	double sum = 0.0;
	for (int y = 0; y < h; ++y)
		for (int x = 0; x < w; ++x)
			sum += std::norm(wf[y * w + x]);
	for (int y = 0; y < h; ++y)
		for (int x = 0; x < w; ++x)
			wf[y * w + x] /= static_cast<flt_t>(sum);
}

bool LoadWFScenario(cplx_t* wf, flt_t** p_x, flt_t** p_y, int* p_num, flt_t* p_m, flt_t* bec_as)
{
	std::ifstream scenario("scenario.txt");
	std::string line, identifierStr;
	std::vector<flt_t> l_p_x;
	std::vector<flt_t> l_p_y;
	int p_n;
	while (std::getline(scenario, line))
	{
		if (line[0] == '%' || line.length() < 3)
			continue;
		std::stringstream ss(line);
		ss >> identifierStr;
		if (identifierStr == "bec_coupling_constant")
		{
			flt_t subBECValue;
			ss >> subBECValue;
			*bec_as = subBECValue;
		}
		else if (identifierStr == "gauss")
		{
			flt_t posX, posY, sigmaX, sigmaY, impulseX, impulseY;
			ss >> posX;
			ss >> posY;
			ss >> sigmaX;
			ss >> sigmaY;
			ss >> impulseX;
			ss >> impulseY;
			std::cout << "load gauss with x:" << posX << " y:" << posY << '\n';
			std::cout << "sx: " << sigmaX << " sy:" << sigmaY << " px:" << impulseX << " py:" << impulseY << std::endl;
			for (int y = 0; y < h; ++y)
			{
				flt_t yy = static_cast<flt_t>(y);
				for (int x = 0; x < w; ++x)
				{
					flt_t xx = static_cast<flt_t>(x);
					wf[y * w + x] += std::exp(
						li * (impulseX * xx + impulseY * yy) //impulse
						-((xx - posX) * (xx - posX) / sigmaX + (yy - posY) * (yy - posY) / sigmaY) //amplitude
					);
				}
			}
		}
		else if (identifierStr == "plainwave_horizontal")
		{
			flt_t posX, posY, width, height, impulse;
			ss >> posX;
			ss >> posY;
			ss >> width;
			ss >> height;
			ss >> impulse;
			std::cout << "load horizontal cos wave with x:" << posX << " y:" << posY << '\n';
			std::cout << "w: " << width << " h:" << height << " p:" << impulse << std::endl;
			flt_t height2 = static_cast<flt_t>(height) / 2.0;
			for (int y = posY; y < (posY + height); ++y)
			{
				flt_t realY = (y - posY) / height - height2;
				flt_t amplitude = (1 - std::cos(realY * M_PI * 2.0)) / 2.0;
				cplx_t value = amplitude * std::exp(li * (impulse * y));
				for (int x = posX; x < (posX + width); ++x)
					wf[y * w + x] += value;
			}
		}
		else if (identifierStr == "plainwave_vertical")
		{
			flt_t posX, posY, width, height, impulse;
			ss >> posX;
			ss >> posY;
			ss >> width;
			ss >> height;
			ss >> impulse;
			std::cout << "load vertical cos wave with x:" << posX << " y:" << posY << '\n';
			std::cout << "w: " << width << " h:" << height << " p:" << impulse << std::endl;
			flt_t width2 = static_cast<flt_t>(width) / 2.0;
			for (int x = posX; x < (posX + width); ++x) //arghhh cache coherence horror
			{
				flt_t realX = (x - posX) / width - width2;
				flt_t amplitude = (1 - std::cos(realX * M_PI * 2.0)) / 2.0;
				cplx_t value = amplitude * std::exp(li * (impulse * x));
				for (int y = posY; y < (posY + height); ++y)
					wf[y * w + x] += value;
			}
		}
		else if (identifierStr == "constant")
		{
			flt_t posX, posY, width, height;
			ss >> posX;
			ss >> posY;
			ss >> width;
			ss >> height;
			std::cout << "load constant x:" << posX << " y:" << posY << '\n';
			std::cout << "w: " << width << " h:" << height << std::endl;
			for (int y = posY; y < (posY + height); ++y)
			{
				for (int x = posX; x < (posX + width); ++x)
					wf[y * w + x] += 1;
			}
		}
		else if (identifierStr == "grayscaleblit")
		{
			flt_t posX, posY;
			ss >> posX;
			ss >> posY;
			std::string fileName;
			ss >> fileName;
			std::cout << "load grayscale image x:" << posX << " y:" << posY << '\n';
			std::cout << fileName << std::endl;
			std::vector<uint8_t> imBuffer, image; //for image loadings
			unsigned long pngW, pngH;
			loadFile2(imBuffer, fileName);
			int error = decodePNG(image, pngW, pngH, imBuffer.empty() ? 0 : &imBuffer[0], (unsigned long)imBuffer.size());
			if(error != 0)
				std::cout << "error while decoding png" << std::endl;
			else
			{
				for (int y = posY; y < h; ++y)
				{
					int realY = (y - posY) * 4 * pngW;
					for (int x = posX; x < w; ++x)
					{
						flt_t sum = static_cast<flt_t>(image[realY + (x - posX) * 4])
							+ static_cast<flt_t>(image[realY + (x - posX) * 4 + 1])
							+ static_cast<flt_t>(image[realY + (x - posX) * 4 + 2]);
						sum /= (3 * 255);
						wf[y * w + x] += sum;
					}
				}
			}
		}
		else if (identifierStr == "particle")
		{
			flt_t posX, posY;
			ss >> posX;
			ss >> posY;
			l_p_x.push_back(posX);
			l_p_y.push_back(posY);
			++(*p_num);
		}
		else if (identifierStr == "particlefield")
		{
			flt_t posX, posY;
			flt_t width, height;
			int num_x, num_y;
			ss >> posX;
			ss >> posY;
			ss >> width;
			ss >> height;
			ss >> num_x;
			ss >> num_y;

			for (int i = 0; i < num_y; ++i)
			{
				flt_t p_yPos = (static_cast<flt_t>(i) / static_cast<flt_t>(num_y - 1) * height + posY);
				for (int j = 0; j < num_x; ++j)
				{
					flt_t p_xPos = (static_cast<flt_t>(j) / static_cast<flt_t>(num_x - 1) * width + posX);
					l_p_x.push_back(p_xPos);
					l_p_y.push_back(p_yPos);
				}
			}

			(*p_num) += num_x * num_y;
		}
		else if (identifierStr == "particlemass")
		{
			flt_t p_mass;
			ss >> p_mass;
			(*p_m) = p_mass;
		}
	}
	(*p_x) = new flt_t[(*p_num)];
	(*p_y) = new flt_t[(*p_num)];
	for (int p = 0; p < (*p_num); ++p)
	{
		(*p_x)[p] = l_p_x[p] / static_cast<flt_t>(w);
		(*p_y)[p] = l_p_y[p] / static_cast<flt_t>(h);
	}
	NormalizeWF(wf);
	return true;
}

void Calculate_Particle_Step(cplx_t* wf, flt_t* x, flt_t* y, const flt_t& dt, const flt_t& dx, const flt_t& m);
void SaveBuffer(const uint8_t* wf);
void FindModes(const cplx_t* wf, int& modesX, int& modesY);
