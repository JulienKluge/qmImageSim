#include <complex>
#include <iostream>
#include <fstream>
#include <math.h>

#include "picopng.cpp"
#include "SFML/Graphics.hpp"

using uint8_t = unsigned char;
using flt_t = float;
using cplx_t = std::complex<flt_t>;

#define POW2(x) ((x) * (x))

const int imidiateSteps = 1;
const int w = 500, h = 300;
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

bool LoadWFScenario(cplx_t* wf, flt_t* bec_as)
{
	std::ifstream scenario("scenario.txt");
	std::string line, identifierStr;
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
			loadFile(imBuffer, fileName.c_str());
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
	}
	NormalizeWF(wf);
	return true;
}
