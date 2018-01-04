#include "main.h"

void CalculateStep(cplx_t* wf, const flt_t* v, const bool* mask, const flt_t dt, const flt_t dx, const flt_t& bec_as);

cplx_t* k1 = new cplx_t[w * h]; //runge kutta steps
cplx_t* k2 = new cplx_t[w * h];
cplx_t* k3 = new cplx_t[w * h];
cplx_t* k4 = new cplx_t[w * h];
cplx_t* kSub = new cplx_t[w * h]; //runge kutta substep

const flt_t potentialScaling = 10;

int main()
{
	cplx_t* wf = new cplx_t[w * h]; //wave function

	flt_t* v = new flt_t[w * h]; //potential
	bool* mask = new bool[w * h]; //dirichlet condition mask. false=infinity

	constexpr flt_t dx = 0.1;
	constexpr flt_t dt = dx * dx / 2.0;

	flt_t bec_as = 0;

	for (int y = 0; y < h; ++y)
		for (int x = 0; x < w; ++x)
		{
			wf[y * w + x] = cplx_t(static_cast<flt_t>(0.0), static_cast<flt_t>(0.0));
			k1[y * w + x] = cplx_t();
			k2[y * w + x] = cplx_t();
			k3[y * w + x] = cplx_t();
			k4[y * w + x] = cplx_t();
			kSub[y * w + x] = cplx_t();
		}
	if (!LoadWFScenario(wf, &bec_as))
	{
		std::cout << "Error while reading wf scenario file" << std::endl;
		std::system("pause");
		return 0;
	}

	std::vector<uint8_t> imBuffer, image; //for image loadings
	const char* maskFile = "mask.png";
	const char* potentialFile = "potential.png";
	unsigned long pngW, pngH;

	//load dirichlet image
	loadFile(imBuffer, maskFile);
	int error = decodePNG(image, pngW, pngH, imBuffer.empty() ? 0 : &imBuffer[0], (unsigned long)imBuffer.size());
	if(error != 0 || pngW != w || pngH != h)
	{
		std::cout << "error: " << error << std::endl;
		std::system("pause");
		return error;
	}
	for (int i = 0; i < w * h * 4; i += 4)
	{
		mask[i / 4] = (static_cast<unsigned int>(image[i])
			+ static_cast<unsigned int>(image[i + 1])
			+ static_cast<unsigned int>(image[i + 2])) == 0;
		if (mask[i / 4])
			wf[i / 4] = 0;
	}
	
	//load potential image
	loadFile(imBuffer, potentialFile);
	error = decodePNG(image, pngW, pngH, imBuffer.empty() ? 0 : &imBuffer[0], (unsigned long)imBuffer.size());
	if(error != 0 || pngW != w || pngH != h)
	{
		std::cout << "error: " << error << std::endl;
		std::system("pause");
		return error;
	}
	for (int i = 0; i < w * h * 4; i += 4)
		v[i / 4] = (static_cast<flt_t>(static_cast<unsigned int>(image[i])
			+ static_cast<unsigned int>(image[i + 1])
			+ static_cast<unsigned int>(image[i + 2])) / static_cast<flt_t>(255.0 * 3.0)) * potentialScaling;

	sf::RenderWindow window(sf::VideoMode(w, h), "Qm Sim");

	sf::Texture texture;
	texture.create(w, h);
	texture.setSmooth(false);
	texture.setRepeated(false);
	uint8_t* buffer = new uint8_t[w * h * 4];
	for (int y = 0; y < h; ++y)
		for (int x = 3; x < w * 4; x += 4)
			buffer[y * w * 4 + x] = 255;

	sf::Sprite sprite;

	CalculateStep(wf, v, mask, dt, dx, bec_as);
	bool inCalculation = false;
	while (window.isOpen())
    {
        sf::Event event;
		bool windowWasClosed = false;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
			{
                window.close();
				windowWasClosed = true;
				break;
			}
			if (event.type == sf::Event::KeyReleased)
			{
				if (event.key.code == sf::Keyboard::Escape)
				{
					window.close();
					windowWasClosed = true;
					break;
				}
				inCalculation = !inCalculation;
			}
		}
		if (windowWasClosed)
			break;

		if (inCalculation)
			for (int st = 0; st < imidiateSteps; ++st)
				CalculateStep(wf, v, mask, dt, dx, bec_as);

		flt_t maxReFactor = static_cast<flt_t>(0.0), maxImFactor = static_cast<flt_t>(0.0), maxNormFactor = static_cast<flt_t>(0.0);
		for (int y = 0; y < h; ++y)
		{
			int offset = y * w;
			for (int x = 0; x < w; ++x)
			{
				if (std::abs(wf[offset + x].real()) > maxReFactor)
					maxReFactor = std::abs(wf[offset + x].real());
				if (std::abs(wf[offset + x].imag()) > maxImFactor)
					maxImFactor = std::abs(wf[offset + x].imag());
				if (std::norm(wf[offset + x]) > maxNormFactor)
					maxNormFactor = std::norm(wf[offset + x]);
			}
		}
		maxReFactor = static_cast<flt_t>(255.0) / maxReFactor;
		maxImFactor = static_cast<flt_t>(255.0) / maxImFactor;
		maxNormFactor = static_cast<flt_t>(255.0) / maxNormFactor;
		#pragma omp parallel for
		for (int y = 0; y < h; ++y)
		{
			for (int x = 0; x < w; ++x)
			{
				buffer[(y * w + x) * 4    ] = static_cast<uint8_t>(std::abs(wf[y * w + x].real()) * maxReFactor);
				buffer[(y * w + x) * 4 + 1] = static_cast<uint8_t>(std::norm(wf[y * w + x]) * maxNormFactor);
				buffer[(y * w + x) * 4 + 2] = static_cast<uint8_t>(std::abs(wf[y * w + x].imag()) * maxImFactor);
			}
		}

		texture.update(buffer);
		sprite.setTexture(texture);
		window.draw(sprite);
		
		window.display();
	}

	delete[] buffer;
	delete[] wf;
	delete[] k1;
	delete[] k2;
	delete[] k3;
	delete[] k4;
	delete[] kSub;
	delete[] v;
	delete[] mask;
}

void inline KStep(cplx_t* dst, const cplx_t* src, const flt_t* v, const bool* mask, const flt_t dx, const flt_t& as)
{
	flt_t spaceDerivativeQuotient = static_cast<flt_t>(1.0) / (dx * dx * static_cast<flt_t>(2.0));
	#pragma omp parallel for
	for (int y = 1; y < h - 1; ++y)
	{
		for (int x = 1; x < w - 1; ++x)
		{
			dst[y * w + x] = li * (
				((src[y * w + x + 1] + src[y * w + x - 1] - static_cast<flt_t>(2.0) * src[y * w + x])
				+ (src[(y + 1) * w + x] + src[(y - 1) * w + x] - static_cast<flt_t>(2.0) * src[y * w + x])) * spaceDerivativeQuotient
				- v[y * w + x] * src[y * w + x]
				- std::norm(src[y * w + x]) * src[y * w + x] * static_cast<flt_t>(2.0 * M_PI * as)
				);
		}
	}
	#pragma omp parallel for
	for (int y = 1; y < h - 1; ++y)
		for (int x = 1; x < w - 1; ++x)
			if (mask[y * w + x])
				dst[y * w + x] = 0;
}

void inline CalculateStep(cplx_t* wf, const flt_t* v, const bool* mask, const flt_t dt, const flt_t dx, const flt_t& bec_as)
{
	KStep(k1, wf, v, mask, dx, bec_as);
	#pragma omp parallel for
	for (int y = 0; y < h; ++y)
		for (int x = 0; x < w; ++x)
			kSub[y * w + x] = (dt / static_cast<flt_t>(2.0)) * k1[y * w + x] + wf[y * w + x];
	KStep(k2, kSub, v, mask, dx, bec_as);
	#pragma omp parallel for
	for (int y = 0; y < h; ++y)
		for (int x = 0; x < w; ++x)
			kSub[y * w + x] = (dt / static_cast<flt_t>(2.0)) * k2[y * w + x] + wf[y * w + x];
	KStep(k3, kSub, v, mask, dx, bec_as);
	#pragma omp parallel for
	for (int y = 0; y < h; ++y)
		for (int x = 0; x < w; ++x)
			kSub[y * w + x] = dt * k3[y * w + x] + wf[y * w + x];
	KStep(k4, kSub, v, mask, dx, bec_as);
	#pragma omp parallel for
	for (int y = 0; y < h; ++y)
		for (int x = 0; x < w; ++x)
			wf[y * w + x] += (dt / static_cast<flt_t>(6.0))
				* (k1[y * w + x]
				+ static_cast<flt_t>(2.0) * k2[y * w + x]
				+ static_cast<flt_t>(2.0) * k3[y * w + x]
				+ k4[y * w + x]);
}