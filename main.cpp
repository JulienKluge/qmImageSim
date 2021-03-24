#include "main.h"
#define M_Pi 3.14159265358f

void CalculateStep(cplx_t* wf, const flt_t* v, const bool* mask, const flt_t dt, const flt_t dx, const flt_t& bec_as);

cplx_t* k1 = new cplx_t[w * h]; //runge kutta steps
cplx_t* k2 = new cplx_t[w * h];
cplx_t* k3 = new cplx_t[w * h];
cplx_t* k4 = new cplx_t[w * h];
cplx_t* kSub = new cplx_t[w * h]; //runge kutta substep

const flt_t potentialScaling = 100;

int main()
{
	bool render_frames = true;
	int frame = 0;
	int old_frame = -1;
	sf::RenderWindow window(sf::VideoMode(w, h), "Qm Sim");
reload: //sue me...
	std::ios_base::sync_with_stdio(false);

	cplx_t* wf = new cplx_t[w * h]; //wave function

	flt_t* v = new flt_t[w * h]; //potential
	bool* mask = new bool[w * h]; //dirichlet condition mask. false=infinity

	constexpr flt_t dx = 0.1;
	constexpr flt_t dt = dx * dx / 2.0;

	flt_t bec_as = 1.0;

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
	
	flt_t* p_x;
	flt_t* p_y;
	int p_num = 0;
	flt_t p_m = 40.0;
	if (!LoadWFScenario(wf, &p_x, &p_y, &p_num, &p_m, &bec_as))
	{
		std::cout << "Error while reading wf scenario file" << std::endl;
		std::system("pause");
		return 0;
	}

	std::vector<uint8_t> imBuffer, image; //for image loadings
	std::string maskFile = "mask.png";
	std::string potentialFile = "potential.png";
	unsigned long pngW, pngH;

	//load dirichlet image

	//loadFile2(imBuffer, maskFile);
	imBuffer = load_file(maskFile);
	
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
	//loadFile2(imBuffer, potentialFile);
	imBuffer = load_file(potentialFile);
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
				else if (event.key.code == sf::Keyboard::Delete)
				{
					inCalculation = false;
					goto reload; //I SAID SUE ME!!
				}
				inCalculation = !inCalculation;
			}
		}
		if (windowWasClosed)
			break;

		if (inCalculation)
		{
			for (int st = 0; st < imidiateSteps; ++st)
			{
				CalculateStep(wf, v, mask, dt, dx, bec_as);
				for (int p = 0; p < p_num; ++p)
				{
					flt_t temp_px = p_x[p];
					flt_t temp_py = p_y[p];
					Calculate_Particle_Step(wf, &temp_px, &temp_py, dt, dx, p_m);
					p_x[p] = temp_px;
					p_y[p] = temp_py;
				}
			}
			++frame;
			/*int mX = 1, mY = 1;
			FindModes(wf, mX, mY);*/
		}
		

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
		for (int p = 0; p < p_num; ++p)
		{
			const int x = static_cast<int>((w - 1) * p_x[p]);
			const int y = static_cast<int>((h - 1) * p_y[p]);
			buffer[(y * w + x) * 4    ] = 255;
			buffer[(y * w + x) * 4 + 1] = 255;
			buffer[(y * w + x) * 4 + 2] = 255;
		}

		texture.update(buffer);
		sprite.setTexture(texture);
		window.draw(sprite);
		
		window.display();
		
		if (render_frames && inCalculation && frame != old_frame)
		{
			old_frame = frame;
			SaveBuffer(buffer);
		}
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

void inline Calculate_Particle_Step(cplx_t* wf, flt_t* x, flt_t* y, const flt_t& dt, const flt_t& dx, const flt_t& m)
{
	const int x_idx = static_cast<int>((w - 1) * (*x));
	const int y_idx = static_cast<int>((h - 1) * (*y));

	if (isnan(*x) || isnan(*y))
		return;
	//std::cout << (*x) << "  " << (*y) << std::endl;
	
	cplx_t der;
	flt_t acc;
	if (x_idx == 0)
		der = -li * (wf[y_idx * w + x_idx + 1] - wf[y_idx * w + x_idx]) / cplx_t{dx, 0};
	else if (x_idx == (w - 1))
		der = -li * (wf[y_idx * w + x_idx] - wf[y_idx * w + x_idx - 1]) / cplx_t{dx, 0};
	else
		der = -li * (wf[y_idx * w + x_idx + 1] - wf[y_idx * w + x_idx - 1]) / cplx_t{2 * dx, 0};
	acc = (der / wf[y_idx * w + x_idx]).real() * dt / m;
	//std::cout << "acc " << acc << "  isnan = " << (acc != acc) << std::endl;
	if (!isnan(acc))
		(*x) = (*x) + acc;

	if (y_idx == 0)
		der = -li * (wf[(y_idx + 1) * w + x_idx] - wf[y_idx * w + x_idx]) / cplx_t{dx, 0};
	else if (y_idx == (h - 1))
		der = -li * (wf[y_idx * w + x_idx] - wf[(y_idx - 1) * w + x_idx]) / cplx_t{dx, 0};
	else
		der = -li * (wf[(y_idx + 1) * w + x_idx] - wf[(y_idx - 1) * w + x_idx]) / cplx_t{2 * dx, 0};
	acc = (der / wf[y_idx * w + x_idx]).real() * dt / m;
	if (!isnan(acc))
		(*y) += acc;


	if ((*x) > 0.9999)
		(*x) = 0.9999;
	else if ((*x) < 0.0001)
		(*x) = 0.0001;
	
	if ((*y) > 0.9999)
		(*y) = 0.9999;
	else if ((*y) < 0.0001)
		(*y) = 0.0001;
	
	//std::cout << (*x) << "  " << (*y) << std::endl;
}

void SaveBuffer(const uint8_t* frame_buffer)
{
	static std::ofstream frameFile;
	frameFile.open("frameFile.bin", std::ios_base::app | std::ios_base::binary);
	frameFile.write((char*)frame_buffer, w * h * 4);
	//for (int i = 0; i < w * h * 4; ++i)
	//	frameFile << frame_buffer[i];
	frameFile.flush();
	frameFile.close();
}

void inline FindModes(const cplx_t* wf, int& modesX, int& modesY)
{
	static flt_t* xVals = new flt_t[w];
	static flt_t* yVals = new flt_t[h];
	static int lastXMode = -1;
	static int lastYMode = -1;
	static int lastUpdatedStep = 0;
	static std::ofstream modeFile;

	if (!modeFile.is_open())
		modeFile.open("modefile.csv", std::ios_base::app);
	

	flt_t xMax = 0;
	flt_t yMax = 0;
	for (int x = 0; x < w; ++x)
	{
		xVals[x] = 0;
		for(int y = 0; y < h; ++y)
			xVals[x] += std::abs(wf[y * w + x]);
		if (xVals[x] > xMax)
			xMax = xVals[x];
	}
	for (int y = 0; y < h; ++y)
	{
		yVals[y] = 0;
		for(int x = 0; x < w; ++x)
			yVals[y] += std::abs(wf[y * w + x]);
		if (yVals[y] > yMax)
			yMax = yVals[y];
	}
	const flt_t xUpperLimit = 0.9 * xMax;
	const flt_t xLowerLimit = 0.7 * xMax;

	const flt_t yUpperLimit = 0.9 * yMax;
	const flt_t yLowerLimit = 0.7 * yMax;
	
	modesX = 0;
	modesY = 0;
	bool gotXStreak = false, gotYStreak = false;
	for (int i = 0; i < w; ++i)
	{
		if (gotXStreak)
		{
			if (xVals[i] < xLowerLimit)
				gotXStreak = false;
		}
		else
		{
			if (xVals[i] > xUpperLimit)
			{
				gotXStreak = true;
				modesX += 1;
			}
		}
	}
	for (int i = 0; i < h; ++i)
	{
		if (gotYStreak)
		{
			if (yVals[i] < yLowerLimit)
				gotYStreak = false;
		}
		else
		{
			if (yVals[i] > yUpperLimit)
			{
				gotYStreak = true;
				++modesY;
			}
		}
	}

	if (lastXMode != modesX || lastYMode != modesY)
	{
		modeFile << modesX << ',' << modesY << ',' << lastUpdatedStep << '\n';
		lastXMode = modesX;
		lastYMode = modesY;
		lastUpdatedStep = 0;
	}
	else
	{
		++lastUpdatedStep;
	}
}