// Homework.cpp : Defines the entry point for the console application.
// Canny edge detector and template matching
// Author: Aditya Bettadapura

#include <afxwin.h>  // necessary for MFC to work properly
#include "Homework.h"
#include "../../src/blepo.h"
#include "math.h"
#include <stack>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#define PI 3.14

using namespace blepo;

void GetChamferDistance(ImgBinary hysteresisImage, ImgInt *chamferImage)
{
	int maxvalue = hysteresisImage.Width() * hysteresisImage.Height() + 1;

	//first pass
	for (int y = 0; y < hysteresisImage.Height(); y++)
	{
		for (int x = 0; x < hysteresisImage.Width(); x++)
		{
			if (hysteresisImage(x, y)) (*chamferImage)(x, y) = 0;
			else
			{
				int distance = maxvalue;
				if (y > 0)  distance = min(distance, (*chamferImage)(x, y - 1) + 1);
				if (x > 0)  distance = min(distance, (*chamferImage)(x - 1, y) + 1);
				(*chamferImage)(x, y) = distance;
			}
		}
	}

	//second pass
	for (int y = hysteresisImage.Height() - 1; y >= 0; y--)
	{
		for (int x = hysteresisImage.Width() - 1; x >= 0; x--)
		{
			if (hysteresisImage(x, y))  (*chamferImage)(x, y) = 0;
			else
			{
				int distance = (*chamferImage)(x, y);
				if (y < hysteresisImage.Height() - 1)  distance = min(distance, (*chamferImage)(x, y + 1) + 1);
				if (x < hysteresisImage.Width() - 1)  distance = min(distance, (*chamferImage)(x + 1, y) + 1);
				(*chamferImage)(x, y) = distance;
			}
		}
	}
}


void GetGradient(ImgGray img, int width, int halfwidth, ImgFloat *outx, ImgFloat *outy, double dervkern[], double kern[])
{
	double sum = 0;
	int x, y, i, j;
	double *temp;
	temp = (double *)calloc(img.Height()*img.Width(), sizeof(double));
	for (y = 0; y < img.Height(); y++)
	{
		for (x = halfwidth; x < img.Width() - halfwidth - 1; x++)
		{
			sum = 0;
			for (i = 0; i < width; i++)
			{
				sum += ((img(x + halfwidth - i, y))*(dervkern[i]));
			}
			temp[y*img.Width() + x] = sum;
		}
	}
	for (y = halfwidth; y < img.Height() - halfwidth - 1; y++)
	{
		for (x = 0; x < img.Width(); x++)
		{
			sum = 0;
			for (i = 0; i < width; i++)
			{
				sum += ((temp[(y + halfwidth - i)* img.Width() + x])*(kern[i]));
			}
			(*outx)(x, y) = sum;
		}
	}

	for (y = 0; y < img.Height(); y++)
	{
		for (x = halfwidth; x < img.Width() - halfwidth - 1; x++)
		{
			sum = 0;
			for (i = 0; i < width; i++)
			{
				sum += ((img(x + halfwidth - i, y))*(kern[i]));
			}
			temp[y*img.Width() + x] = sum;
		}
	}
	for (y = halfwidth; y < img.Height() - halfwidth - 1; y++)
	{
		for (x = 0; x < img.Width(); x++)
		{
			sum = 0;
			for (i = 0; i < width; i++)
			{
				sum += ((temp[(y + halfwidth - i)*img.Width() + x])*(dervkern[i]));
			}
			(*outy)(x, y) = sum;
		}
	}

}

void DoubleThresholding(ImgGray *inputImage, ImgBinary *highImage, ImgBinary *lowImage, ImgBinary *hysteresisImage, ImgFloat *outputGrad, int highThreshold, int lowThreshold)
{

	for (int y = 1; y < ((*inputImage).Height()); y++)
	{
		for (int x = 1; x < (*inputImage).Width(); x++)
		{
			if ((*outputGrad)(x, y) >= highThreshold)
			{
				(*highImage)(x, y) = 1;
				(*hysteresisImage)(x, y) = 1;
			}
			else if (((*outputGrad)(x, y) > lowThreshold) && ((*outputGrad)(x, y) < highThreshold))
			{
				(*lowImage)(x, y) = 1;
			}
		}
	}

	for (int y = 1; y < ((*inputImage).Height()-1); y++)
	{
		for (int x = 1; x < (*inputImage).Width()-1; x++)
		{
			if ((*lowImage)(x, y) == 1)
			{
				
					if (((*highImage)(x - 1, y - 1) == 1) ||
					((*highImage)(x, y - 1) == 1) ||
					((*highImage)(x + 1, y - 1) == 1) ||
					((*highImage)(x - 1, y) == 1) ||
					((*highImage)(x + 1, y) == 1) ||
					((*highImage)(x - 1, y + 1) == 1) ||
					((*highImage)(x, y + 1) == 1) ||
					((*highImage)(x + 1, y + 1) == 1))
				{
					(*hysteresisImage)(x, y) = 1;
				}
			}
		}
	}

	for (int y = 1; y < ((*inputImage).Height()-1); y++)
	{
		for (int x = 1; x < (*inputImage).Width()-1; x++)
		{
			if ((*hysteresisImage)(x, y) == 1)
			{
				if (((*highImage)(x - 1, y - 1) == 1) ||
					((*highImage)(x, y - 1) == 1) ||
					((*highImage)(x + 1, y - 1) == 1) ||
					((*highImage)(x - 1, y) == 1) ||
					((*highImage)(x + 1, y) == 1) ||
					((*highImage)(x - 1, y + 1) == 1) ||
					((*highImage)(x, y + 1) == 1) ||
					((*highImage)(x + 1, y + 1) == 1))
				{
					(*hysteresisImage)(x, y) = 1;
				}
			}
		}
	}

	for (int y = 1; y < ((*inputImage).Height()-1); y++)
	{
		for (int x = 1; x < (*inputImage).Width()-1; x++)
		{
			if (((*hysteresisImage)(x, y) == 0) && ((*lowImage)(x, y) == 1))
			{
				if (((*hysteresisImage)(x - 1, y - 1) == 1) ||
					((*hysteresisImage)(x, y - 1) == 1) ||
					((*hysteresisImage)(x + 1, y - 1) == 1) ||
					((*hysteresisImage)(x - 1, y) == 1) ||
					((*hysteresisImage)(x + 1, y) == 1) ||
					((*hysteresisImage)(x - 1, y + 1) == 1) ||
					((*hysteresisImage)(x, y + 1) == 1) ||
					((*hysteresisImage)(x + 1, y + 1) == 1))
				{
					(*hysteresisImage)(x, y) = 1;
				}
			}
		}
	}
}

int main(int argc, const char* argv[], const char* envp[])
{
	// Initialize MFC and return if failure
	HMODULE hModule = ::GetModuleHandle(NULL);
	if (hModule == NULL || !AfxWinInit(hModule, NULL, ::GetCommandLine(), 0))
	{
		printf("Fatal Error: MFC initialization failed (hModule = %x)\n", hModule);
		return 1;
	}


	try
	{
		//Read sigma value
		float sigma = atof(argv[1]);

		//Read first image name
		string argument2 = argv[2];
		string path = "../../images/";
		string filename1 = path + argument2;
		const char *file1 = filename1.c_str();

		ImgGray inputImage;
		Load(file1, &inputImage);

		ImgBgr finalImage;
		Load(file1, &finalImage);

		Figure fig1;
		fig1.SetTitle("Original Image");
		fig1.Draw(inputImage);

		//If third argument exists

		ImgGray templateImage;
		ImgBinary temphysteresisImage;


		//calculate width
		int mu = (int)(2.5*sigma - 0.5);
		int w = 2 * mu + 1;
		double *gauss;
		double *gaussDeriv;
		double sumGauss = 0;
		double sumGaussDeriv = 0;

		gauss = (double*)calloc(w, sizeof(double));
		gaussDeriv = (double*)calloc(w, sizeof(double));

		printf("width = %d\n", w);
		printf("mu = %d\n", mu);

		//create gaussian kernel
		for (int i = 0; i < w; i++)
		{
			gauss[i] = exp((-(i - mu)*(i - mu)) / (2 * sigma*sigma));
			sumGauss += gauss[i];

			gaussDeriv[i] = (i - mu)*gauss[i];
			sumGaussDeriv += gaussDeriv[i] * i;
		}

		for (int j = 0; j < w; j++)
		{
			gauss[j] = gauss[j] / sumGauss;
			gaussDeriv[j] = gaussDeriv[j] / sumGaussDeriv;
		}

		//Print Gaussian Kernel
		printf("Gaussian Kernel = [");
		for (int i = 0; i < w; i++)
		{
			printf("%f ", gauss[i]);
		}
		printf("]\n");

		//Print Gaussian Derivative kernel
		printf("Gaussian Derivative = [");
		for (int i = 0; i < w; i++)
		{
			printf("%f ", gaussDeriv[i]);
		}
		printf("]\n");

		//Compute gradient of the image
		ImgFloat grad_x, grad_y;

		grad_x.Reset(inputImage.Width(), inputImage.Height());
		grad_y.Reset(inputImage.Width(), inputImage.Height());

		Set(&grad_x, 0);
		Set(&grad_y, 0);

		// Calculate Gradient
		GetGradient(inputImage, w, mu, &grad_x, &grad_y, gaussDeriv, gauss);

		Figure fig2, fig3;
		fig2.SetTitle("X-Gradient");
		fig2.Draw(grad_x);
		fig3.SetTitle("Y-Gradient");
		fig3.Draw(grad_y);

		//Calculate gradient magnitude and phase
		ImgFloat gradMagnitude, gradPhase, outputGrad;

		outputGrad.Reset(inputImage.Width(), inputImage.Height());
		gradMagnitude.Reset(inputImage.Width(), inputImage.Height());
		gradPhase.Reset(inputImage.Width(), inputImage.Height());

		Set(&gradMagnitude, 0);
		Set(&gradPhase, 0);
		Set(&outputGrad, 0);

		//Compute Gradient magnitude for input
		for (int y = mu; y < inputImage.Height() - mu; y++)
		{
			for (int x = mu; x < inputImage.Width() - mu; x++)
			{
				gradMagnitude(x, y) = max(abs(grad_x(x, y)), abs(grad_y(x, y)));
			}

		}

		//Compute Gradient phase for input
		for (int y = mu; y < inputImage.Height() - mu; y++)
		{
			for (int x = mu; x < inputImage.Width() - mu; x++)
			{
				float theta = atan2(grad_y(x, y), grad_x(x, y));
				gradPhase(x, y) = theta*180/PI;

				while (theta > 157.5)
				{
					theta -= 180;
				}

				if (theta >= -22.5 && theta < 22.5)
				{
					if (gradMagnitude(x, y) < gradMagnitude(x - 1, y) || (gradMagnitude(x, y) < gradMagnitude(x + 1, y)))
					{
						outputGrad(x, y) = 0;
					}
					else
					{
						outputGrad(x, y) = gradMagnitude(x, y);
					}
				}
				else if (theta >= 22.5 && theta < 67.5)
				{
					if (gradMagnitude(x, y) < gradMagnitude(x - 1, y - 1) || (gradMagnitude(x, y) < gradMagnitude(x + 1, y + 1)))
					{
						outputGrad(x, y) = 0;
					}
					else
					{
						outputGrad(x, y) = gradMagnitude(x, y);
					}
				}
				else if (theta >= 67.5 && theta < 112.5)
				{
					if (gradMagnitude(x, y) < gradMagnitude(x, y - 1) || (gradMagnitude(x, y) < gradMagnitude(x, y + 1)))
					{
						outputGrad(x, y) = 0;
					}
					else
					{
						outputGrad(x, y) = gradMagnitude(x, y);
					}
				}
				else if (theta >= 112.5 && theta < 157.5)
				{
					if (gradMagnitude(x, y) < gradMagnitude(x - 1, y + 1) || (gradMagnitude(x, y) < gradMagnitude(x + 1, y - 1)))
					{
						outputGrad(x, y) = 0;
					}
					else
					{
						outputGrad(x, y) = gradMagnitude(x, y) ;
					}
				}
			}
		}

		
		//Find threshold for input image
		int highThreshold = 0, lowThreshold = 0;

		double *temporaryArray;
		temporaryArray = (double*)calloc(inputImage.Height()*inputImage.Width(), sizeof(double));
		for (int y = 0; y < (inputImage.Height()); y++)
		{
			for (int x = 0; x < inputImage.Width(); x++)
			{
				temporaryArray[y*inputImage.Width() + x] = gradMagnitude(x, y);
			}
		}

		std::vector<double> sortedmagArray(temporaryArray, temporaryArray + (inputImage.Height()*inputImage.Width()));
		std::sort(sortedmagArray.begin(), sortedmagArray.end());

		highThreshold = sortedmagArray[0.9*(inputImage.Height()*inputImage.Width())];

		lowThreshold = highThreshold / 5;

		ImgBinary highImage, lowImage, hysteresisImage;
		highImage.Reset(inputImage.Width(), inputImage.Height());
		lowImage.Reset(inputImage.Width(), inputImage.Height());
		hysteresisImage.Reset(inputImage.Width(), inputImage.Height());

		Set(&highImage, 0);
		Set(&lowImage, 0);
		Set(&hysteresisImage, 0);


		DoubleThresholding(&inputImage, &highImage, &lowImage, &hysteresisImage, &outputGrad, highThreshold, lowThreshold);
		
		ImgInt chamferImage;
		chamferImage.Reset(inputImage.Width(), inputImage.Height());
		Set(&chamferImage, 100);
		GetChamferDistance(hysteresisImage, &chamferImage);

		//Edge template matching
		ImgFloat probabilityMap;
		probabilityMap.Reset(inputImage.Width(), inputImage.Height());
		Set(&probabilityMap, 255);



		if (argc == 4)
		{
			string argument3 = argv[3];
			string filename2 = path + argument3;
			const char *file2 = filename2.c_str();

			Load(file2, &templateImage);

			Figure fig8;
			fig8.SetTitle("Template Image");
			fig8.Draw(templateImage);

			ImgFloat tmplate_x, tmplate_y;
			tmplate_x.Reset(templateImage.Width(), templateImage.Height());
			tmplate_y.Reset(templateImage.Width(), templateImage.Height());

			Set(&tmplate_x, 0);
			Set(&tmplate_y, 0);

			GetGradient(templateImage, w, mu, &tmplate_x, &tmplate_y, gaussDeriv, gauss);

			ImgFloat templateMagnitude, templatePhase, templateOutputGrad;

			templateMagnitude.Reset(templateImage.Width(), templateImage.Height());
			templatePhase.Reset(templateImage.Width(), templateImage.Height());
			templateOutputGrad.Reset(templateImage.Width(), templateImage.Height());

			Set(&templateMagnitude, 0);
			Set(&templatePhase, 0);
			Set(&templateOutputGrad, 0);

			//Compute Gradient magnitude for template
			for (int y = mu; y < templateImage.Height() - mu; y++)
			{
				for (int x = mu; x < templateImage.Width() - mu; x++)
				{
					templateMagnitude(x, y) = max(abs(tmplate_x(x, y)), abs(tmplate_y(x, y)));
				}
			}

			//Compute Gradient phase for template
			for (int y = mu; y < templateImage.Height() - mu; y++)
			{
				for (int x = mu; x < templateImage.Width() - mu; x++)
				{
					float theta = atan2(tmplate_y(x, y), tmplate_x(x, y));
					templatePhase(x, y) = theta * 180 / PI;

					while (theta > 157.5)
					{
						theta -= 180;
					}

					if (theta >= -22.5 && theta < 22.5)
					{
						if (templateMagnitude(x, y) < templateMagnitude(x - 1, y) || (templateMagnitude(x, y) < templateMagnitude(x + 1, y)))
						{
							templateOutputGrad(x, y) = 0;
						}
						else
						{
							templateOutputGrad(x, y) = templateMagnitude(x, y);
						}
					}
					else if (theta >= 22.5 && theta < 67.5)
					{
						if (templateMagnitude(x, y) < templateMagnitude(x - 1, y - 1) || (templateMagnitude(x, y) < templateMagnitude(x + 1, y + 1)))
						{
							templateOutputGrad(x, y) = 0;
						}
						else
						{
							templateOutputGrad(x, y) = templateMagnitude(x, y);
						}
					}
					else if (theta >= 67.5 && theta < 112.5)
					{
						if (templateMagnitude(x, y) < templateMagnitude(x, y - 1) || (templateMagnitude(x, y) < templateMagnitude(x, y + 1)))
						{
							templateOutputGrad(x, y) = 0;
						}
						else
						{
							templateOutputGrad(x, y) = templateMagnitude(x, y);
						}
					}
					else if (theta >= 112.5 && theta < 157.5)
					{
						if (templateMagnitude(x, y) < templateMagnitude(x - 1, y + 1) || (templateMagnitude(x, y) < templateMagnitude(x + 1, y - 1)))
						{
							templateOutputGrad(x, y) = 0;
						}
						else
						{
							templateOutputGrad(x, y) = templateMagnitude(x, y);
						}
					}
				}
			}

			//Find threshold for template image
			int tempHighThreshold = 0, tempLowThreshold = 0;
			double *tempArray;
			tempArray = (double*)calloc(templateImage.Height()*templateImage.Width(), sizeof(double));
			for (int y = 0; y < (templateImage.Height()); y++)
			{
				for (int x = 0; x < templateImage.Width(); x++)
				{
					tempArray[y*templateImage.Width() + x] = templateMagnitude(x, y);
				}
			}

			std::vector<double> sortedtempArray(tempArray, tempArray + (templateImage.Height()*templateImage.Width()));
			std::sort(sortedtempArray.begin(), sortedtempArray.end());

			tempHighThreshold = sortedtempArray[0.9*(templateImage.Height()*templateImage.Width())];
			tempLowThreshold = tempHighThreshold / 5;

			ImgBinary temphighImage, templowImage;
			temphighImage.Reset(templateImage.Width(), templateImage.Height());
			templowImage.Reset(templateImage.Width(), templateImage.Height());
			temphysteresisImage.Reset(templateImage.Width(), templateImage.Height());

			Set(&temphighImage, 0);
			Set(&templowImage, 0);
			Set(&temphysteresisImage, 0);

			DoubleThresholding(&templateImage, &temphighImage, &templowImage, &temphysteresisImage, &templateOutputGrad, tempHighThreshold, tempLowThreshold);

		}


		//Create Inverse probability map
		for (int y = 0; y < (inputImage.Height() - templateImage.Height()); y++)
		{
			for (int x = 0; x < (inputImage.Width() - templateImage.Width()); x++)
			{
				int matchSum = 0;
				for (int j = 0; j < templateImage.Height(); j++)
				{
					for (int i = 0; i < templateImage.Width(); i++)
					{
						matchSum += temphysteresisImage(i, j)*hysteresisImage(x + i, y + j);
					}
				}
				probabilityMap(x, y) = 255-matchSum;
			}
		}

		//Find match in inverse probability map
		int probability = probabilityMap(0, 0);
		int col = 0, row = 0;

		for (int y = 0; y < inputImage.Height(); y++)
		{
			for (int x = 1; x < inputImage.Width(); x++)
			{
				if (probabilityMap(x, y) < probability)
				{
					probability = probabilityMap(x, y);
					col = x;
					row = y;
				}
			}
		}

		//Draw a rectangle around detected object
		Rect objectFound(col, row, col + templateImage.Width(), row + templateImage.Height());
		DrawRect(objectFound, &finalImage, Bgr(0,255,0));

		//Output Images
		Figure fig4, fig5, fig6, fig7, fig9;
		fig4.SetTitle("Magnitude of Gradient");
		fig4.Draw(gradMagnitude);
		fig5.SetTitle("Phase of Gradient");
		fig5.Draw(gradPhase);
		fig9.SetTitle("Non-maximal suppression image");
		fig9.Draw(outputGrad);
		fig6.SetTitle("Canny edge detection of input image");
		fig6.Draw(hysteresisImage);
		fig7.SetTitle("Chamfer distance image");
		fig7.Draw(chamferImage);

		if (argc == 4)
		{
			Figure  fig10, fig11, fig12;
			fig10.SetTitle("Canny edge detection of template image");
			fig10.Draw(temphysteresisImage);
			fig12.SetTitle("Inverse Probability Map");
			fig12.Draw(probabilityMap);
			fig11.SetTitle("Final image - Object found!");
			fig11.Draw(finalImage);
		}
		EventLoop();
	}
	catch (const Exception& e)
	{
		e.Display();    // display exception to user in a popup window 
	}
	return 0;
}
