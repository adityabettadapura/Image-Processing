// Homework.cpp : Defines the entry point for the console application.
// Marker and Non-Marker based Watersheding algorithm.
// Author: Aditya Bettadapura

#include <afxwin.h>  // necessary for MFC to work properly
#include "Homework.h"
#include "../../src/blepo.h"
#include "math.h"
#include <stack>
#include <queue>
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


void PerformFloodFill(ImgInt input, int x, int y, ImgInt *output, int new_label)
{
	stack<int> xpos;
	stack<int> ypos;
	int x_temp = x;
	int y_temp = y;
	int pixel = input(x, y);
	xpos.push(x_temp);
	ypos.push(y_temp);
	while (!xpos.empty() && !ypos.empty())
	{
		pixel = input(xpos.top(), ypos.top());
		x_temp = xpos.top();
		y_temp = ypos.top();
		xpos.pop();
		ypos.pop();
		for (int i = -1; i <= 1; i++)
		{
			for (int j = -1; j <= 1; j++)
			{
				if (x_temp + j < 0 || x_temp + j >= input.Width())
					continue;
				if (y_temp + i < 0 || y_temp + i >= input.Height())
					continue;
				if (input(x_temp + j, y_temp + i) == input(x_temp, y_temp) && (*output)(x_temp + j, y_temp + i) != new_label)
				{
					xpos.push(x_temp + j);
					ypos.push(y_temp + i);
					(*output)(x_temp + j, y_temp + i) = new_label;
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
		printf("Recommended threshold for holes.pgm = 70\n");
		printf("Recommended threshold for cells_small.pgm = 30\n");

		if (argc != 3)
		{
			printf("\nError! Input correct number of input parameters!\n");
			printf("\nCorrect usage: Homework <filename> <threshold>\n");
			exit(0);
		}

		//Read sigma value
		float sigma = 1;

		//Read image name
		string argument1 = argv[1];
		string path = "../../images/";
		string filename1 = path + argument1;
		const char *file1 = filename1.c_str();

		const char* argument2 = argv[2];
		int threshold = atoi(argument2);
		printf("Threshold entered = %d\n", threshold);

		ImgGray inputImage;
		Load(file1, &inputImage);

		ImgBgr finalImage;
		Load(file1, &finalImage);

		Figure fig1;
		fig1.SetTitle("Original Image");
		fig1.Draw(inputImage);

		//calculate width
		int mu = round(2.5*sigma - 0.5);
		int w = 2 * mu + 1;
		double *gauss;
		double *gaussDeriv;
		double sumGauss = 0;
		double sumGaussDeriv = 0;

		gauss = (double*)calloc(w, sizeof(double));
		gaussDeriv = (double*)calloc(w, sizeof(double));

		//	printf("width = %d\n", w);
		//	printf("mu = %d\n", mu);

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
		ImgFloat gradMagnitude;
		gradMagnitude.Reset(inputImage.Width(), inputImage.Height());
		Set(&gradMagnitude, 0);

		//Compute Gradient magnitude for input
		for (int y = mu; y < inputImage.Height() - mu; y++)
		{
			for (int x = mu; x < inputImage.Width() - mu; x++)
			{
				gradMagnitude(x, y) = max(abs(grad_x(x, y)), abs(grad_y(x, y)));
			}

		}

		//Normalize gradient image to 0-255 range
		ImgGray normalGrad;
		normalGrad.Reset(inputImage.Width(), inputImage.Height());
		normalGrad.Reset(inputImage.Width(), inputImage.Height());
		Set(&normalGrad, 0);

		float min = gradMagnitude(0, 0);
		float max = gradMagnitude(0, 0);

		for (int y = 0; y < inputImage.Height(); y++)
		{
			for (int x = 1; x < inputImage.Width(); x++)
			{
				if (gradMagnitude(x, y) > max)
				{
					max = gradMagnitude(x, y);
				}
				else if (gradMagnitude(x, y) < min)
				{
					min = gradMagnitude(x, y);
				}
			}
		}

		for (int y = 0; y < inputImage.Height(); y++)
		{
			for (int x = 1; x < inputImage.Width(); x++)
			{
				normalGrad(x, y) = (gradMagnitude(x, y) - min) * 255 / (max - min);
			}
		}


		//Find threshold for input image
		int highThreshold = threshold;

		ImgBinary highThresholdImg;
		highThresholdImg.Reset(inputImage.Width(), inputImage.Height());
		Set(&highThresholdImg, 0);


		for (int j = 0; j < inputImage.Height(); j++)
		{
			for (int i = 0; i < inputImage.Width(); i++)
			{
				if (inputImage(i, j) < highThreshold)
					highThresholdImg(i, j) = 1;
			}
		}

		//Chamfer image
		ImgInt chamferImage;
		chamferImage.Reset(inputImage.Width(), inputImage.Height());
		Set(&chamferImage, 0);
		GetChamferDistance(highThresholdImg, &chamferImage);

		//Find edges from watershed image
		ImgBinary edgeImage;
		edgeImage.Reset(inputImage.Width(), inputImage.Height());
		Set(&edgeImage, 0);

		ImgInt labelImage;
		labelImage.Reset(inputImage.Width(), inputImage.Height());
		Set(&labelImage, -1);

		//Watershed algorithm
		std::vector<Point> hist[256];
		for (int y = 0; y < inputImage.Height(); y++)
		{
			for (int x = 0; x < inputImage.Width(); x++)
			{
				//hist[inputImage(x, y)] += 1;
				hist[chamferImage(x, y)].push_back(Point(x, y));
			}
		}

		Point point;
		int globallabel = 0;
		std::queue<Point> frontier;
		for (int i = 0; i < 256; i++)
		{
			//GROW
			for (int j = 0; j < hist[i].size(); j++)
			{
				point = hist[i][j];
				{
					if (point.x < chamferImage.Width() - 1)
					{
						if (labelImage(point.x + 1, point.y) >= 0)
						{
							labelImage(point.x, point.y) = labelImage(point.x + 1, point.y);
							frontier.push(point);
						}
					}
					if (point.x >0)
					{
						if (labelImage(point.x - 1, point.y) >= 0)
						{
							labelImage(point.x, point.y) = labelImage(point.x - 1, point.y);
							frontier.push(point);
						}
					}
					if (point.y < chamferImage.Height() - 1)
					{
						if (labelImage(point.x, point.y + 1) >= 0)
						{
							labelImage(point.x, point.y) = labelImage(point.x, point.y + 1);
							frontier.push(point);
						}
					}
					if (point.y > 0)
					{
						if (labelImage(point.x, point.y - 1) >= 0)
						{
							labelImage(point.x, point.y) = labelImage(point.x, point.y - 1);
							frontier.push(point);
						}
					}
				}
			}

			//EXPAND
			while (!frontier.empty())
			{
				point = frontier.front();
				frontier.pop();
				{
					if (point.y > 0)
					{
						if (labelImage(point.x, point.y - 1) == -1 && chamferImage(point.x, point.y) == i)
						{
							labelImage(point.x, point.y - 1) = labelImage(point.x, point.y);
							frontier.push(point);
						}
					}
					if (point.y < chamferImage.Height() - 1)
					{
						if (labelImage(point.x, point.y + 1) == -1 && chamferImage(point.x, point.y) == i)
						{
							labelImage(point.x, point.y + 1) = labelImage(point.x, point.y);
							frontier.push(point);
						}
					}
					if (point.x > 0)
					{
						if (labelImage(point.x - 1, point.y) == -1 && chamferImage(point.x, point.y) == i)
						{
							labelImage(point.x - 1, point.y) = labelImage(point.x, point.y);
							frontier.push(point);
						}
					}
					if (point.x < chamferImage.Width() - 1)
					{
						if (labelImage(point.x + 1, point.y) == -1 && chamferImage(point.x, point.y) == i)
						{
							labelImage(point.x + 1, point.y) = labelImage(point.x, point.y);
							frontier.push(point);
						}
					}
				}
			}

			//CREATE	
			Point point1;
			for (int j = 0; j < hist[i].size(); j++)
			{
				point1 = hist[i][j];
				if (labelImage(point1.x, point1.y) < 0)
				{
					PerformFloodFill(chamferImage, point1.x, point1.y, &labelImage, globallabel);
					globallabel++;
				}
			}

		}

		for (int i = 0; i < globallabel; i++)
		{
			for (int y = 1; y < inputImage.Height() - 1; y++)
			{
				for (int x = 1; x < inputImage.Width() - 1; x++)
				{
					if ((labelImage(x - 1, y - 1) != labelImage(x, y)) ||
						(labelImage(x, y - 1) != labelImage(x, y)) ||
						(labelImage(x + 1, y - 1) != labelImage(x, y)) ||
						(labelImage(x - 1, y) != labelImage(x, y)) ||
						(labelImage(x + 1, y) != labelImage(x, y)) ||
						(labelImage(x - 1, y + 1) != labelImage(x, y)) ||
						(labelImage(x, y + 1) != labelImage(x, y)) ||
						(labelImage(x + 1, y + 1) != labelImage(x, y)))
					{
						edgeImage(x, y) = 1;
					}
				}
			}

		}


		ImgBinary markerImage;
		markerImage.Reset(inputImage.Width(), inputImage.Height());
		Set(&markerImage, 0);

		for (int y = 0; y < inputImage.Height(); y++)
		{
			for (int x = 0; x < inputImage.Width(); x++)
			{
				if ((edgeImage(x, y) == 1) || (highThresholdImg(x, y) == 1))
				{
					markerImage(x, y) = 1;
				}
			}
		}


		//Marker based watersheding
		ImgInt labels;
		labels.Reset(inputImage.Width(), inputImage.Height());
		Set(&labels, -1);

		ImgInt gradMagInt;
		gradMagInt.Reset(inputImage.Width(), inputImage.Height());
		Convert(normalGrad, &gradMagInt);

		std::vector<Point> list[256];
		for (int y = 0; y < inputImage.Height(); y++)
		{
			for (int x = 0; x < inputImage.Width(); x++)
			{
				//list[inputImage(x, y)] += 1;
				list[normalGrad(x, y)].push_back(Point(x, y));
			}
		}



		Point pixel;
		int global_label = 0;
		std::queue<Point> frontier2;
		for (int i = 0; i < 256; i++)
		{
			//GROW
			for (int j = 0; j < list[i].size(); j++)
			{
				pixel = list[i][j];
				{
					if (pixel.x < labels.Width() - 1)
					{
						if (labels(pixel.x + 1, pixel.y) >= 0)
						{
							labels(pixel.x, pixel.y) = labels(pixel.x + 1, pixel.y);
							frontier2.push(pixel);
						}
					}
					if (pixel.x > 0)
					{
						if (labels(pixel.x - 1, pixel.y) >= 0)
						{
							labels(pixel.x, pixel.y) = labels(pixel.x - 1, pixel.y);
							frontier2.push(pixel);
						}
					}
					if (pixel.y < labels.Height() - 1)
					{
						if (labels(pixel.x, pixel.y + 1) >= 0)
						{
							labels(pixel.x, pixel.y) = labels(pixel.x, pixel.y + 1);
							frontier2.push(pixel);
						}
					}
					if (pixel.y > 0)
					{
						if (labels(pixel.x, pixel.y - 1) >= 0)
						{
							labels(pixel.x, pixel.y) = labels(pixel.x, pixel.y - 1);
							frontier2.push(pixel);
						}
					}
				}
			}

			//EXPAND
			while (!frontier2.empty())
			{
				pixel = frontier2.front();
				frontier2.pop();
				{
					if (pixel.y > 0)
					{
						if (labels(pixel.x, pixel.y - 1) == -1 && gradMagInt(pixel.x, pixel.y) <= i)
						{
							labels(pixel.x, pixel.y - 1) = labels(pixel.x, pixel.y);
							frontier2.push(pixel);
						}
					}
					if (pixel.y < labels.Height() - 1)
					{
						if (labels(pixel.x, pixel.y + 1) == -1 && gradMagInt(pixel.x, pixel.y) <= i)
						{
							labels(pixel.x, pixel.y + 1) = labels(pixel.x, pixel.y);
							frontier2.push(pixel);
						}
					}
					if (pixel.x > 0)
					{
						if (labels(pixel.x - 1, pixel.y) == -1 && gradMagInt(pixel.x, pixel.y) <= i)
						{
							labels(pixel.x - 1, pixel.y) = labels(pixel.x, pixel.y);
							frontier2.push(pixel);
						}
					}
					if (pixel.x < labels.Width() - 1)
					{
						if (labels(pixel.x + 1, pixel.y) == -1 && gradMagInt(pixel.x, pixel.y) <= i)
						{
							labels(pixel.x + 1, pixel.y) = labels(pixel.x, pixel.y);
							frontier2.push(pixel);
						}
					}
				}
			}

			ImgInt marker;
			marker.Reset(inputImage.Width(), inputImage.Height());
			Convert(markerImage, &marker);

			//CREATE
			Point pix;
			for (int j = 0; j < list[i].size(); j++)
			{
				pix = list[i][j];
				if (markerImage(pix.x, pix.y) == 1)
				{
					PerformFloodFill(marker, pix.x, pix.y, &labels, global_label);
					global_label++;
				}
			}

		}



		for (int g = 0; g < global_label; g++)
		{
			for (int y = 1; y < inputImage.Height() - 1; y++)
			{
				for (int x = 1; x < inputImage.Width() - 1; x++)
				{
					if (labels(x, y) == g)
					{
						if ((labels(x, y - 1) != labels(x, y)) ||
							(labels(x - 1, y) != labels(x, y)) ||
							(labels(x + 1, y) != labels(x, y)) ||
							(labels(x, y + 1) != labels(x, y)))
						{
							if (highThresholdImg(x, y) != 0)
							{
								finalImage(x, y) = Bgr(0, 255, 0);
							}
						}
					}

				}
			}
		}

		//Output Images
		Figure fig4, fig7, fig10, fig12, fig13, fig15, fig16, fig18;

		fig10.SetTitle("high threshold image");
		fig10.Draw(highThresholdImg);

		fig4.SetTitle("Magnitude of Gradient");
		fig4.Draw(gradMagnitude);

		fig13.SetTitle("normalized gradient image");
		fig13.Draw(normalGrad);

		fig7.SetTitle("Chamfer distance image");
		fig7.Draw(chamferImage);

		fig16.SetTitle("Non-marker watershed image");
		fig16.Draw(labelImage);

		fig15.SetTitle("edge image");
		fig15.Draw(edgeImage);

		fig12.SetTitle("marker image");
		fig12.Draw(markerImage);

		fig18.SetTitle("final image");
		fig18.Draw(finalImage);

		EventLoop();
	}
	catch (const Exception& e)
	{
		e.Display();    // display exception to user in a popup window 
	}
	return 0;
}
