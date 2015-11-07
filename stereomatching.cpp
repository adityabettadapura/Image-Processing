// Homework.cpp : Defines the entry point for the console application.
// Stereo matching algorithm
// Author: Aditya Bettadapura

#include <afxwin.h>  // necessary for MFC to work properly
#include "Homework.h"
#include "../../src/blepo.h"
#include "math.h"
#include <stack>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

using namespace blepo;

void Convolution(ImgInt* input, ImgInt* output, int w, int h, int size)
{
	//Convolution with separable kernels
	for (int y = 0; y < h; y++)
	{
		for (int x = size; x < w - size; x++)
		{
			int sum = 0;
			for (int i = -size; i <= size; i++)
			{
				sum += (*input)(x + i, y);
			}
			(*output)(x, y) = sum;
		}
	}
	for (int y = size; y < h - size; y++)
	{
		for (int x = 0; x < w; x++)
		{
			int sum = 0;
			for (int j = -size; j <= size; j++)
			{
				sum += (*output)(x, y + j);
			}
			(*input)(x, y) = sum / ((2 * size + 1)*(2 * size + 1));
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
		CString filename1, filename2, path;

		path = "../../images/";

		filename1 = argv[1];
		filename1 = path + filename1;

		filename2 = argv[2];
		filename2 = path + filename2;

		int max_disp = atoi(argv[3]);
		int min_disp = 0;

		ImgBgr leftinput, rightinput;
		Load(filename1, &leftinput);
		Load(filename2, &rightinput);

		if ((leftinput.Height() != rightinput.Height()) || (leftinput.Width() != rightinput.Width()))
		{
			printf("ERROR: Images are of different sizes!\n");
			exit(0);
		}

		Figure fig1;
		fig1.SetTitle("Original Left Image");
		fig1.Draw(leftinput);

		Figure fig2;
		fig2.SetTitle("Original Right Image");
		fig2.Draw(rightinput);

		ImgGray leftimg, rightimg;

		Convert(leftinput, &leftimg);
		Convert(rightinput, &rightimg);

		Figure fig3;
		fig3.SetTitle("Left Image");
		fig3.Draw(leftimg);

		Figure fig4;
		fig4.SetTitle("Right Image");
		fig4.Draw(rightimg);

		printf("Max disparity entered = %d\n", max_disp);

		//Take vector of Integer images
		std::vector<ImgInt> dbar(max_disp + 1);
		for (int i = 0; i <= max_disp; i++)
		{
			dbar[i].Reset(leftinput.Width(), leftinput.Height());
			Set(&dbar[i], 0);
		}

		ImgInt tempimg;
		tempimg.Reset(leftinput.Width(), leftinput.Height());

		ImgInt dleft, dright, dmap;
		dleft.Reset(leftinput.Width(), leftinput.Height());
		Set(&dleft, 0);
		dright.Reset(leftinput.Width(), leftinput.Height());
		Set(&dright, 0);
		dmap.Reset(leftinput.Width(), leftinput.Height());
		Set(&dmap, 0);

		//For win_size=4, Convolution window size = 9x9
		int win_size = 4;
		
		//BlockMatch
		{
			//PreCompute dissimilarities
			for (int d = min_disp; d <= max_disp; d++)
			{
				for (int y = 0; y < leftimg.Height(); y++)
				{
					for (int x = d; x < leftimg.Width(); x++)
					{
						dbar[d](x, y) = abs(leftimg(x, y) - rightimg(x - d, y));
					}
				}

				Convolution(&dbar[d], &tempimg, leftinput.Width(), leftinput.Height(), win_size);

			}

			int ghat, dhat;

			//dleft - disparity map for left image
			for (int y = 0; y < leftimg.Height(); y++)
			{
				for (int x = max_disp; x < leftimg.Width() - max_disp; x++)
				{
					ghat = 1000;
					for (int d = min_disp; d < max_disp; d++)
					{
						if (dbar[d](x, y) < ghat)
						{
							ghat = dbar[d](x, y);
							dhat = d;
						}
					}
					dleft(x, y) = dhat;
				}
			}

			//dright - disparity map for right image
			for (int y = 0; y < leftimg.Height(); y++)
			{
				for (int x = max_disp; x < leftimg.Width() - max_disp; x++)
				{
					ghat = 1000;
					for (int d = min_disp; d < max_disp; d++)
					{
						if (dbar[d](x - dleft(x, y) + d, y) < ghat)
						{
							ghat = dbar[d](x - dleft(x, y) + d, y);
							dhat = d;
						}
					}
					dright(x, y) = dhat;
					dmap(x, y) = (dleft(x, y) == dright(x, y)) ? dleft(x, y):0;

				}
			}
		}

		ImgFloat depth;
		depth.Reset(leftinput.Width(), leftinput.Height());
		Set(&depth, 0);

		//Calculate depth and output ply file
		FILE *fptr;
		fptr = fopen("output.ply", "wb");
		fprintf(fptr, "ply\nformat ascii 1.0\nelement vertex %d\nproperty float x\nproperty float y\nproperty float z\nproperty uchar diffuse_red\nproperty uchar diffuse_green\nproperty uchar diffuse_blue\nend_header\n", leftinput.Width()*leftinput.Height());

		for (int y = 0; y < leftinput.Height(); y++)
		{
			for (int x = 0; x < leftinput.Width(); x++)
			{
				if (dmap(x,y) != 0)
				{
					depth(x, y) = (float)(1000 /dmap(x, y));
				}
				else
				{
					depth(x, y) = 1000;
				}
				fprintf(fptr, "%d %d %f %d %d %d\n", x, y, depth(x, y), leftinput(x, y).r, leftinput(x, y).g, leftinput(x, y).b);
			}
		}
		fclose(fptr);

		Figure fig5;
		fig5.SetTitle("Disparity Map");
		fig5.Draw(dleft);

		Figure fig6;
		fig6.SetTitle("left_right consistency check");
		fig6.Draw(dmap);

		EventLoop();
	}
	catch (const Exception& e)
	{
		e.Display();    // display exception to user in a popup window 
	}
	return 0;
}
