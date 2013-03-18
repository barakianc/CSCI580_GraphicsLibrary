/*   CS580 HW   */
#include    "stdafx.h"  
#include	"Gz.h"
#include	"disp.h"


int GzNewFrameBuffer(char** framebuffer, int width, int height)
{
/* create a framebuffer:
 -- allocate memory for framebuffer : (sizeof)GzPixel x width x height
 -- pass back pointer 
*/
	if(width <= MAXXRES && height <= MAXYRES){
		*framebuffer = new char[width*height*3];
	}
	else{
		return GZ_FAILURE;
	}

	return GZ_SUCCESS;
}

int GzNewDisplay(GzDisplay	**display, GzDisplayClass dispClass, int xRes, int yRes)
{

/* create a display:
  -- allocate memory for indicated class and resolution
  -- pass back pointer to GzDisplay object in display
*/
	if(xRes <= MAXXRES && yRes <= MAXYRES){
		*display = (GzDisplay*)malloc(sizeof(GzDisplay));
		(*display)->xres = xRes;
		(*display)->yres = yRes;
		(*display)->dispClass = dispClass;
		(*display)->open = 0;	//open set to 1 when properly intialize
		(*display)->fbuf = new GzPixel[xRes*yRes];	//needs to be initialized
		
	}
	else{
		return GZ_FAILURE;
	}

	return GZ_SUCCESS;
}


int GzFreeDisplay(GzDisplay	*display)
{
/* clean up, free memory */
	if( display != NULL){
		delete [] display->fbuf;	//free the array
		free (display);				//now free the struct
	}
	else{
		return GZ_FAILURE;
	}
	return GZ_SUCCESS;
}


int GzGetDisplayParams(GzDisplay *display, int *xRes, int *yRes, GzDisplayClass	*dispClass)
{
/* pass back values for an open display */
	if(display != NULL){
		*xRes = display->xres;
		*yRes = display->yres;
		*dispClass = display->dispClass;
	}
	else{
		return GZ_FAILURE;
	}
	return GZ_SUCCESS;
}


int GzInitDisplay(GzDisplay	*display)
{
/* set everything to some default values - start a new frame */
	if(display != NULL){
		display->open = 1;
		int xres = display->xres;
		int yres = display->yres;
		for(int i = 0; i<xres*yres;i++){
			display->fbuf[i].red = 4095;
			display->fbuf[i].green = 0;
			display->fbuf[i].blue = 0;
			display->fbuf[i].alpha = 0;
			display->fbuf[i].z = MAXINT;
		}
	}
	else{
		return GZ_FAILURE;
	}
	return GZ_SUCCESS;
}


int GzPutDisplay(GzDisplay *display, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
/* write pixel values into the display */
	if(display != NULL && display->open == 1){
		if(0<= i && i< display->xres && 0<= j && j< display->yres){
			int pixel = i+(j*display->xres);
			//Set the red value
			if(r<0){
				display->fbuf[pixel].red = 0;
			}
			else if(r>4095){
				display->fbuf[pixel].red = 4095;
			}
			else{
				display->fbuf[pixel].red = r;
			}
			//set the green value
			if(g<0){
				display->fbuf[pixel].green = 0;
			}
			else if(g>4095){
				display->fbuf[pixel].green = 4095;
			}
			else{
				display->fbuf[pixel].green = g;
			}
			//set the blue value
			if(b<0){
				display->fbuf[pixel].blue = 0;
			}
			else if(b > 4095){
				display->fbuf[pixel].blue = 4095;
			}
			else{
				display->fbuf[pixel].blue = b;
			}
			//set alpha
			display->fbuf[pixel].alpha = a;
			//set z depth
			display->fbuf[pixel].z = z;
		}
	}
	else{
		return GZ_FAILURE;
	}
	return GZ_SUCCESS;
}


int GzGetDisplay(GzDisplay *display, int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
	/* pass back pixel value in the display */
	/* check display class to see what vars are valid */
	if(display != NULL && display->open == 1){
		if(0<= i && i < display->xres && 0<= j && j< display->yres){
			int pixel = i+(j*display->xres);
			*r = display->fbuf[pixel].red;
			*g = display->fbuf[pixel].green;
			*b = display->fbuf[pixel].blue;
			*a = display->fbuf[pixel].alpha;
			*z = display->fbuf[pixel].z;
		}
		else{
			*r = 0;
			*g = 0;
			*b = 0;
			*a = 0;
			*z = 0;
		}
	}
	else{
		return GZ_FAILURE;
	}
	return GZ_SUCCESS;
}


int GzFlushDisplay2File(FILE* outfile, GzDisplay *display)
{

	/* write pixels to ppm file based on display class -- "P6 %d %d 255\r" */
	if(outfile != NULL && display != NULL && display->open == 1){
		char head[20];
		int n;
		//write the header for the PPM file
		n = sprintf_s(head,"P6 %d %d 255\n",display->xres,display->yres);
		fputs(head,outfile);
		GzIntensity tempR,tempG,tempB;
		//Now go through the displays frame buffer and write the pixels to the file
		for(int i = 0; i < display->xres*display->yres;i++){
			tempR = display->fbuf[i].red;
			tempG = display->fbuf[i].green;
			tempB = display->fbuf[i].blue;
			tempR = tempR>>4;
			tempG = tempG>>4;
			tempB = tempB>>4;
			char r = (char)tempR;
			fputc((char)tempR,outfile);
			fputc((char)tempG,outfile);
			fputc((char)tempB,outfile);
		}
	}
	else{
		return GZ_FAILURE;
	}

	return GZ_SUCCESS;
}

int GzFlushDisplay2FrameBuffer(char* framebuffer, GzDisplay *display)
{

	/* write pixels to framebuffer: 
		- Put the pixels into the frame buffer
		- Caution: store the pixel to the frame buffer as the order of blue, green, and red 
		- Not red, green, and blue !!!
	*/
	if(framebuffer != NULL && display != NULL && display->open == 1){
		GzIntensity tempR,tempG,tempB;
		int j =0;
		for(int i = 0; i< display->xres*display->yres;i++){
			tempR = display->fbuf[i].red;
			tempG = display->fbuf[i].green;
			tempB = display->fbuf[i].blue;
			tempR = tempR>>4;
			tempG = tempG>>4;
			tempB = tempB>>4;
			framebuffer[j] = (char)tempB;
			framebuffer[j+1] = (char)tempG;
			framebuffer[j+2] = (char)tempR;
			j+=3;
		}
	}
	else{
		return GZ_FAILURE;
	}

	return GZ_SUCCESS;
}