#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"


//DDA struct for edges
struct DDA{
	GzCoord start;
	GzCoord end;
	GzCoord current;
	float Xslope;
	float Zslope;
};

//span DDA struct
struct SPANDDA{
	GzCoord start;
	GzCoord end;
	GzCoord current;
	float Zslope;
};


/* NOT part of API - just for general assistance */

short	ctoi(float color)		/* convert float color to GzIntensity short */
{
  return(short)((int)(color * ((1 << 12) - 1)));
}

//Z interpolation helper function takes in (x,y) A,B,C,D and spits out a z
float	interpZ(float x,float y,float A,float B,float C,float D){
	float intZ = 0;
	intZ = ((-D)+(-A*x)+(-B*y))/C;
	return intZ;
}

int GzNewRender(GzRender **render, GzRenderClass renderClass, GzDisplay *display)
{
/* 
- malloc a renderer struct
- keep closed until BeginRender inits are done
- span interpolator needs pointer to display for pixel writes
- check for legal class GZ_Z_BUFFER_RENDER
*/

	if(display != NULL && renderClass == GZ_Z_BUFFER_RENDER ){
		*render = (GzRender*)malloc(sizeof(GzRender));
		(*render)->renderClass = renderClass;
		(*render)->display = display;
		(*render)->open = 0;
	}
	else{
		return GZ_FAILURE;
	}
	return GZ_SUCCESS;
}


int GzFreeRender(GzRender *render)
{
/* 
-free all renderer resources
*/
	if(render != NULL){
		render->display = NULL;
		free(render);
	}
	else{
		return GZ_FAILURE;
	}

	return GZ_SUCCESS;
}


int GzBeginRender(GzRender	*render)
{
/* 
- set up for start of each frame - init frame buffer
*/

	
	if(render != NULL){
		render->open = 1;
		//init Ximage array and Xnorm array
		for(int i = 0; i < MATLEVELS; i++){
			for(int j = 0; j < 4; j++){
				for(int k = 0; k < 4; k++){
					(render->Ximage[i])[j][k] = 0;
					(render->Xnorm[i])[j][k] = 0;
				}
			}
		}

		//init the Xsp matrix
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				(render->Xsp)[i][j] = 0;
			}
		}

		//init the flat render color to blue
		render->flatcolor[0] = 0;
		render->flatcolor[1] = 0;
		render->flatcolor[2] = 1.0;
		
		render->interp_mode = 0;
		render->numlights = 0;

		//init Ka to 0
		render->Ka[0] = 0;
		render->Ka[1] = 0;
		render->Ka[2] = 0;

		//init Kd to 0
		render->Kd[0] = 0;
		render->Kd[1] = 0;
		render->Kd[2] = 0;

		//init Ks to 0
		render->Ks[0] = 0;
		render->Ks[1] = 0;
		render->Ks[2] = 0;

		//init spec to 0
		render->spec = 0;
	}
	else{
		return GZ_FAILURE;
	}
	return GZ_SUCCESS;
}


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer *valueList) /* void** valuelist */
{
/*
- set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
- later set shaders, interpolaters, texture maps, and lights
*/
	if(render != NULL && render->open == 1 && numAttributes > 0 && nameList != NULL && valueList != NULL){
		//loop through each of the attributes and set the corresponding attributes
		for(int i = 0; i< numAttributes; i++){
			switch(nameList[i]){
			case(GZ_RGB_COLOR):
				{
					GzPointer tempPointer = valueList[i];
					GzColor tempColor = {((float*)tempPointer)[0],((float*)tempPointer)[1],((float*)tempPointer)[2]};

					//set R color
					if( tempColor[0] < 0.0){
						render->flatcolor[0] = 0.0;
					}
					else if(tempColor[0] > 1.0){
						render->flatcolor[0] = 1.0;
					}
					else{
						render->flatcolor[0] = tempColor[0];
					}

					//set G color
					if( tempColor[1] < 0.0){
						render->flatcolor[1] = 0.0;
					}
					else if(tempColor[1] > 1.0){
						render->flatcolor[1] = 1.0;
					}
					else{
						render->flatcolor[1] = tempColor[1];
					}

					//set B color
					if( tempColor[2] < 0.0){
						render->flatcolor[2] = 0.0;
					}
					else if(tempColor[2] > 1.0){
						render->flatcolor[2] = 1.0;
					}
					else{
						render->flatcolor[2] = tempColor[2];
					}

					break;
				}
			}
		}
	}
	else{
		return GZ_FAILURE;
	}

	return GZ_SUCCESS;
}


int GzPutTriangle(GzRender *render, int	numParts, GzToken *nameList,
	GzPointer *valueList) 
/* numParts - how many names and values */
{
/* 
- pass in a triangle description with tokens and values corresponding to
      GZ_NULL_TOKEN:		do nothing - no values
      GZ_POSITION:		3 vert positions in model space
- Invoke the scan converter and return an error code
*/
	if(render != NULL && numParts >= 0 && nameList != NULL && valueList != NULL && render->open == 1){
		for(int i = 0; i< numParts; i++){
			switch(nameList[i]){
			case(GZ_POSITION):
				{
					//for each triangle we to determine the pixel that make up that triangle and set the color accordingly
					//We also need to interpolate z for each pixel to do z-buffer checking

					GzPointer tempPointer = valueList[i];
					//get the coords for the 3 verts
					GzCoord V1 = {((GzCoord*)tempPointer)[0][0],((GzCoord*)tempPointer)[0][1],((GzCoord*)tempPointer)[0][2]};
					GzCoord V2 = {((GzCoord*)tempPointer)[1][0],((GzCoord*)tempPointer)[1][1],((GzCoord*)tempPointer)[1][2]};
					GzCoord V3 = {((GzCoord*)tempPointer)[2][0],((GzCoord*)tempPointer)[2][1],((GzCoord*)tempPointer)[2][2]};

					//sort verts by Y coord, least to greatest
					GzCoord vertsList[3];
					memcpy(vertsList[0],V1,3*sizeof(float));
					memcpy(vertsList[1],V2,3*sizeof(float));
					memcpy(vertsList[2],V3,3*sizeof(float));
					bool changes = true;
					GzCoord tempCoord;
					while(changes){
						changes = false;
						for(int j = 0; j < 2 ; j++){
							if(vertsList[j][1] > vertsList[j+1][1]){
								memcpy(tempCoord,vertsList[j+1],3*sizeof(float));
								memcpy(vertsList[j+1],vertsList[j],3*sizeof(float));
								memcpy(vertsList[j],tempCoord,3*sizeof(float));
								changes = true;
							}
						}
					}

					if(abs(vertsList[2][1] - vertsList[1][1]) < .001){
						int test = 0;
					}
					DDA* Lverts[2];
					DDA* Rverts[2];
					//if the triangle is L(1-3) R(1-2 2-3)
					float s1 = (vertsList[1][0]-vertsList[0][0])/(vertsList[1][1]-vertsList[0][1]);
					float s2 = (vertsList[2][0]-vertsList[0][0])/(vertsList[2][1]-vertsList[0][1]);
					if((((vertsList[1][0]-vertsList[0][0])/(vertsList[1][1]-vertsList[0][1]) < (vertsList[2][0]-vertsList[0][0])/(vertsList[2][1]-vertsList[0][1])) && vertsList[1][0] < vertsList[2][0]) && (vertsList[1][0]-vertsList[0][0])/(vertsList[1][1]-vertsList[0][1])> 0 && (vertsList[2][0]-vertsList[0][0])/(vertsList[2][1]-vertsList[0][1]) > 0){
						int test = 0;
					}
					if(vertsList[1][1] == vertsList[2][1]){
						int test = 0;
					}
					else if (((vertsList[1][0]-vertsList[0][0])/(vertsList[1][1]-vertsList[0][1]) > (vertsList[2][0]-vertsList[0][0])/(vertsList[2][1]-vertsList[0][1]))){
						//L edge DDA
						Lverts[0] = new DDA();
						memcpy(Lverts[0]->start,vertsList[0],3*sizeof(float));
						memcpy(Lverts[0]->end,vertsList[2],3*sizeof(float));
						memcpy(Lverts[0]->current,vertsList[0],3*sizeof(float));
						Lverts[0]->Xslope = (Lverts[0]->end[0]-Lverts[0]->start[0])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[0]->Zslope = (Lverts[0]->end[2]-Lverts[0]->start[2])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[1] = NULL;

						//R edges DDAs
						Rverts[0] = new DDA();
						memcpy(Rverts[0]->start,vertsList[0],3*sizeof(float));
						memcpy(Rverts[0]->end,vertsList[1],3*sizeof(float));
						memcpy(Rverts[0]->current, vertsList[0],3*sizeof(float));
						Rverts[0]->Xslope = (Rverts[0]->end[0]-Rverts[0]->start[0])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->Zslope = (Rverts[0]->end[2]-Rverts[0]->start[2])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[1] = new DDA();
						memcpy(Rverts[1]->start,vertsList[1],3*sizeof(float));
						memcpy(Rverts[1]->end,vertsList[2],3*sizeof(float));
						memcpy(Rverts[1]->current,vertsList[1],3*sizeof(float));
						Rverts[1]->Xslope = (Rverts[1]->end[0]-Rverts[1]->start[0])/(Rverts[1]->end[1]-Rverts[1]->start[1]);
						Rverts[1]->Zslope = (Rverts[1]->end[2]-Rverts[1]->start[2])/(Rverts[1]->end[1]-Rverts[1]->start[1]);

					}
					//if  the triangle is L(1-2 2-3) R(1-3)
					else if(((vertsList[1][0]-vertsList[0][0])/(vertsList[1][1]-vertsList[0][1]) < (vertsList[2][0]-vertsList[0][0])/(vertsList[2][1]-vertsList[0][1]) )){
						//L edges DDAS
						Lverts[0] = new DDA();
						memcpy(Lverts[0]->start,vertsList[0],3*sizeof(float));
						memcpy(Lverts[0]->end,vertsList[1],3*sizeof(float));
						memcpy(Lverts[0]->current,vertsList[0],3*sizeof(float));
						Lverts[0]->Xslope = (Lverts[0]->end[0]-Lverts[0]->start[0])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[0]->Zslope = (Lverts[0]->end[2]-Lverts[0]->start[2])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[1] = new DDA();
						memcpy(Lverts[1]->start,vertsList[1],3*sizeof(float));
						memcpy(Lverts[1]->end,vertsList[2],3*sizeof(float));
						memcpy(Lverts[1]->current,vertsList[1],3*sizeof(float));
						//float tempfl = Lverts[1]->end;//-*Lverts[1]->start[1]);///(*Lverts[1]->end[1]-*Lverts[1]->start[1]);
						Lverts[1]->Xslope = (Lverts[1]->end[0]-Lverts[1]->start[0])/(Lverts[1]->end[1]-Lverts[1]->start[1]);
						Lverts[1]->Zslope = (Lverts[1]->end[2]-Lverts[1]->start[2])/(Lverts[1]->end[1]-Lverts[1]->start[1]);

						//R edge DDA
						Rverts[0] = new DDA();
						memcpy(Rverts[0]->start,vertsList[0],3*sizeof(float));
						memcpy(Rverts[0]->end,vertsList[2],3*sizeof(float));
						memcpy(Rverts[0]->current,vertsList[0],3*sizeof(float));
						Rverts[0]->Xslope = (Rverts[0]->end[0]-Rverts[0]->start[0])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->Zslope = (Rverts[0]->end[2]-Rverts[0]->start[2])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[1] = NULL;

					}

					//interp z need to find A B C D
					float A =0 ,B = 0,C = 0,D = 0;
					GzCoord ed1,ed2;
					ed1[0] = vertsList[1][0]-vertsList[0][0];
					ed1[1] = vertsList[1][1]-vertsList[0][1];
					ed1[2] = vertsList[1][2]-vertsList[0][2];
					ed2[0] = vertsList[2][0]-vertsList[0][0];
					ed2[1] = vertsList[2][1]-vertsList[0][1];
					ed2[2] = vertsList[2][2]-vertsList[0][2];

					A = (ed1[1]*ed2[2])-(ed1[2]*ed2[1]);
					B = (ed1[2]*ed2[0])-(ed1[0]*ed2[2]);
					C = (ed1[0]*ed2[1])-(ed1[1]*ed2[0]);
					D = (-A*vertsList[0][0])+(-B*vertsList[0][1])+(-C*vertsList[0][2]);

					//Now perform the Scanline algorithm
					int minY = 0 , MaxY = 0,midY,minX = 0,MaxX = 0;
					if(Lverts[1] == NULL){
						minY = Lverts[0]->start[1];
						MaxY = Lverts[0]->end[1];
						midY = Rverts[1]->start[1];
						//first move down the scanline
						float dY;
						for(int yLine = minY+1; yLine < MaxY+1; yLine++){
							//for first line calc delta Y after that Y changes by one each time
							if(yLine == minY+1 && minY != midY){
								
								dY = minY+1 - Lverts[0]->start[1];
								Lverts[0]->current[0] = Lverts[0]->current[0]+(Lverts[0]->Xslope*dY);
								Lverts[0]->current[1] = Lverts[0]->current[1]+dY;
								Lverts[0]->current[2] = Lverts[0]->current[2]+(Lverts[0]->Zslope*dY);
								Rverts[0]->current[0] = Rverts[0]->current[0]+(Rverts[0]->Xslope*dY);
								Rverts[0]->current[1] = Rverts[0]->current[1]+dY;
								Rverts[0]->current[2] = Rverts[0]->current[2]+(Rverts[0]->Zslope*dY);

								//set up span DDA
								minX = Lverts[0]->current[0];
								MaxX = Rverts[0]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[0]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[0]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[0]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);
								float dX;
								if(minX < MaxX){
									for(int xLine = minX; xLine < MaxX + 1;xLine++){
										//first time 
										//if start falls on line
										if(xLine == minX && spanX.current[0] == xLine){
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0]+1 != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								/*else{
											float  intZ = interpZ(MaxX,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,MaxX,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,MaxX,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
								}*/

							}
							else if(yLine <= midY){
								dY = 1;
								Lverts[0]->current[0] = Lverts[0]->current[0]+(Lverts[0]->Xslope*dY);
								Lverts[0]->current[1] = Lverts[0]->current[1]+dY;
								Lverts[0]->current[2] = Lverts[0]->current[2]+(Lverts[0]->Zslope*dY);
								Rverts[0]->current[0] = Rverts[0]->current[0]+(Rverts[0]->Xslope*dY);
								Rverts[0]->current[1] = Rverts[0]->current[1]+dY;
								Rverts[0]->current[2] = Rverts[0]->current[2]+(Rverts[0]->Zslope*dY);

								//set up span DDA
								minX = Lverts[0]->current[0];
								MaxX = Rverts[0]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[0]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[0]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[0]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);
								float dX;
								if(minX < MaxX){
									for(int xLine = minX; xLine < MaxX + 1;xLine++){
										//first time 
										if(xLine == minX && spanX.current[0] == xLine){
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0]+1 != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								/*else{
											float  intZ = interpZ(MaxX,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,MaxX,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,MaxX,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
								}*/
							}
							else if(minY == midY && yLine == minY+1){
								float dYL = midY - Rverts[1]->start[1]+1;
								float dYR = midY - Lverts[0]->start[1]+1;
								Lverts[0]->current[0] = Lverts[0]->current[0]+(Lverts[0]->Xslope*dYL);
								Lverts[0]->current[1] = Lverts[0]->current[1]+dYL;
								Lverts[0]->current[2] = Lverts[0]->current[2]+(Lverts[0]->Zslope*dYL);
								Rverts[1]->current[0] = Rverts[1]->current[0]+(Rverts[1]->Xslope*dYR);
								Rverts[1]->current[1] = Rverts[1]->current[1]+dYR;
								Rverts[1]->current[2] = Rverts[1]->current[2]+(Rverts[1]->Zslope*dYR);

								//set up span DDA
								minX = Lverts[0]->current[0];
								MaxX = Rverts[1]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[0]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[1]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[0]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);
								float dX;
								if(minX < MaxX){
									for(int xLine = minX; xLine < MaxX + 1;xLine++){
										//first time
										if(xLine == minX && spanX.current[0] == xLine){
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0] != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								/*else{
											float  intZ = interpZ(MaxX,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,MaxX,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,MaxX,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
								}*/
							}
							else if((yLine == midY+1 && midY != MaxY)){
								dY = midY - Rverts[1]->start[1]+1;
								Lverts[0]->current[0] = Lverts[0]->current[0]+(Lverts[0]->Xslope*1);
								Lverts[0]->current[1] = Lverts[0]->current[1]+1;
								Lverts[0]->current[2] = Lverts[0]->current[2]+(Lverts[0]->Zslope*1);
								Rverts[1]->current[0] = Rverts[1]->current[0]+(Rverts[1]->Xslope*dY);
								Rverts[1]->current[1] = Rverts[1]->current[1]+dY;
								Rverts[1]->current[2] = Rverts[1]->current[2]+(Rverts[1]->Zslope*dY);

								//set up span DDA
								minX = Lverts[0]->current[0];
								MaxX = Rverts[1]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[0]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[1]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[0]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);
								float dX;
								if(minX < MaxX){
									for(int xLine = minX; xLine < MaxX + 1;xLine++){
										//first time
										if(xLine == minX && spanX.current[0] == xLine){
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0]+1 != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								/*else{
											float  intZ = interpZ(MaxX,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,MaxX,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,MaxX,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
								}*/
							}
							else if(midY != MaxY && yLine > midY){
								dY = 1;
								Lverts[0]->current[0] = Lverts[0]->current[0]+(Lverts[0]->Xslope*dY);
								Lverts[0]->current[1] = Lverts[0]->current[1]+dY;
								Lverts[0]->current[2] = Lverts[0]->current[2]+(Lverts[0]->Zslope*dY);
								Rverts[1]->current[0] = Rverts[1]->current[0]+(Rverts[1]->Xslope*dY);
								Rverts[1]->current[1] = Rverts[1]->current[1]+dY;
								Rverts[1]->current[2] = Rverts[1]->current[2]+(Rverts[1]->Zslope*dY);

								//set up span DDA
								minX = Lverts[0]->current[0];
								MaxX = Rverts[1]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[0]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[1]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[0]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);
								float dX;
								if(minX < MaxX){
									for(int xLine = minX; xLine < MaxX +1;xLine++){
										//first time
										if(xLine == minX && spanX.current[0] == xLine){
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0]+1 != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								/*else{
											float  intZ = interpZ(MaxX,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,MaxX,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,MaxX,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
								}*/
							}

						}
					}
					else{
						minY = Lverts[0]->start[1];
						MaxY = Rverts[0]->end[1];
						midY = Lverts[1]->start[1];
						//first move down the scanline
						float dY;
						for(int yLine = minY+1; yLine < MaxY+1; yLine++){
							//for first line calc delta Y after that Y changes by one each time
							if(yLine == minY+1 && minY != midY){
								dY = minY+1 - Lverts[0]->start[1];
								Lverts[0]->current[0] = Lverts[0]->current[0]+(Lverts[0]->Xslope*dY);
								Lverts[0]->current[1] = Lverts[0]->current[1]+dY;
								Lverts[0]->current[2] = Lverts[0]->current[2]+(Lverts[0]->Zslope*dY);
								Rverts[0]->current[0] = Rverts[0]->current[0]+(Rverts[0]->Xslope*dY);
								Rverts[0]->current[1] = Rverts[0]->current[1]+dY;
								Rverts[0]->current[2] = Rverts[0]->current[2]+(Rverts[0]->Zslope*dY);

								//set up span DDA
								minX = Lverts[0]->current[0];
								MaxX = Rverts[0]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[0]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[0]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[0]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);
								float dX;
								if(minX < MaxX){
									for(int xLine = minX; xLine < MaxX + 1;xLine++){
										//first time
										if(xLine == minX && spanX.current[0] == xLine){
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0]+1 != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								/*else{
											float  intZ = interpZ(MaxX,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,MaxX,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,MaxX,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
								}*/

							}
							else if(yLine <= midY){
								dY = 1;
								Lverts[0]->current[0] = Lverts[0]->current[0]+(Lverts[0]->Xslope*dY);
								Lverts[0]->current[1] = Lverts[0]->current[1]+dY;
								Lverts[0]->current[2] = Lverts[0]->current[2]+(Lverts[0]->Zslope*dY);
								Rverts[0]->current[0] = Rverts[0]->current[0]+(Rverts[0]->Xslope*dY);
								Rverts[0]->current[1] = Rverts[0]->current[1]+dY;
								Rverts[0]->current[2] = Rverts[0]->current[2]+(Rverts[0]->Zslope*dY);

								//set up span DDA
								minX = Lverts[0]->current[0];
								MaxX = Rverts[0]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[0]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[0]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[0]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);
								float dX;
								if(minX < MaxX){
									for(int xLine = minX; xLine < MaxX + 1;xLine++){
										//first time
										if(xLine == minX && spanX.current[0] == xLine){
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0]+1 != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								/*else{
											float  intZ = interpZ(MaxX,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,MaxX,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,MaxX,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
								}*/
							}
							else if(minY == midY && yLine == minY+1){
								float dYL = midY - Lverts[1]->start[1]+1;
								float dYR = midY - Rverts[0]->start[1]+1;
								Lverts[1]->current[0] = Lverts[1]->current[0]+(Lverts[1]->Xslope*dYL);
								Lverts[1]->current[1] = Lverts[1]->current[1]+dYL;
								Lverts[1]->current[2] = Lverts[1]->current[2]+(Lverts[1]->Zslope*dYL);
								Rverts[0]->current[0] = Rverts[0]->current[0]+(Rverts[0]->Xslope*dYR);
								Rverts[0]->current[1] = Rverts[0]->current[1]+dYR;
								Rverts[0]->current[2] = Rverts[0]->current[2]+(Rverts[0]->Zslope*dYR);

								//set up span DDA
								minX = Lverts[1]->current[0];
								MaxX = Rverts[0]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[1]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[0]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[1]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);
								float dX;
								if(minX < MaxX){
									for(int xLine = minX; xLine < MaxX + 1;xLine++){
										//first time
										if(xLine == minX && spanX.current[0] == xLine){
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0] != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								/*else{
											float  intZ = interpZ(MaxX,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,MaxX,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,MaxX,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
								}*/
							}
							else if((yLine == midY+1 && midY != MaxY)){
								
								dY = midY - Lverts[1]->start[1]+1;
								Lverts[1]->current[0] = Lverts[1]->current[0]+(Lverts[1]->Xslope*dY);
								Lverts[1]->current[1] = Lverts[1]->current[1]+dY;
								Lverts[1]->current[2] = Lverts[1]->current[2]+(Lverts[1]->Zslope*dY);
								Rverts[0]->current[0] = Rverts[0]->current[0]+(Rverts[0]->Xslope*1);
								Rverts[0]->current[1] = Rverts[0]->current[1]+1;
								Rverts[0]->current[2] = Rverts[0]->current[2]+(Rverts[0]->Zslope*1);

								//set up span DDA
								minX = Lverts[1]->current[0];
								MaxX = Rverts[0]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[1]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[0]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[1]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);
								float dX;
								if(minX < MaxX){
									for(int xLine = minX; xLine < MaxX + 1;xLine++){
										//first time
										if(xLine == minX && spanX.current[0] == xLine){
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0] != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								/*else{
											float  intZ = interpZ(MaxX,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,MaxX,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,MaxX,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
								}*/
							}
							else if(midY != MaxY && yLine > midY){
								dY = 1;
								Lverts[1]->current[0] = Lverts[1]->current[0]+(Lverts[1]->Xslope*dY);
								Lverts[1]->current[1] = Lverts[1]->current[1]+dY;
								Lverts[1]->current[2] = Lverts[1]->current[2]+(Lverts[1]->Zslope*dY);
								Rverts[0]->current[0] = Rverts[0]->current[0]+(Rverts[0]->Xslope*dY);
								Rverts[0]->current[1] = Rverts[0]->current[1]+dY;
								Rverts[0]->current[2] = Rverts[0]->current[2]+(Rverts[0]->Zslope*dY);

								//set up span DDA
								minX = Lverts[1]->current[0];
								MaxX = Rverts[0]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[1]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[0]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[1]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);
								float dX;
								if(minX < MaxX){
									for(int xLine = minX; xLine < MaxX+1;xLine++){
										//first time
										if(xLine == minX && spanX.current[0] == xLine){
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0] != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);
											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								/*else{
											float  intZ = interpZ(MaxX,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,MaxX,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												r = ctoi(render->flatcolor[0]);
												g = ctoi(render->flatcolor[1]);
												b = ctoi(render->flatcolor[2]);
												a = 4095;
												z = intZ;
												GzPutDisplay(render->display,MaxX,yLine,r,g,b,a,z);
											}
											else{
												//skip the pixel
											}
								}*/
							}

						}
					}





					delete Lverts[0];
					delete Lverts[1];
					delete Rverts[0];
					delete Rverts[1];
					break;
				}
			case(GZ_NULL_TOKEN):
				{
					break;
				}
			}
		}
	}
	else{
		return GZ_FAILURE;
	}

	return GZ_SUCCESS;
}


