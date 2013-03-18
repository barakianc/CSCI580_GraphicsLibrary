/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#define _USE_MATH_DEFINES
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

//DDA struct for edges
struct DDA{
	GzCoord start;
	GzCoord end;
	GzCoord current;
	GzColor	vert1Color;
	GzColor currentColor;
	GzColor vert2Color;
	GzCoord	vert1Norm;
	GzCoord currentNorm;
	GzCoord	vert2Norm;
	float Xslope;
	float Zslope;
	float Rslope;
	float Gslope;
	float Bslope;
	float NXslope;
	float NYslope;
	float NZslope;
	float Ustart,Uend,Vstart,Vend,Ucurrent,Vcurrent;
	float Uslope,Vslope;
};

//span DDA struct
struct SPANDDA{
	GzCoord start;
	GzCoord end;
	GzCoord current;
	GzColor	vert1Color;
	GzColor currentColor;
	GzColor vert2Color;
	GzCoord	vert1Norm;
	GzCoord currentNorm;
	GzCoord	vert2Norm;
	float Zslope;
	float Rslope;
	float Gslope;
	float Bslope;
	float NXslope;
	float NYslope;
	float NZslope;
	float Ustart,Uend,Vstart,Vend,Ucurrent,Vcurrent;
	float Uslope,Vslope;
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

//calculate the color given the normal and the lights
bool calculateColor(GzRender* render,GzCoord norm,GzColor resColor){
	//First we need to calulate N and R, we already know L and E
	//each light has an L as its direction, and E is (0,0,-1) in image space
	GzCoord resN, resN2;
	GzMatrix res1, res2;
	//calulate the Xn transformation
	multMat4x4(render->Xnorm[render->matlevel+1],render->Xnorm[render->matlevel],res1);
	memcpy(res2,res1,16*sizeof(float));
	for(int i = render->matlevel+2; i < MATLEVELS; i++){
		multMat4x4(render->Xnorm[i],res2,res1);
		memcpy(res2,res1,16*sizeof(float));
	}
	float w;

	

	//apply the transformation to the normal
	multMat4x1(res2,norm,resN,&w);
	resN[0] = resN[0]/w;
	resN[1] = resN[1]/w;
	resN[2] = resN[2]/w;

	//now normalize the normal
	float normN = sqrt((resN[0]*resN[0])+(resN[1]*resN[1])+(resN[2]*resN[2]));
	resN[0] = resN[0]/normN;
	resN[1] = resN[1]/normN;
	resN[2] = resN[2]/normN;

	//Now we Loop through each light calculating the R for each light
	GzColor difClr,specClr,ambClr;
	difClr[0] = 0;
	difClr[1] = 0;
	difClr[2] = 0;
	specClr[0] = 0;
	specClr[1] = 0;
	specClr[2] = 0;
	//calculate ambient color contribution
	ambClr[0] = render->ambientlight.color[0]*render->Ka[0];
	ambClr[1] = render->ambientlight.color[1]*render->Ka[1];
	ambClr[2] = render->ambientlight.color[2]*render->Ka[2];

	for(int i = 0; i < render->numlights; i++){
		//First compare N.L and N.E
		float NL = (resN[0]*render->lights[i].direction[0])+(resN[1]*render->lights[i].direction[1])+(resN[2]*render->lights[i].direction[2]);
		float NE = (resN[2]*-1);

		if( NL < 0 && NE < 0){
			resN[0] *= -1;
			resN[1] *= -1;
			resN[2] *= -1;
		}

		if((NL < 0 && NE > 0)|| (NL > 0 && NE < 0)){
			//skip This light
		}
		else{
			//calculate R for the given color
			NL = (resN[0]*render->lights[i].direction[0])+(resN[1]*render->lights[i].direction[1])+(resN[2]*render->lights[i].direction[2]);
			GzCoord R;
			R[0] = (2*NL*resN[0])-render->lights[i].direction[0];
			R[1] = (2*NL*resN[1])-render->lights[i].direction[1];
			R[2] = (2*NL*resN[2])-render->lights[i].direction[2];
			//normalize R
			float normR = sqrt((R[0]*R[0])+(R[1]*R[1])+(R[2]*R[2]));
			R[0] /= normR;
			R[1] /= normR;
			R[2] /= normR;

			//added in the diffuse component
			difClr[0] += (render->lights[i].color[0] *NL);
			difClr[1] += (render->lights[i].color[1] *NL);
			difClr[2] += (render->lights[i].color[2] *NL);

			//calculate R.E and add the spec component
			float RE = (R[2]*-1);
			//clamp R.E to [0,1]
			if(RE > 1.0){
				RE = 1.0;
			}
			else if(RE < 0.0){
				RE = 0.0;
			}
			//calc spec color
			specClr[0] += (render->lights[i].color[0]*(pow(RE,render->spec)));
			specClr[1] += (render->lights[i].color[1]*(pow(RE,render->spec)));
			specClr[2] += (render->lights[i].color[2]*(pow(RE,render->spec)));
		}
		

	}
	//now multiply the spec and diffuse sums by there respective components
	difClr[0] *= render->Kd[0];
	difClr[1] *= render->Kd[1];
	difClr[2] *= render->Kd[2];
	specClr[0] *= render->Ks[0];
	specClr[1] *= render->Ks[1];
	specClr[2] *= render->Ks[2];

	//now sum all the color components
	resColor[0] = specClr[0] + difClr[0] + ambClr[0];
	resColor[1] = specClr[1] + difClr[1] + ambClr[1];
	resColor[2] = specClr[2] + difClr[2] + ambClr[2];

	//check for color overflow and clamp color to 1.0
	if(resColor[0] > 1.0){
		resColor[0] = 1.0;
	}

	if(resColor[1] > 1.0){
		resColor[1] = 1.0;
	}

	if(resColor[2] > 1.0){
		resColor[2] = 1.0;
	}


	return true;
}

//calculate color for  Gouraud Shading not applying Kt and just returning the sums
bool calculateColorGouraud(GzRender* render,GzCoord norm,GzColor resColor){
	//First we need to calulate N and R, we already know L and E
	//each light has an L as its direction, and E is (0,0,-1) in image space
	GzCoord resN, resN2;
	GzMatrix res1, res2;
	//calulate the Xn transformation
	multMat4x4(render->Xnorm[render->matlevel+1],render->Xnorm[render->matlevel],res1);
	memcpy(res2,res1,16*sizeof(float));
	for(int i = render->matlevel+2; i < MATLEVELS; i++){
		multMat4x4(render->Xnorm[i],res2,res1);
		memcpy(res2,res1,16*sizeof(float));
	}
	float w;

	

	//apply the transformation to the normal
	multMat4x1(res2,norm,resN,&w);
	resN[0] = resN[0]/w;
	resN[1] = resN[1]/w;
	resN[2] = resN[2]/w;

	//now normalize the normal
	float normN = sqrt((resN[0]*resN[0])+(resN[1]*resN[1])+(resN[2]*resN[2]));
	resN[0] = resN[0]/normN;
	resN[1] = resN[1]/normN;
	resN[2] = resN[2]/normN;

	//Now we Loop through each light calculating the R for each light
	GzColor difClr,specClr,ambClr;
	difClr[0] = 0;
	difClr[1] = 0;
	difClr[2] = 0;
	specClr[0] = 0;
	specClr[1] = 0;
	specClr[2] = 0;
	//calculate ambient color contribution
	ambClr[0] = render->ambientlight.color[0];
	ambClr[1] = render->ambientlight.color[1];
	ambClr[2] = render->ambientlight.color[2];

	for(int i = 0; i < render->numlights; i++){
		//First compare N.L and N.E
		float NL = (resN[0]*render->lights[i].direction[0])+(resN[1]*render->lights[i].direction[1])+(resN[2]*render->lights[i].direction[2]);
		float NE = (resN[2]*-1);

		if( NL < 0 && NE < 0){
			resN[0] *= -1;
			resN[1] *= -1;
			resN[2] *= -1;
		}

		if((NL < 0 && NE > 0)|| (NL > 0 && NE < 0)){
			//skip This light
		}
		else{
			//calculate R for the given color
			NL = (resN[0]*render->lights[i].direction[0])+(resN[1]*render->lights[i].direction[1])+(resN[2]*render->lights[i].direction[2]);
			GzCoord R;
			R[0] = (2*NL*resN[0])-render->lights[i].direction[0];
			R[1] = (2*NL*resN[1])-render->lights[i].direction[1];
			R[2] = (2*NL*resN[2])-render->lights[i].direction[2];
			//normalize R
			float normR = sqrt((R[0]*R[0])+(R[1]*R[1])+(R[2]*R[2]));
			R[0] /= normR;
			R[1] /= normR;
			R[2] /= normR;

			//added in the diffuse component
			difClr[0] += (render->lights[i].color[0] *NL);
			difClr[1] += (render->lights[i].color[1] *NL);
			difClr[2] += (render->lights[i].color[2] *NL);

			//calculate R.E and add the spec component
			float RE = (R[2]*-1);
			//clamp R.E to [0,1]
			if(RE > 1.0){
				RE = 1.0;
			}
			else if(RE < 0.0){
				RE = 0.0;
			}
			//calc spec color
			specClr[0] += (render->lights[i].color[0]*(pow(RE,render->spec)));
			specClr[1] += (render->lights[i].color[1]*(pow(RE,render->spec)));
			specClr[2] += (render->lights[i].color[2]*(pow(RE,render->spec)));
		}
		

	}
	//now multiply the spec and diffuse sums by there respective components
	//we dont do this for Gouraud Shading
	/*
	difClr[0] *= render->Kd[0];
	difClr[1] *= render->Kd[1];
	difClr[2] *= render->Kd[2];
	specClr[0] *= render->Ks[0];
	specClr[1] *= render->Ks[1];
	specClr[2] *= render->Ks[2];*/

	//now sum all the color components
	resColor[0] = specClr[0] + difClr[0] + ambClr[0];
	resColor[1] = specClr[1] + difClr[1] + ambClr[1];
	resColor[2] = specClr[2] + difClr[2] + ambClr[2];

	//check for color overflow and clamp color to 1.0
	if(resColor[0] > 1.0){
		resColor[0] = 1.0;
	}

	if(resColor[1] > 1.0){
		resColor[1] = 1.0;
	}

	if(resColor[2] > 1.0){
		resColor[2] = 1.0;
	}


	return true;
}


//Matrix multiplier for a 4x4 and 4x4 matrix returns teh resulting matrix
bool multMat4x4(GzMatrix m1,GzMatrix m2,GzMatrix result){
	for(int i = 0; i < 4; i++){
		result[i][0] = (m1[i][0]*m2[0][0])+(m1[i][1]*m2[1][0])+(m1[i][2]*m2[2][0])+(m1[i][3]*m2[3][0]);
		result[i][1] = (m1[i][0]*m2[0][1])+(m1[i][1]*m2[1][1])+(m1[i][2]*m2[2][1])+(m1[i][3]*m2[3][1]);
		result[i][2] = (m1[i][0]*m2[0][2])+(m1[i][1]*m2[1][2])+(m1[i][2]*m2[2][2])+(m1[i][3]*m2[3][2]);
		result[i][3] = (m1[i][0]*m2[0][3])+(m1[i][1]*m2[1][3])+(m1[i][2]*m2[2][3])+(m1[i][3]*m2[3][3]);

	}
	return true;
}

//Matrix multiplier for a 4x4 and 4x1
bool multMat4x1(GzMatrix m1,GzCoord c1,GzCoord result,float* w){
	for(int i = 0; i < 3; i++){
		result[i] = (m1[i][0]*c1[0])+(m1[i][1]*c1[1])+(m1[i][2]*c1[2])+(m1[i][3]*1);
	}
	*w = (m1[3][0]*c1[0])+(m1[3][1]*c1[1])+(m1[3][2]*c1[2])+(m1[3][3]*1);
	/*result[0] = result[0]/w;
	result[1] = result[1]/w;
	result[2] = result[2]/w;*/
	return true;
}


int GzRotXMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value
	GzMatrix tempmat = 
	{
		1.0,	0.0,					0.0,						0.0,
		0.0,	cos(degree*(M_PI/180)),	-sin(degree*(M_PI/180)),	0.0,
		0.0,	sin(degree*(M_PI/180)),	cos(degree*(M_PI/180)),		0.0,
		0.0,	0.0,					0.0,						1.0
	};
	memcpy(mat,tempmat,16*sizeof(float));
	return GZ_SUCCESS;
}


int GzRotYMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value
	GzMatrix tempmat = 
	{
		cos(degree*(M_PI/180)),		0.0,	sin(degree*(M_PI/180)),	0.0,
		0.0,						1.0,	0.0,					0.0,
		-sin(degree*(M_PI/180)),	0.0,	cos(degree*(M_PI/180)),	0.0,
		0.0,						0.0,	0.0,					1.0
	};
	memcpy(mat,tempmat,16*sizeof(float));

	return GZ_SUCCESS;
}


int GzRotZMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value

	GzMatrix tempmat = 
	{
		cos(degree*(M_PI/180)),	-sin(degree*(M_PI/180)),	0.0,	0.0,
		sin(degree*(M_PI/180)),	cos(degree*(M_PI/180)),		0.0,	0.0,
		0.0,					0.0,						1.0,	0.0,
		0.0,					0.0,						0.0,	1.0
	};
	memcpy(mat,tempmat,16*sizeof(float));
	return GZ_SUCCESS;
}


int GzTrxMat(GzCoord translate, GzMatrix mat)
{
// Create translation matrix
// Pass back the matrix using mat value
	GzMatrix tempmat = 
	{
		1.0,	0.0,	0.0,	translate[0],
		0.0,	1.0,	0.0,	translate[1],
		0.0,	0.0,	1.0,	translate[2],
		0.0,	0.0,	0.0,	1.0
	};
	memcpy(mat,tempmat,16*sizeof(float));
	return GZ_SUCCESS;
}


int GzScaleMat(GzCoord scale, GzMatrix mat)
{
// Create scaling matrix
// Pass back the matrix using mat value
	GzMatrix tempmat = 
	{
		scale[0],	0.0,		0.0,		0.0,
		0.0,		scale[1],	0.0,		0.0,
		0.0,		0.0,		scale[2],	0.0,
		0.0,		0.0,		0.0,		1.0
	};
	memcpy(mat,tempmat,16*sizeof(float));
	return GZ_SUCCESS;
}


//----------------------------------------------------------
// Begin main functions

int GzNewRender(GzRender **render, GzRenderClass renderClass, GzDisplay	*display)
{
/*  
- malloc a renderer struct 
- keep closed until all inits are done 
- setup Xsp and anything only done once 
- span interpolator needs pointer to display 
- check for legal class GZ_Z_BUFFER_RENDER 
- init default camera 
*/ 

	if(display != NULL && renderClass == GZ_Z_BUFFER_RENDER ){
		*render = (GzRender*)malloc(sizeof(GzRender));
		(*render)->renderClass = renderClass;
		(*render)->display = display;
		(*render)->open = 0;
		(*render)->matlevel = MATLEVELS;
		//init the default camera
		(*render)->camera.position[0] = DEFAULT_IM_X;
		(*render)->camera.position[1] = DEFAULT_IM_Y;
		(*render)->camera.position[2] = DEFAULT_IM_Z;
		(*render)->camera.lookat[0] = 0.0;
		(*render)->camera.lookat[1] = 0.0;
		(*render)->camera.lookat[2] = 0.0;
		(*render)->camera.worldup[0] = 0.0;
		(*render)->camera.worldup[1] = 1.0;
		(*render)->camera.worldup[2] = 0.0;
		(*render)->camera.FOV = DEFAULT_FOV;
		
		//initialize Xsp
		GzMatrix tempXsp = 
		{
			(display->xres/2),	0.0,					0.0,										(display->xres/2),
			0.0,				(-display->yres/2),		0.0,										(display->yres/2),
			0.0,				0.0,					(MAXINT*(tan((((*render)->camera.FOV)/2)*(M_PI/180)))),	0.0,
			0.0,				0.0,					0.0,										1.0
		};
		memcpy((*render)->Xsp,tempXsp,16*sizeof(float));
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


int GzBeginRender(GzRender *render)
{
/*  
- set up for start of each frame - clear frame buffer 
- compute Xiw and projection xform Xpi from camera definition 
- init Ximage - put Xsp at base of stack, push on Xpi and Xiw 
- now stack contains Xsw and app can push model Xforms if it want to. 
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

		//Push Xsp at base of the stack
		GzPushMatrix(render,render->Xsp);
		
		//Calculate Xpi for this render frame
		GzMatrix tempXpi =
		{
			1.0,	0.0,	0.0,							0.0,
			0.0,	1.0,	0.0,							0.0,
			0.0,	0.0,	1.0,							0.0,
			0.0,	0.0,	(tan(((render->camera.FOV)/2)*(M_PI/180))),	1.0
		};

		GzMatrix res;
		//calculate Xsi
		multMat4x4(render->Xsp,tempXpi,res);

		//pop Xsp off
		GzPopMatrix(render);
		//Push Xsi for this render frame
		GzPushMatrix(render,res);

		//calculate Xiw for this render frame
		GzCoord cl =
		{
			(render->camera.lookat[0]-render->camera.position[0]),	(render->camera.lookat[1]-render->camera.position[1]),	(render->camera.lookat[2]-render->camera.position[2])
		};

		float normcl = sqrt((cl[0]*cl[0])+(cl[1]*cl[1])+(cl[2]*cl[2]));

		GzCoord z =
		{
			(cl[0]/normcl),	(cl[1]/normcl),	(cl[2]/normcl)
		};

		float UpdotZ = (render->camera.worldup[0]*z[0]) + (render->camera.worldup[1]*z[1]) + (render->camera.worldup[2]*z[2]);

		GzCoord UpPrime =
		{
			(render->camera.worldup[0]-(UpdotZ*z[0])),	(render->camera.worldup[1]-(UpdotZ*z[1])), (render->camera.worldup[2]-(UpdotZ*z[2]))
		};

		float normUp = sqrt((UpPrime[0]*UpPrime[0])+(UpPrime[1]*UpPrime[1])+(UpPrime[2]*UpPrime[2]));

		GzCoord y =
		{
			(UpPrime[0]/normUp),	(UpPrime[1]/normUp),	(UpPrime[2]/normUp)
		};

		GzCoord x =
		{
			((y[1]*z[2])-(y[2]*z[1])),	((y[2]*z[0])-(y[0]*z[2])),	((y[0]*z[1])-(y[1]*z[0]))
		};

		float XdC = (x[0]*render->camera.position[0])+(x[1]*render->camera.position[1])+(x[2]*render->camera.position[2]);
		float YdC = (y[0]*render->camera.position[0])+(y[1]*render->camera.position[1])+(y[2]*render->camera.position[2]);
		float ZdC = (z[0]*render->camera.position[0])+(z[1]*render->camera.position[1])+(z[2]*render->camera.position[2]);

		GzMatrix tempXiw =
		{
			x[0],	x[1],	x[2],	-XdC,
			y[0],	y[1],	y[2],	-YdC,
			z[0],	z[1],	z[2],	-ZdC,
			0.0,	0.0,	0.0,	1.0
		};

		//calculate Xsw
		multMat4x4(render->Ximage[render->matlevel],tempXiw,res);
		//pop off Xsi
		GzPopMatrix(render);
		//Push Xiw onto the stack
		GzPushMatrix(render,res);
		GzMatrix tempMat;
		memcpy(tempMat,tempXiw,16*sizeof(float));
		tempMat[0][3] = 0;
		tempMat[1][3] = 0;
		tempMat[2][3] = 0;

		//compute k
		float k = 1 / (sqrt((pow(tempMat[0][0],2))+(pow(tempMat[0][1],2))+(pow(tempMat[0][2],2))));

		//apply to the 3x3 matrix
		for(int i = 0; i < 3; i++){
			tempMat[i][0] = tempMat[i][0]*k;
			tempMat[i][1] = tempMat[i][1]*k;
			tempMat[i][2] = tempMat[i][2]*k;
		}
		memcpy(render->Xnorm[99],tempMat,sizeof(float)*16);

		//init the flat render color to blue
		render->flatcolor[0] = 0;
		render->flatcolor[1] = 0;
		render->flatcolor[2] = 1.0;
		
		render->interp_mode = 0;
		render->numlights = 0;

		//init Ka to 0
		render->Ka[0] = 0.1;
		render->Ka[1] = 0.1;
		render->Ka[2] = 0.1;

		//init Kd to 0
		render->Kd[0] = 0.7;
		render->Kd[1] = 0.6;
		render->Kd[2] = 0.5;

		//init Ks to 0
		render->Ks[0] = 0.2;
		render->Ks[1] = 0.3;
		render->Ks[2] = 0.4;

		//init spec to 0
		render->spec = DEFAULT_SPEC;

		//init aa offset
		render->aaX = 0.0;
		render->aaY = 0.0;
	}
	else{
		return GZ_FAILURE;
	}
	return GZ_SUCCESS;
}

int GzPutCamera(GzRender *render, GzCamera *camera)
{
/*
- overwrite renderer camera structure with new camera definition
*/
	if(render != NULL && camera != NULL){
		memcpy(&(render->camera),camera,sizeof(GzCamera));

		//create new Xsp
		GzMatrix tempXsp =
		{
			((render->display->xres)/2),		0.0,							0.0,									((render->display->xres)/2),
			0.0,								(-(render->display->yres)/2),	0.0,									((render->display->yres)/2),
			0.0,								0.0,							(MAXINT*(tan((render->camera.FOV*(M_PI/180))/2))),	0.0,
			0.0,								0.0,							0.0,									1.0
		};

		memcpy(render->Xsp,tempXsp,16*sizeof(float));
	}
	else{
		return GZ_FAILURE;
	}
	return GZ_SUCCESS;	
}

int GzPushMatrix(GzRender *render, GzMatrix	matrix)
{
/*
- push a matrix onto the Ximage stack
- check for stack overflow
*/
	if(render != NULL){
		if(render->matlevel > 0 && render->matlevel < MATLEVELS){
			//compute the new matrix to put onto Xnorm
			GzMatrix tempMat;
			memcpy(tempMat,matrix,16*sizeof(float));
			//make translation zeroed out
			tempMat[0][3] = 0;
			tempMat[1][3] = 0;
			tempMat[2][3] = 0;

			//compute k
			float k = 1 / (sqrt((pow(tempMat[0][0],2))+(pow(tempMat[0][1],2))+(pow(tempMat[0][2],2))));

			//apply to the 3x3 matrix
			for(int i = 0; i < 3; i++){
				tempMat[i][0] = tempMat[i][0]*k;
				tempMat[i][1] = tempMat[i][1]*k;
				tempMat[i][2] = tempMat[i][2]*k;
			}

			//now push the matrix onto Xn
			memcpy(render->Xnorm[render->matlevel-1],tempMat,16*sizeof(float));
		}

		if(render->matlevel > 0){
			render->matlevel--;
			memcpy(render->Ximage[render->matlevel],matrix,16*sizeof(float));
		}

		
	}
	else{
		return GZ_FAILURE;
	}
	return GZ_SUCCESS;
}

int GzPopMatrix(GzRender *render)
{
/*
- pop a matrix off the Ximage stack
- check for stack underflow
*/
	if(render != NULL){
		if(render->matlevel < MATLEVELS){
			render->matlevel++;
		}
	}
	else{
		return GZ_FAILURE;
	}
	return GZ_SUCCESS;
}


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer	*valueList) /* void** valuelist */
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
			case(GZ_INTERPOLATE):
				{
					//get the interpolation mode
					render->interp_mode = *(int*)valueList[i];
					break;
				}
			case(GZ_DIRECTIONAL_LIGHT):
				{
					if(render->numlights < MAX_LIGHTS){
						memcpy(render->lights[render->numlights].color,((GzLight*)valueList[i])->color,sizeof(float)*3);
						memcpy(render->lights[render->numlights].direction,((GzLight*)valueList[i])->direction,sizeof(float)*3);
						render->numlights++;
					}
					break;
				}
			case(GZ_AMBIENT_LIGHT):
				{
					memcpy(&render->ambientlight,(GzLight*)valueList[i],sizeof(GzLight));
					break;
				}
			case(GZ_AMBIENT_COEFFICIENT):
				{
					memcpy(render->Ka,(GzLight*)valueList[i],sizeof(float)*3);
					break;
				}
			case(GZ_DIFFUSE_COEFFICIENT):
				{
					memcpy(render->Kd,(GzLight*)valueList[i],sizeof(float)*3);
					break;
				}
			case(GZ_SPECULAR_COEFFICIENT):
				{
					memcpy(render->Ks,(GzLight*)valueList[i],sizeof(float)*3);
					break;
				}
			case(GZ_DISTRIBUTION_COEFFICIENT):
				{
					render->spec = *(float*)valueList[i];
					break;
				}
			case(GZ_TEXTURE_MAP):
				{
					render->tex_fun = (GzTexture)valueList[i];
					break;
				}
			case(GZ_AASHIFTX):
				{
					render->aaX = *(float*)valueList[i];
					break;
				}
			case(GZ_AASHIFTY):
				{
					render->aaY = *(float*)valueList[i];
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

int GzPutTriangle(GzRender	*render, int numParts, GzToken *nameList, 
				  GzPointer	*valueList)
/* numParts : how many names and values */
{
/*  
- pass in a triangle description with tokens and values corresponding to 
      GZ_POSITION:3 vert positions in model space 
- Xform positions of verts  
- Clip - just discard any triangle with verts behind view plane 
       - test for triangles with all three verts off-screen 
- invoke triangle rasterizer  
*/ 
	if(render != NULL && numParts >= 0 && nameList != NULL && valueList != NULL && render->open == 1){
		for(int i = 0; i< numParts; i++){
			switch(nameList[i]){
			case(GZ_POSITION):
				{
					//for each triangle we to determine the pixel that make up that triangle and set the color accordingly
					//We also need to interpolate z for each pixel to do z-buffer checking

					GzPointer tempPointer = valueList[i];
					GzPointer tempPointer2 = valueList[i+1];
					GzPointer tempPointer3 = valueList[i+2];
					//get the coords for the 3 verts and the 3 normals
					GzCoord V1 = {((GzCoord*)tempPointer)[0][0],((GzCoord*)tempPointer)[0][1],((GzCoord*)tempPointer)[0][2]};
					GzCoord V2 = {((GzCoord*)tempPointer)[1][0],((GzCoord*)tempPointer)[1][1],((GzCoord*)tempPointer)[1][2]};
					GzCoord V3 = {((GzCoord*)tempPointer)[2][0],((GzCoord*)tempPointer)[2][1],((GzCoord*)tempPointer)[2][2]};
					GzCoord N1 = {((GzCoord*)tempPointer2)[0][0],((GzCoord*)tempPointer2)[0][1],((GzCoord*)tempPointer2)[0][2]};
					GzCoord N2 = {((GzCoord*)tempPointer2)[1][0],((GzCoord*)tempPointer2)[1][1],((GzCoord*)tempPointer2)[1][2]};
					GzCoord N3 = {((GzCoord*)tempPointer2)[2][0],((GzCoord*)tempPointer2)[2][1],((GzCoord*)tempPointer2)[2][2]};
					GzTextureIndex uv1 ={((GzTextureIndex*)tempPointer3)[0][0],((GzTextureIndex*)tempPointer3)[0][1]};
					GzTextureIndex uv2 ={((GzTextureIndex*)tempPointer3)[1][0],((GzTextureIndex*)tempPointer3)[1][1]};
					GzTextureIndex uv3 ={((GzTextureIndex*)tempPointer3)[2][0],((GzTextureIndex*)tempPointer3)[2][1]};

					GzCoord resV1, resV2, resV3;
					GzMatrix res1, res2;
					multMat4x4(render->Ximage[render->matlevel+1],render->Ximage[render->matlevel],res1);
					memcpy(res2,res1,16*sizeof(float));
					for(int i = render->matlevel+2; i < MATLEVELS; i++){
						multMat4x4(render->Ximage[i],res2,res1);
						memcpy(res2,res1,16*sizeof(float));
					}

					float w1 = 0.0f,w2 = 0.0f,w3 = 0.0f;
					//transform points
					multMat4x1(res2,V1,resV1,&w1);
					memcpy(V1,resV1,3*sizeof(float));
					multMat4x1(res2,V2,resV2,&w2);
					memcpy(V2,resV2,3*sizeof(float));
					multMat4x1(res2,V3,resV3,&w3);
					memcpy(V3,resV3,3*sizeof(float));

					//we need to discard all triangles that are behind the camera
					if(V1[2] < 0 && V2[2] < 0 && V3[2] < 0){
						return GZ_SUCCESS;
					}
					/*Now devide by each vertexs repective w*/
					V1[0] = V1[0]/w1;
					V1[1] = V1[1]/w1;
					V1[2] = V1[2]/w1;

					V2[0] = V2[0]/w2;
					V2[1] = V2[1]/w2;
					V2[2] = V2[2]/w2;

					V3[0] = V3[0]/w3;
					V3[1] = V3[1]/w3;
					V3[2] = V3[2]/w3;


					//sort verts by Y coord, least to greatest
					GzCoord vertsList[3];
					GzCoord normList[3];
					GzTextureIndex uvList[3];
					//verts
					memcpy(vertsList[0],V1,3*sizeof(float));
					memcpy(vertsList[1],V2,3*sizeof(float));
					memcpy(vertsList[2],V3,3*sizeof(float));
					//normals
					memcpy(normList[0],N1,3*sizeof(float));
					memcpy(normList[1],N2,3*sizeof(float));
					memcpy(normList[2],N3,3*sizeof(float));
					//uv's
					memcpy(uvList[0],uv1,2*sizeof(float));
					memcpy(uvList[1],uv2,2*sizeof(float));
					memcpy(uvList[2],uv3,2*sizeof(float));

					bool changes = true;
					GzCoord tempCoord;
					GzTextureIndex tempUV;
					while(changes){
						changes = false;
						for(int j = 0; j < 2 ; j++){
							if(vertsList[j][1] > vertsList[j+1][1]){
								//verts
								memcpy(tempCoord,vertsList[j+1],3*sizeof(float));
								memcpy(vertsList[j+1],vertsList[j],3*sizeof(float));
								memcpy(vertsList[j],tempCoord,3*sizeof(float));
								//norms
								memcpy(tempCoord,normList[j+1],3*sizeof(float));
								memcpy(normList[j+1],normList[j],3*sizeof(float));
								memcpy(normList[j],tempCoord,3*sizeof(float));
								//uv's
								memcpy(tempUV,uvList[j+1],2*sizeof(float));
								memcpy(uvList[j+1],uvList[j],2*sizeof(float));
								memcpy(uvList[j],tempUV,2*sizeof(float));
								changes = true;
							}
						}
					}

					/*apply aa offset*/
					vertsList[0][0] -= render->aaX;
					vertsList[0][1] -= render->aaY;
					vertsList[1][0] -= render->aaX;
					vertsList[1][1] -= render->aaY;
					vertsList[2][0] -= render->aaX;
					vertsList[2][1] -= render->aaY;
					/*.....................................*/

					/*if(abs(vertsList[2][1] - vertsList[1][1]) < .001){
						int test = 0;
					}*/
					DDA* Lverts[2];
					DDA* Rverts[2];
					//if the triangle is L(1-3) R(1-2 2-3)
					float s1 = (vertsList[1][0]-vertsList[0][0])/(vertsList[1][1]-vertsList[0][1]);
					float s2 = (vertsList[2][0]-vertsList[0][0])/(vertsList[2][1]-vertsList[0][1]);
					/*if((((vertsList[1][0]-vertsList[0][0])/(vertsList[1][1]-vertsList[0][1]) < (vertsList[2][0]-vertsList[0][0])/(vertsList[2][1]-vertsList[0][1])) && vertsList[1][0] < vertsList[2][0]) && (vertsList[1][0]-vertsList[0][0])/(vertsList[1][1]-vertsList[0][1])> 0 && (vertsList[2][0]-vertsList[0][0])/(vertsList[2][1]-vertsList[0][1]) > 0){
						int test = 0;
					}*/
					if(vertsList[1][1] == vertsList[2][1]){
						//int test = 0;
						//L edge DDA
						Lverts[0] = new DDA();
						memcpy(Lverts[0]->start,vertsList[0],3*sizeof(float));
						memcpy(Lverts[0]->end,vertsList[2],3*sizeof(float));
						memcpy(Lverts[0]->current,vertsList[0],3*sizeof(float));
						Lverts[0]->Xslope = (Lverts[0]->end[0]-Lverts[0]->start[0])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[0]->Zslope = (Lverts[0]->end[2]-Lverts[0]->start[2])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						memcpy(Lverts[0]->vert1Norm,normList[0],3*sizeof(float));
						memcpy(Lverts[0]->currentNorm,normList[0],3*sizeof(float));
						memcpy(Lverts[0]->vert2Norm,normList[2],3*sizeof(float));
						Lverts[0]->NXslope = (Lverts[0]->vert2Norm[0]-Lverts[0]->vert1Norm[0])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[0]->NYslope = (Lverts[0]->vert2Norm[1]-Lverts[0]->vert1Norm[1])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[0]->NZslope = (Lverts[0]->vert2Norm[2]-Lverts[0]->vert1Norm[2])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[0]->Ustart = uvList[0][0]/((Lverts[0]->start[2]/(MAXINT-Lverts[0]->start[2]))+1);
						Lverts[0]->Uend = uvList[2][0]/((Lverts[0]->end[2]/(MAXINT-Lverts[0]->end[2]))+1);
						Lverts[0]->Ucurrent = uvList[0][0]/((Lverts[0]->start[2]/(MAXINT-Lverts[0]->start[2]))+1);
						Lverts[0]->Vstart = uvList[0][1]/((Lverts[0]->start[2]/(MAXINT-Lverts[0]->start[2]))+1);
						Lverts[0]->Vend = uvList[2][1]/((Lverts[0]->end[2]/(MAXINT-Lverts[0]->end[2]))+1);
						Lverts[0]->Vcurrent = uvList[0][1]/((Lverts[0]->start[2]/(MAXINT-Lverts[0]->start[2]))+1);
						Lverts[0]->Uslope = (Lverts[0]->Uend-Lverts[0]->Ustart)/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[0]->Vslope = (Lverts[0]->Vend-Lverts[0]->Vstart)/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[1] = NULL;

						//R edges DDAs
						Rverts[0] = new DDA();
						memcpy(Rverts[0]->start,vertsList[0],3*sizeof(float));
						memcpy(Rverts[0]->end,vertsList[1],3*sizeof(float));
						memcpy(Rverts[0]->current, vertsList[0],3*sizeof(float));
						Rverts[0]->Xslope = (Rverts[0]->end[0]-Rverts[0]->start[0])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->Zslope = (Rverts[0]->end[2]-Rverts[0]->start[2])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						memcpy(Rverts[0]->vert1Norm,normList[0],3*sizeof(float));
						memcpy(Rverts[0]->currentNorm,normList[0],3*sizeof(float));
						memcpy(Rverts[0]->vert2Norm,normList[1],3*sizeof(float));
						Rverts[0]->NXslope = (Rverts[0]->vert2Norm[0]-Rverts[0]->vert1Norm[0])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->NYslope = (Rverts[0]->vert2Norm[1]-Rverts[0]->vert1Norm[1])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->NZslope = (Rverts[0]->vert2Norm[2]-Rverts[0]->vert1Norm[2])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->Ustart = uvList[0][0]/((Rverts[0]->start[2]/(MAXINT-Rverts[0]->start[2]))+1);
						Rverts[0]->Uend = uvList[1][0]/((Rverts[0]->end[2]/(MAXINT-Rverts[0]->end[2]))+1);
						Rverts[0]->Ucurrent = uvList[0][0]/((Rverts[0]->start[2]/(MAXINT-Rverts[0]->start[2]))+1);
						Rverts[0]->Vstart = uvList[0][1]/((Rverts[0]->start[2]/(MAXINT-Rverts[0]->start[2]))+1);
						Rverts[0]->Vend = uvList[1][1]/((Rverts[0]->end[2]/(MAXINT-Rverts[0]->end[2]))+1);
						Rverts[0]->Vcurrent = uvList[0][1]/((Rverts[0]->start[2]/(MAXINT-Rverts[0]->start[2]))+1);
						Rverts[0]->Uslope = (Rverts[0]->Uend-Rverts[0]->Ustart)/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->Vslope = (Rverts[0]->Vend-Rverts[0]->Vstart)/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[1] = new DDA();
						memcpy(Rverts[1]->start,vertsList[1],3*sizeof(float));
						memcpy(Rverts[1]->end,vertsList[2],3*sizeof(float));
						memcpy(Rverts[1]->current,vertsList[1],3*sizeof(float));
						Rverts[1]->Xslope = (Rverts[1]->end[0]-Rverts[1]->start[0])/(Rverts[1]->end[1]-Rverts[1]->start[1]);
						Rverts[1]->Zslope = (Rverts[1]->end[2]-Rverts[1]->start[2])/(Rverts[1]->end[1]-Rverts[1]->start[1]);
						memcpy(Rverts[1]->vert1Norm,normList[1],3*sizeof(float));
						memcpy(Rverts[1]->currentNorm,normList[1],3*sizeof(float));
						memcpy(Rverts[1]->vert2Norm,normList[2],3*sizeof(float));
						Rverts[1]->NXslope = (Rverts[1]->vert2Norm[0]-Rverts[1]->vert1Norm[0])/(Rverts[1]->end[1]-Rverts[1]->start[1]);
						Rverts[1]->NYslope = (Rverts[1]->vert2Norm[1]-Rverts[1]->vert1Norm[1])/(Rverts[1]->end[1]-Rverts[1]->start[1]);
						Rverts[1]->NZslope = (Rverts[1]->vert2Norm[2]-Rverts[1]->vert1Norm[2])/(Rverts[1]->end[1]-Rverts[1]->start[1]);
						Rverts[1]->Ustart = uvList[1][0]/((Rverts[1]->start[2]/(MAXINT-Rverts[1]->start[2]))+1);
						Rverts[1]->Uend = uvList[2][0]/((Rverts[1]->end[2]/(MAXINT-Rverts[1]->end[2]))+1);
						Rverts[1]->Ucurrent = uvList[1][0]/((Rverts[1]->start[2]/(MAXINT-Rverts[1]->start[2]))+1);
						Rverts[1]->Vstart = uvList[1][1]/((Rverts[1]->start[2]/(MAXINT-Rverts[1]->start[2]))+1);
						Rverts[1]->Vend = uvList[2][1]/((Rverts[1]->end[2]/(MAXINT-Rverts[1]->end[2]))+1);
						Rverts[1]->Vcurrent = uvList[1][1]/((Rverts[1]->start[2]/(MAXINT-Rverts[1]->start[2]))+1);
						Rverts[1]->Uslope = (Rverts[1]->Uend-Rverts[1]->Ustart)/(Rverts[1]->end[1]-Rverts[1]->start[1]);
						Rverts[1]->Vslope = (Rverts[1]->Vend-Rverts[1]->Vstart)/(Rverts[1]->end[1]-Rverts[1]->start[1]);
					}
					else if (((vertsList[1][0]-vertsList[0][0])/(vertsList[1][1]-vertsList[0][1]) > (vertsList[2][0]-vertsList[0][0])/(vertsList[2][1]-vertsList[0][1]))){
						//L edge DDA
						Lverts[0] = new DDA();
						memcpy(Lverts[0]->start,vertsList[0],3*sizeof(float));
						memcpy(Lverts[0]->end,vertsList[2],3*sizeof(float));
						memcpy(Lverts[0]->current,vertsList[0],3*sizeof(float));
						Lverts[0]->Xslope = (Lverts[0]->end[0]-Lverts[0]->start[0])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[0]->Zslope = (Lverts[0]->end[2]-Lverts[0]->start[2])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						memcpy(Lverts[0]->vert1Norm,normList[0],3*sizeof(float));
						memcpy(Lverts[0]->currentNorm,normList[0],3*sizeof(float));
						memcpy(Lverts[0]->vert2Norm,normList[2],3*sizeof(float));
						Lverts[0]->NXslope = (Lverts[0]->vert2Norm[0]-Lverts[0]->vert1Norm[0])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[0]->NYslope = (Lverts[0]->vert2Norm[1]-Lverts[0]->vert1Norm[1])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[0]->NZslope = (Lverts[0]->vert2Norm[2]-Lverts[0]->vert1Norm[2])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[0]->Ustart = uvList[0][0]/((Lverts[0]->start[2]/(MAXINT-Lverts[0]->start[2]))+1);
						Lverts[0]->Uend = uvList[2][0]/((Lverts[0]->end[2]/(MAXINT-Lverts[0]->end[2]))+1);
						Lverts[0]->Ucurrent = uvList[0][0]/((Lverts[0]->start[2]/(MAXINT-Lverts[0]->start[2]))+1);
						Lverts[0]->Vstart = uvList[0][1]/((Lverts[0]->start[2]/(MAXINT-Lverts[0]->start[2]))+1);
						Lverts[0]->Vend = uvList[2][1]/((Lverts[0]->end[2]/(MAXINT-Lverts[0]->end[2]))+1);
						Lverts[0]->Vcurrent = uvList[0][1]/((Lverts[0]->start[2]/(MAXINT-Lverts[0]->start[2]))+1);
						Lverts[0]->Uslope = (Lverts[0]->Uend-Lverts[0]->Ustart)/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[0]->Vslope = (Lverts[0]->Vend-Lverts[0]->Vstart)/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[1] = NULL;

						//R edges DDAs
						Rverts[0] = new DDA();
						memcpy(Rverts[0]->start,vertsList[0],3*sizeof(float));
						memcpy(Rverts[0]->end,vertsList[1],3*sizeof(float));
						memcpy(Rverts[0]->current, vertsList[0],3*sizeof(float));
						Rverts[0]->Xslope = (Rverts[0]->end[0]-Rverts[0]->start[0])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->Zslope = (Rverts[0]->end[2]-Rverts[0]->start[2])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						memcpy(Rverts[0]->vert1Norm,normList[0],3*sizeof(float));
						memcpy(Rverts[0]->currentNorm,normList[0],3*sizeof(float));
						memcpy(Rverts[0]->vert2Norm,normList[1],3*sizeof(float));
						Rverts[0]->NXslope = (Rverts[0]->vert2Norm[0]-Rverts[0]->vert1Norm[0])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->NYslope = (Rverts[0]->vert2Norm[1]-Rverts[0]->vert1Norm[1])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->NZslope = (Rverts[0]->vert2Norm[2]-Rverts[0]->vert1Norm[2])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->Ustart = uvList[0][0]/((Rverts[0]->start[2]/(MAXINT-Rverts[0]->start[2]))+1);
						Rverts[0]->Uend = uvList[1][0]/((Rverts[0]->end[2]/(MAXINT-Rverts[0]->end[2]))+1);
						Rverts[0]->Ucurrent = uvList[0][0]/((Rverts[0]->start[2]/(MAXINT-Rverts[0]->start[2]))+1);
						Rverts[0]->Vstart = uvList[0][1]/((Rverts[0]->start[2]/(MAXINT-Rverts[0]->start[2]))+1);
						Rverts[0]->Vend = uvList[1][1]/((Rverts[0]->end[2]/(MAXINT-Rverts[0]->end[2]))+1);
						Rverts[0]->Vcurrent = uvList[0][1]/((Rverts[0]->start[2]/(MAXINT-Rverts[0]->start[2]))+1);
						Rverts[0]->Uslope = (Rverts[0]->Uend-Rverts[0]->Ustart)/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->Vslope = (Rverts[0]->Vend-Rverts[0]->Vstart)/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[1] = new DDA();
						memcpy(Rverts[1]->start,vertsList[1],3*sizeof(float));
						memcpy(Rverts[1]->end,vertsList[2],3*sizeof(float));
						memcpy(Rverts[1]->current,vertsList[1],3*sizeof(float));
						Rverts[1]->Xslope = (Rverts[1]->end[0]-Rverts[1]->start[0])/(Rverts[1]->end[1]-Rverts[1]->start[1]);
						Rverts[1]->Zslope = (Rverts[1]->end[2]-Rverts[1]->start[2])/(Rverts[1]->end[1]-Rverts[1]->start[1]);
						memcpy(Rverts[1]->vert1Norm,normList[1],3*sizeof(float));
						memcpy(Rverts[1]->currentNorm,normList[1],3*sizeof(float));
						memcpy(Rverts[1]->vert2Norm,normList[2],3*sizeof(float));
						Rverts[1]->NXslope = (Rverts[1]->vert2Norm[0]-Rverts[1]->vert1Norm[0])/(Rverts[1]->end[1]-Rverts[1]->start[1]);
						Rverts[1]->NYslope = (Rverts[1]->vert2Norm[1]-Rverts[1]->vert1Norm[1])/(Rverts[1]->end[1]-Rverts[1]->start[1]);
						Rverts[1]->NZslope = (Rverts[1]->vert2Norm[2]-Rverts[1]->vert1Norm[2])/(Rverts[1]->end[1]-Rverts[1]->start[1]);
						Rverts[1]->Ustart = uvList[1][0]/((Rverts[1]->start[2]/(MAXINT-Rverts[1]->start[2]))+1);
						Rverts[1]->Uend = uvList[2][0]/((Rverts[1]->end[2]/(MAXINT-Rverts[1]->end[2]))+1);
						Rverts[1]->Ucurrent = uvList[1][0]/((Rverts[1]->start[2]/(MAXINT-Rverts[1]->start[2]))+1);
						Rverts[1]->Vstart = uvList[1][1]/((Rverts[1]->start[2]/(MAXINT-Rverts[1]->start[2]))+1);
						Rverts[1]->Vend = uvList[2][1]/((Rverts[1]->end[2]/(MAXINT-Rverts[1]->end[2]))+1);
						Rverts[1]->Vcurrent = uvList[1][1]/((Rverts[1]->start[2]/(MAXINT-Rverts[1]->start[2]))+1);
						Rverts[1]->Uslope = (Rverts[1]->Uend-Rverts[1]->Ustart)/(Rverts[1]->end[1]-Rverts[1]->start[1]);
						Rverts[1]->Vslope = (Rverts[1]->Vend-Rverts[1]->Vstart)/(Rverts[1]->end[1]-Rverts[1]->start[1]);

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
						memcpy(Lverts[0]->vert1Norm,normList[0],3*sizeof(float));
						memcpy(Lverts[0]->currentNorm,normList[0],3*sizeof(float));
						memcpy(Lverts[0]->vert2Norm,normList[1],3*sizeof(float));
						Lverts[0]->NXslope = (Lverts[0]->vert2Norm[0]-Lverts[0]->vert1Norm[0])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[0]->NYslope = (Lverts[0]->vert2Norm[1]-Lverts[0]->vert1Norm[1])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[0]->NZslope = (Lverts[0]->vert2Norm[2]-Lverts[0]->vert1Norm[2])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[0]->Ustart = uvList[0][0]/((Lverts[0]->start[2]/(MAXINT-Lverts[0]->start[2]))+1);
						Lverts[0]->Uend = uvList[1][0]/((Lverts[0]->end[2]/(MAXINT-Lverts[0]->end[2]))+1);
						Lverts[0]->Ucurrent = uvList[0][0]/((Lverts[0]->start[2]/(MAXINT-Lverts[0]->start[2]))+1);
						Lverts[0]->Vstart = uvList[0][1]/((Lverts[0]->start[2]/(MAXINT-Lverts[0]->start[2]))+1);
						Lverts[0]->Vend = uvList[1][1]/((Lverts[0]->end[2]/(MAXINT-Lverts[0]->end[2]))+1);
						Lverts[0]->Vcurrent = uvList[0][1]/((Lverts[0]->start[2]/(MAXINT-Lverts[0]->start[2]))+1);
						Lverts[0]->Uslope = (Lverts[0]->Uend-Lverts[0]->Ustart)/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[0]->Vslope = (Lverts[0]->Vend-Lverts[0]->Vstart)/(Lverts[0]->end[1]-Lverts[0]->start[1]);
						Lverts[1] = new DDA();
						memcpy(Lverts[1]->start,vertsList[1],3*sizeof(float));
						memcpy(Lverts[1]->end,vertsList[2],3*sizeof(float));
						memcpy(Lverts[1]->current,vertsList[1],3*sizeof(float));
						//float tempfl = Lverts[1]->end;//-*Lverts[1]->start[1]);///(*Lverts[1]->end[1]-*Lverts[1]->start[1]);
						Lverts[1]->Xslope = (Lverts[1]->end[0]-Lverts[1]->start[0])/(Lverts[1]->end[1]-Lverts[1]->start[1]);
						Lverts[1]->Zslope = (Lverts[1]->end[2]-Lverts[1]->start[2])/(Lverts[1]->end[1]-Lverts[1]->start[1]);
						memcpy(Lverts[1]->vert1Norm,normList[1],3*sizeof(float));
						memcpy(Lverts[1]->currentNorm,normList[1],3*sizeof(float));
						memcpy(Lverts[1]->vert2Norm,normList[2],3*sizeof(float));
						Lverts[1]->NXslope = (Lverts[1]->vert2Norm[0]-Lverts[1]->vert1Norm[0])/(Lverts[1]->end[1]-Lverts[1]->start[1]);
						Lverts[1]->NYslope = (Lverts[1]->vert2Norm[1]-Lverts[1]->vert1Norm[1])/(Lverts[1]->end[1]-Lverts[1]->start[1]);
						Lverts[1]->NZslope = (Lverts[1]->vert2Norm[2]-Lverts[1]->vert1Norm[2])/(Lverts[1]->end[1]-Lverts[1]->start[1]);
						Lverts[1]->Ustart = uvList[1][0]/((Lverts[1]->start[2]/(MAXINT-Lverts[1]->start[2]))+1);
						Lverts[1]->Uend = uvList[2][0]/((Lverts[1]->end[2]/(MAXINT-Lverts[1]->end[2]))+1);
						Lverts[1]->Ucurrent = uvList[1][0]/((Lverts[1]->start[2]/(MAXINT-Lverts[1]->start[2]))+1);
						Lverts[1]->Vstart = uvList[1][1]/((Lverts[1]->start[2]/(MAXINT-Lverts[1]->start[2]))+1);
						Lverts[1]->Vend = uvList[2][1]/((Lverts[1]->end[2]/(MAXINT-Lverts[1]->end[2]))+1);
						Lverts[1]->Vcurrent = uvList[1][1]/((Lverts[1]->start[2]/(MAXINT-Lverts[1]->start[2]))+1);
						Lverts[1]->Uslope = (Lverts[1]->Uend-Lverts[1]->Ustart)/(Lverts[1]->end[1]-Lverts[1]->start[1]);
						Lverts[1]->Vslope = (Lverts[1]->Vend-Lverts[1]->Vstart)/(Lverts[1]->end[1]-Lverts[1]->start[1]);

						//R edge DDA
						Rverts[0] = new DDA();
						memcpy(Rverts[0]->start,vertsList[0],3*sizeof(float));
						memcpy(Rverts[0]->end,vertsList[2],3*sizeof(float));
						memcpy(Rverts[0]->current,vertsList[0],3*sizeof(float));
						Rverts[0]->Xslope = (Rverts[0]->end[0]-Rverts[0]->start[0])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->Zslope = (Rverts[0]->end[2]-Rverts[0]->start[2])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						memcpy(Rverts[0]->vert1Norm,normList[0],3*sizeof(float));
						memcpy(Rverts[0]->currentNorm,normList[0],3*sizeof(float));
						memcpy(Rverts[0]->vert2Norm,normList[2],3*sizeof(float));
						Rverts[0]->NXslope = (Rverts[0]->vert2Norm[0]-Rverts[0]->vert1Norm[0])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->NYslope = (Rverts[0]->vert2Norm[1]-Rverts[0]->vert1Norm[1])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->NZslope = (Rverts[0]->vert2Norm[2]-Rverts[0]->vert1Norm[2])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->Ustart = uvList[0][0]/((Rverts[0]->start[2]/(MAXINT-Rverts[0]->start[2]))+1);
						Rverts[0]->Uend = uvList[2][0]/((Rverts[0]->end[2]/(MAXINT-Rverts[0]->end[2]))+1);
						Rverts[0]->Ucurrent = uvList[0][0]/((Rverts[0]->start[2]/(MAXINT-Rverts[0]->start[2]))+1);
						Rverts[0]->Vstart = uvList[0][1]/((Rverts[0]->start[2]/(MAXINT-Rverts[0]->start[2]))+1);
						Rverts[0]->Vend = uvList[2][1]/((Rverts[0]->end[2]/(MAXINT-Rverts[0]->end[2]))+1);
						Rverts[0]->Vcurrent = uvList[0][1]/((Rverts[0]->start[2]/(MAXINT-Rverts[0]->start[2]))+1);
						Rverts[0]->Uslope = (Rverts[0]->Uend-Rverts[0]->Ustart)/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						Rverts[0]->Vslope = (Rverts[0]->Vend-Rverts[0]->Vstart)/(Rverts[0]->end[1]-Rverts[0]->start[1]);
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
					//if our interpolation type is Gz_Color comput the color at each vert for interpolation
					if(render->interp_mode == GZ_COLOR){
						if(Lverts[1] == NULL){
							calculateColorGouraud(render,Lverts[0]->vert1Norm,Lverts[0]->vert1Color);
							calculateColorGouraud(render,Lverts[0]->vert2Norm,Lverts[0]->vert2Color);
							memcpy(Lverts[0]->currentColor,Lverts[0]->vert1Color,3*sizeof(float));
							Lverts[0]->Rslope = (Lverts[0]->vert2Color[0]-Lverts[0]->vert1Color[0])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
							Lverts[0]->Gslope = (Lverts[0]->vert2Color[1]-Lverts[0]->vert1Color[1])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
							Lverts[0]->Bslope = (Lverts[0]->vert2Color[2]-Lverts[0]->vert1Color[2])/(Lverts[0]->end[1]-Lverts[0]->start[1]);

							calculateColorGouraud(render,Rverts[0]->vert1Norm,Rverts[0]->vert1Color);
							calculateColorGouraud(render,Rverts[0]->vert2Norm,Rverts[0]->vert2Color);
							memcpy(Rverts[0]->currentColor,Rverts[0]->vert1Color,3*sizeof(float));
							Rverts[0]->Rslope = (Rverts[0]->vert2Color[0]-Rverts[0]->vert1Color[0])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
							Rverts[0]->Gslope = (Rverts[0]->vert2Color[1]-Rverts[0]->vert1Color[1])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
							Rverts[0]->Bslope = (Rverts[0]->vert2Color[2]-Rverts[0]->vert1Color[2])/(Rverts[0]->end[1]-Rverts[0]->start[1]);

							calculateColorGouraud(render,Rverts[1]->vert1Norm,Rverts[1]->vert1Color);
							calculateColorGouraud(render,Rverts[1]->vert2Norm,Rverts[1]->vert2Color);
							memcpy(Rverts[1]->currentColor,Rverts[1]->vert1Color,3*sizeof(float));
							Rverts[1]->Rslope = (Rverts[1]->vert2Color[0]-Rverts[1]->vert1Color[0])/(Rverts[1]->end[1]-Rverts[1]->start[1]);
							Rverts[1]->Gslope = (Rverts[1]->vert2Color[1]-Rverts[1]->vert1Color[1])/(Rverts[1]->end[1]-Rverts[1]->start[1]);
							Rverts[1]->Bslope = (Rverts[1]->vert2Color[2]-Rverts[1]->vert1Color[2])/(Rverts[1]->end[1]-Rverts[1]->start[1]);

						}
						else{
							calculateColorGouraud(render,Lverts[0]->vert1Norm,Lverts[0]->vert1Color);
							calculateColorGouraud(render,Lverts[0]->vert2Norm,Lverts[0]->vert2Color);
							memcpy(Lverts[0]->currentColor,Lverts[0]->vert1Color,3*sizeof(float));
							Lverts[0]->Rslope = (Lverts[0]->vert2Color[0]-Lverts[0]->vert1Color[0])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
							Lverts[0]->Gslope = (Lverts[0]->vert2Color[1]-Lverts[0]->vert1Color[1])/(Lverts[0]->end[1]-Lverts[0]->start[1]);
							Lverts[0]->Bslope = (Lverts[0]->vert2Color[2]-Lverts[0]->vert1Color[2])/(Lverts[0]->end[1]-Lverts[0]->start[1]);

							calculateColorGouraud(render,Lverts[1]->vert1Norm,Lverts[1]->vert1Color);
							calculateColorGouraud(render,Lverts[1]->vert2Norm,Lverts[1]->vert2Color);
							memcpy(Lverts[1]->currentColor,Lverts[1]->vert1Color,3*sizeof(float));
							Lverts[1]->Rslope = (Lverts[1]->vert2Color[0]-Lverts[1]->vert1Color[0])/(Lverts[1]->end[1]-Lverts[1]->start[1]);
							Lverts[1]->Gslope = (Lverts[1]->vert2Color[1]-Lverts[1]->vert1Color[1])/(Lverts[1]->end[1]-Lverts[1]->start[1]);
							Lverts[1]->Bslope = (Lverts[1]->vert2Color[2]-Lverts[1]->vert1Color[2])/(Lverts[1]->end[1]-Lverts[1]->start[1]);

							calculateColorGouraud(render,Rverts[0]->vert1Norm,Rverts[0]->vert1Color);
							calculateColorGouraud(render,Rverts[0]->vert2Norm,Rverts[0]->vert2Color);
							memcpy(Rverts[0]->currentColor,Rverts[0]->vert1Color,3*sizeof(float));
							Rverts[0]->Rslope = (Rverts[0]->vert2Color[0]-Rverts[0]->vert1Color[0])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
							Rverts[0]->Gslope = (Rverts[0]->vert2Color[1]-Rverts[0]->vert1Color[1])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
							Rverts[0]->Bslope = (Rverts[0]->vert2Color[2]-Rverts[0]->vert1Color[2])/(Rverts[0]->end[1]-Rverts[0]->start[1]);
						}
					}

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

								if(render->interp_mode == GZ_COLOR){
									Lverts[0]->currentColor[0] = Lverts[0]->currentColor[0]+(Lverts[0]->Rslope*dY);
									Lverts[0]->currentColor[1] = Lverts[0]->currentColor[1]+(Lverts[0]->Gslope*dY);
									Lverts[0]->currentColor[2] = Lverts[0]->currentColor[2]+(Lverts[0]->Bslope*dY);
									Rverts[0]->currentColor[0] = Rverts[0]->currentColor[0]+(Rverts[0]->Rslope*dY);
									Rverts[0]->currentColor[1] = Rverts[0]->currentColor[1]+(Rverts[0]->Gslope*dY);
									Rverts[0]->currentColor[2] = Rverts[0]->currentColor[2]+(Rverts[0]->Bslope*dY);
								}
								else if(render->interp_mode == GZ_NORMALS){
									Lverts[0]->currentNorm[0] = Lverts[0]->currentNorm[0]+(Lverts[0]->NXslope*dY);
									Lverts[0]->currentNorm[1] = Lverts[0]->currentNorm[1]+(Lverts[0]->NYslope*dY);
									Lverts[0]->currentNorm[2] = Lverts[0]->currentNorm[2]+(Lverts[0]->NZslope*dY);
									Rverts[0]->currentNorm[0] = Rverts[0]->currentNorm[0]+(Rverts[0]->NXslope*dY);
									Rverts[0]->currentNorm[1] = Rverts[0]->currentNorm[1]+(Rverts[0]->NYslope*dY);
									Rverts[0]->currentNorm[2] = Rverts[0]->currentNorm[2]+(Rverts[0]->NZslope*dY);
								}

								Lverts[0]->Ucurrent = Lverts[0]->Ucurrent + (Lverts[0]->Uslope*dY);
								Lverts[0]->Vcurrent = Lverts[0]->Vcurrent + (Lverts[0]->Vslope*dY);
								Rverts[0]->Ucurrent = Rverts[0]->Ucurrent + (Rverts[0]->Uslope*dY);
								Rverts[0]->Vcurrent = Rverts[0]->Vcurrent + (Rverts[0]->Vslope*dY);

								//set up span DDA
								minX = Lverts[0]->current[0];
								MaxX = Rverts[0]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[0]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[0]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[0]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);

								if(render->interp_mode == GZ_COLOR){
									memcpy(spanX.vert1Color,Lverts[0]->currentColor,3*sizeof(float));
									memcpy(spanX.vert2Color,Rverts[0]->currentColor,3*sizeof(float));
									memcpy(spanX.currentColor,Lverts[0]->currentColor,3*sizeof(float));
									spanX.Rslope = (spanX.vert2Color[0]-spanX.vert1Color[0])/(spanX.end[0]-spanX.start[0]);
									spanX.Gslope = (spanX.vert2Color[1]-spanX.vert1Color[1])/(spanX.end[0]-spanX.start[0]);
									spanX.Bslope = (spanX.vert2Color[2]-spanX.vert1Color[2])/(spanX.end[0]-spanX.start[0]);
								}
								else if(render->interp_mode == GZ_NORMALS){
									memcpy(spanX.vert1Norm,Lverts[0]->currentNorm,3*sizeof(float));
									memcpy(spanX.vert2Norm,Rverts[0]->currentNorm,3*sizeof(float));
									memcpy(spanX.currentNorm,Lverts[0]->currentNorm,3*sizeof(float));
									spanX.NXslope = (spanX.vert2Norm[0]-spanX.vert1Norm[0])/(spanX.end[0]-spanX.start[0]);
									spanX.NYslope = (spanX.vert2Norm[1]-spanX.vert1Norm[1])/(spanX.end[0]-spanX.start[0]);
									spanX.NZslope = (spanX.vert2Norm[2]-spanX.vert1Norm[2])/(spanX.end[0]-spanX.start[0]);
								}

								spanX.Ustart = Lverts[0]->Ucurrent;
								spanX.Vstart = Lverts[0]->Vcurrent;
								spanX.Uend = Rverts[0]->Ucurrent;
								spanX.Vend = Rverts[0]->Vcurrent;
								spanX.Ucurrent = Lverts[0]->Ucurrent;
								spanX.Vcurrent = Lverts[0]->Vcurrent;
								spanX.Uslope = (spanX.Uend-spanX.Ustart)/(spanX.end[0]-spanX.start[0]);
								spanX.Vslope = (spanX.Vend-spanX.Vstart)/(spanX.end[0]-spanX.start[0]);

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
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0]+1 != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								

							}
							else if(yLine <= midY){
								dY = 1;
								Lverts[0]->current[0] = Lverts[0]->current[0]+(Lverts[0]->Xslope*dY);
								Lverts[0]->current[1] = Lverts[0]->current[1]+dY;
								Lverts[0]->current[2] = Lverts[0]->current[2]+(Lverts[0]->Zslope*dY);
								Rverts[0]->current[0] = Rverts[0]->current[0]+(Rverts[0]->Xslope*dY);
								Rverts[0]->current[1] = Rverts[0]->current[1]+dY;
								Rverts[0]->current[2] = Rverts[0]->current[2]+(Rverts[0]->Zslope*dY);

								if(render->interp_mode == GZ_COLOR){
									Lverts[0]->currentColor[0] = Lverts[0]->currentColor[0]+(Lverts[0]->Rslope*dY);
									Lverts[0]->currentColor[1] = Lverts[0]->currentColor[1]+(Lverts[0]->Gslope*dY);
									Lverts[0]->currentColor[2] = Lverts[0]->currentColor[2]+(Lverts[0]->Bslope*dY);
									Rverts[0]->currentColor[0] = Rverts[0]->currentColor[0]+(Rverts[0]->Rslope*dY);
									Rverts[0]->currentColor[1] = Rverts[0]->currentColor[1]+(Rverts[0]->Gslope*dY);
									Rverts[0]->currentColor[2] = Rverts[0]->currentColor[2]+(Rverts[0]->Bslope*dY);
								}
								else if(render->interp_mode == GZ_NORMALS){
									Lverts[0]->currentNorm[0] = Lverts[0]->currentNorm[0]+(Lverts[0]->NXslope*dY);
									Lverts[0]->currentNorm[1] = Lverts[0]->currentNorm[1]+(Lverts[0]->NYslope*dY);
									Lverts[0]->currentNorm[2] = Lverts[0]->currentNorm[2]+(Lverts[0]->NZslope*dY);
									Rverts[0]->currentNorm[0] = Rverts[0]->currentNorm[0]+(Rverts[0]->NXslope*dY);
									Rverts[0]->currentNorm[1] = Rverts[0]->currentNorm[1]+(Rverts[0]->NYslope*dY);
									Rverts[0]->currentNorm[2] = Rverts[0]->currentNorm[2]+(Rverts[0]->NZslope*dY);
								}

								Lverts[0]->Ucurrent = Lverts[0]->Ucurrent + (Lverts[0]->Uslope*dY);
								Lverts[0]->Vcurrent = Lverts[0]->Vcurrent + (Lverts[0]->Vslope*dY);
								Rverts[0]->Ucurrent = Rverts[0]->Ucurrent + (Rverts[0]->Uslope*dY);
								Rverts[0]->Vcurrent = Rverts[0]->Vcurrent + (Rverts[0]->Vslope*dY);

								//set up span DDA
								minX = Lverts[0]->current[0];
								MaxX = Rverts[0]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[0]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[0]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[0]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);

								if(render->interp_mode == GZ_COLOR){
									memcpy(spanX.vert1Color,Lverts[0]->currentColor,3*sizeof(float));
									memcpy(spanX.vert2Color,Rverts[0]->currentColor,3*sizeof(float));
									memcpy(spanX.currentColor,Lverts[0]->currentColor,3*sizeof(float));
									spanX.Rslope = (spanX.vert2Color[0]-spanX.vert1Color[0])/(spanX.end[0]-spanX.start[0]);
									spanX.Gslope = (spanX.vert2Color[1]-spanX.vert1Color[1])/(spanX.end[0]-spanX.start[0]);
									spanX.Bslope = (spanX.vert2Color[2]-spanX.vert1Color[2])/(spanX.end[0]-spanX.start[0]);
								}
								else if(render->interp_mode == GZ_NORMALS){
									memcpy(spanX.vert1Norm,Lverts[0]->currentNorm,3*sizeof(float));
									memcpy(spanX.vert2Norm,Rverts[0]->currentNorm,3*sizeof(float));
									memcpy(spanX.currentNorm,Lverts[0]->currentNorm,3*sizeof(float));
									spanX.NXslope = (spanX.vert2Norm[0]-spanX.vert1Norm[0])/(spanX.end[0]-spanX.start[0]);
									spanX.NYslope = (spanX.vert2Norm[1]-spanX.vert1Norm[1])/(spanX.end[0]-spanX.start[0]);
									spanX.NZslope = (spanX.vert2Norm[2]-spanX.vert1Norm[2])/(spanX.end[0]-spanX.start[0]);
								}

								spanX.Ustart = Lverts[0]->Ucurrent;
								spanX.Vstart = Lverts[0]->Vcurrent;
								spanX.Uend = Rverts[0]->Ucurrent;
								spanX.Vend = Rverts[0]->Vcurrent;
								spanX.Ucurrent = Lverts[0]->Ucurrent;
								spanX.Vcurrent = Lverts[0]->Vcurrent;
								spanX.Uslope = (spanX.Uend-spanX.Ustart)/(spanX.end[0]-spanX.start[0]);
								spanX.Vslope = (spanX.Vend-spanX.Vstart)/(spanX.end[0]-spanX.start[0]);

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
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0]+1 != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								
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

								if(render->interp_mode == GZ_COLOR){
									Lverts[0]->currentColor[0] = Lverts[0]->currentColor[0]+(Lverts[0]->Rslope*dYL);
									Lverts[0]->currentColor[1] = Lverts[0]->currentColor[1]+(Lverts[0]->Gslope*dYL);
									Lverts[0]->currentColor[2] = Lverts[0]->currentColor[2]+(Lverts[0]->Bslope*dYL);
									Rverts[1]->currentColor[0] = Rverts[1]->currentColor[0]+(Rverts[1]->Rslope*dYR);
									Rverts[1]->currentColor[1] = Rverts[1]->currentColor[1]+(Rverts[1]->Gslope*dYR);
									Rverts[1]->currentColor[2] = Rverts[1]->currentColor[2]+(Rverts[1]->Bslope*dYR);
								}
								else if(render->interp_mode == GZ_NORMALS){
									Lverts[0]->currentNorm[0] = Lverts[0]->currentNorm[0]+(Lverts[0]->NXslope*dYL);
									Lverts[0]->currentNorm[1] = Lverts[0]->currentNorm[1]+(Lverts[0]->NYslope*dYL);
									Lverts[0]->currentNorm[2] = Lverts[0]->currentNorm[2]+(Lverts[0]->NZslope*dYL);
									Rverts[1]->currentNorm[0] = Rverts[1]->currentNorm[0]+(Rverts[1]->NXslope*dYR);
									Rverts[1]->currentNorm[1] = Rverts[1]->currentNorm[1]+(Rverts[1]->NYslope*dYR);
									Rverts[1]->currentNorm[2] = Rverts[1]->currentNorm[2]+(Rverts[1]->NZslope*dYR);
								}

								Lverts[0]->Ucurrent = Lverts[0]->Ucurrent + (Lverts[0]->Uslope*dYL);
								Lverts[0]->Vcurrent = Lverts[0]->Vcurrent + (Lverts[0]->Vslope*dYL);
								Rverts[1]->Ucurrent = Rverts[1]->Ucurrent + (Rverts[1]->Uslope*dYR);
								Rverts[1]->Vcurrent = Rverts[1]->Vcurrent + (Rverts[1]->Vslope*dYR);

								//set up span DDA
								minX = Lverts[0]->current[0];
								MaxX = Rverts[1]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[0]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[1]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[0]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);

								if(render->interp_mode == GZ_COLOR){
									memcpy(spanX.vert1Color,Lverts[0]->currentColor,3*sizeof(float));
									memcpy(spanX.vert2Color,Rverts[1]->currentColor,3*sizeof(float));
									memcpy(spanX.currentColor,Lverts[0]->currentColor,3*sizeof(float));
									spanX.Rslope = (spanX.vert2Color[0]-spanX.vert1Color[0])/(spanX.end[0]-spanX.start[0]);
									spanX.Gslope = (spanX.vert2Color[1]-spanX.vert1Color[1])/(spanX.end[0]-spanX.start[0]);
									spanX.Bslope = (spanX.vert2Color[2]-spanX.vert1Color[2])/(spanX.end[0]-spanX.start[0]);
								}
								else if(render->interp_mode == GZ_NORMALS){
									memcpy(spanX.vert1Norm,Lverts[0]->currentNorm,3*sizeof(float));
									memcpy(spanX.vert2Norm,Rverts[1]->currentNorm,3*sizeof(float));
									memcpy(spanX.currentNorm,Lverts[0]->currentNorm,3*sizeof(float));
									spanX.NXslope = (spanX.vert2Norm[0]-spanX.vert1Norm[0])/(spanX.end[0]-spanX.start[0]);
									spanX.NYslope = (spanX.vert2Norm[1]-spanX.vert1Norm[1])/(spanX.end[0]-spanX.start[0]);
									spanX.NZslope = (spanX.vert2Norm[2]-spanX.vert1Norm[2])/(spanX.end[0]-spanX.start[0]);
								}

								spanX.Ustart = Lverts[0]->Ucurrent;
								spanX.Vstart = Lverts[0]->Vcurrent;
								spanX.Uend = Rverts[1]->Ucurrent;
								spanX.Vend = Rverts[1]->Vcurrent;
								spanX.Ucurrent = Lverts[0]->Ucurrent;
								spanX.Vcurrent = Lverts[0]->Vcurrent;
								spanX.Uslope = (spanX.Uend-spanX.Ustart)/(spanX.end[0]-spanX.start[0]);
								spanX.Vslope = (spanX.Vend-spanX.Vstart)/(spanX.end[0]-spanX.start[0]);

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
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0] != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								
							}
							else if((yLine == midY+1 && midY != MaxY)){
								dY = midY - Rverts[1]->start[1]+1;
								Lverts[0]->current[0] = Lverts[0]->current[0]+(Lverts[0]->Xslope*1);
								Lverts[0]->current[1] = Lverts[0]->current[1]+1;
								Lverts[0]->current[2] = Lverts[0]->current[2]+(Lverts[0]->Zslope*1);
								Rverts[1]->current[0] = Rverts[1]->current[0]+(Rverts[1]->Xslope*dY);
								Rverts[1]->current[1] = Rverts[1]->current[1]+dY;
								Rverts[1]->current[2] = Rverts[1]->current[2]+(Rverts[1]->Zslope*dY);

								if(render->interp_mode == GZ_COLOR){
									Lverts[0]->currentColor[0] = Lverts[0]->currentColor[0]+(Lverts[0]->Rslope*1);
									Lverts[0]->currentColor[1] = Lverts[0]->currentColor[1]+(Lverts[0]->Gslope*1);
									Lverts[0]->currentColor[2] = Lverts[0]->currentColor[2]+(Lverts[0]->Bslope*1);
									Rverts[1]->currentColor[0] = Rverts[1]->currentColor[0]+(Rverts[1]->Rslope*dY);
									Rverts[1]->currentColor[1] = Rverts[1]->currentColor[1]+(Rverts[1]->Gslope*dY);
									Rverts[1]->currentColor[2] = Rverts[1]->currentColor[2]+(Rverts[1]->Bslope*dY);
								}
								else if(render->interp_mode == GZ_NORMALS){
									Lverts[0]->currentNorm[0] = Lverts[0]->currentNorm[0]+(Lverts[0]->NXslope*1);
									Lverts[0]->currentNorm[1] = Lverts[0]->currentNorm[1]+(Lverts[0]->NYslope*1);
									Lverts[0]->currentNorm[2] = Lverts[0]->currentNorm[2]+(Lverts[0]->NZslope*1);
									Rverts[1]->currentNorm[0] = Rverts[1]->currentNorm[0]+(Rverts[1]->NXslope*dY);
									Rverts[1]->currentNorm[1] = Rverts[1]->currentNorm[1]+(Rverts[1]->NYslope*dY);
									Rverts[1]->currentNorm[2] = Rverts[1]->currentNorm[2]+(Rverts[1]->NZslope*dY);
								}

								Lverts[0]->Ucurrent = Lverts[0]->Ucurrent + (Lverts[0]->Uslope*1);
								Lverts[0]->Vcurrent = Lverts[0]->Vcurrent + (Lverts[0]->Vslope*1);
								Rverts[1]->Ucurrent = Rverts[1]->Ucurrent + (Rverts[1]->Uslope*dY);
								Rverts[1]->Vcurrent = Rverts[1]->Vcurrent + (Rverts[1]->Vslope*dY);

								//set up span DDA
								minX = Lverts[0]->current[0];
								MaxX = Rverts[1]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[0]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[1]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[0]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);

								if(render->interp_mode == GZ_COLOR){
									memcpy(spanX.vert1Color,Lverts[0]->currentColor,3*sizeof(float));
									memcpy(spanX.vert2Color,Rverts[1]->currentColor,3*sizeof(float));
									memcpy(spanX.currentColor,Lverts[0]->currentColor,3*sizeof(float));
									spanX.Rslope = (spanX.vert2Color[0]-spanX.vert1Color[0])/(spanX.end[0]-spanX.start[0]);
									spanX.Gslope = (spanX.vert2Color[1]-spanX.vert1Color[1])/(spanX.end[0]-spanX.start[0]);
									spanX.Bslope = (spanX.vert2Color[2]-spanX.vert1Color[2])/(spanX.end[0]-spanX.start[0]);
								}
								else if(render->interp_mode == GZ_NORMALS){
									memcpy(spanX.vert1Norm,Lverts[0]->currentNorm,3*sizeof(float));
									memcpy(spanX.vert2Norm,Rverts[1]->currentNorm,3*sizeof(float));
									memcpy(spanX.currentNorm,Lverts[0]->currentNorm,3*sizeof(float));
									spanX.NXslope = (spanX.vert2Norm[0]-spanX.vert1Norm[0])/(spanX.end[0]-spanX.start[0]);
									spanX.NYslope = (spanX.vert2Norm[1]-spanX.vert1Norm[1])/(spanX.end[0]-spanX.start[0]);
									spanX.NZslope = (spanX.vert2Norm[2]-spanX.vert1Norm[2])/(spanX.end[0]-spanX.start[0]);
								}

								spanX.Ustart = Lverts[0]->Ucurrent;
								spanX.Vstart = Lverts[0]->Vcurrent;
								spanX.Uend = Rverts[1]->Ucurrent;
								spanX.Vend = Rverts[1]->Vcurrent;
								spanX.Ucurrent = Lverts[0]->Ucurrent;
								spanX.Vcurrent = Lverts[0]->Vcurrent;
								spanX.Uslope = (spanX.Uend-spanX.Ustart)/(spanX.end[0]-spanX.start[0]);
								spanX.Vslope = (spanX.Vend-spanX.Vstart)/(spanX.end[0]-spanX.start[0]);

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
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0]+1 != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								
							}
							else if(midY != MaxY && yLine > midY){
								dY = 1;
								Lverts[0]->current[0] = Lverts[0]->current[0]+(Lverts[0]->Xslope*dY);
								Lverts[0]->current[1] = Lverts[0]->current[1]+dY;
								Lverts[0]->current[2] = Lverts[0]->current[2]+(Lverts[0]->Zslope*dY);
								Rverts[1]->current[0] = Rverts[1]->current[0]+(Rverts[1]->Xslope*dY);
								Rverts[1]->current[1] = Rverts[1]->current[1]+dY;
								Rverts[1]->current[2] = Rverts[1]->current[2]+(Rverts[1]->Zslope*dY);

								if(render->interp_mode == GZ_COLOR){
									Lverts[0]->currentColor[0] = Lverts[0]->currentColor[0]+(Lverts[0]->Rslope*dY);
									Lverts[0]->currentColor[1] = Lverts[0]->currentColor[1]+(Lverts[0]->Gslope*dY);
									Lverts[0]->currentColor[2] = Lverts[0]->currentColor[2]+(Lverts[0]->Bslope*dY);
									Rverts[1]->currentColor[0] = Rverts[1]->currentColor[0]+(Rverts[1]->Rslope*dY);
									Rverts[1]->currentColor[1] = Rverts[1]->currentColor[1]+(Rverts[1]->Gslope*dY);
									Rverts[1]->currentColor[2] = Rverts[1]->currentColor[2]+(Rverts[1]->Bslope*dY);
								}
								else if(render->interp_mode == GZ_NORMALS){
									Lverts[0]->currentNorm[0] = Lverts[0]->currentNorm[0]+(Lverts[0]->NXslope*dY);
									Lverts[0]->currentNorm[1] = Lverts[0]->currentNorm[1]+(Lverts[0]->NYslope*dY);
									Lverts[0]->currentNorm[2] = Lverts[0]->currentNorm[2]+(Lverts[0]->NZslope*dY);
									Rverts[1]->currentNorm[0] = Rverts[1]->currentNorm[0]+(Rverts[1]->NXslope*dY);
									Rverts[1]->currentNorm[1] = Rverts[1]->currentNorm[1]+(Rverts[1]->NYslope*dY);
									Rverts[1]->currentNorm[2] = Rverts[1]->currentNorm[2]+(Rverts[1]->NZslope*dY);
								}

								Lverts[0]->Ucurrent = Lverts[0]->Ucurrent + (Lverts[0]->Uslope*dY);
								Lverts[0]->Vcurrent = Lverts[0]->Vcurrent + (Lverts[0]->Vslope*dY);
								Rverts[1]->Ucurrent = Rverts[1]->Ucurrent + (Rverts[1]->Uslope*dY);
								Rverts[1]->Vcurrent = Rverts[1]->Vcurrent + (Rverts[1]->Vslope*dY);

								//set up span DDA
								minX = Lverts[0]->current[0];
								MaxX = Rverts[1]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[0]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[1]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[0]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);

								if(render->interp_mode == GZ_COLOR){
									memcpy(spanX.vert1Color,Lverts[0]->currentColor,3*sizeof(float));
									memcpy(spanX.vert2Color,Rverts[1]->currentColor,3*sizeof(float));
									memcpy(spanX.currentColor,Lverts[0]->currentColor,3*sizeof(float));
									spanX.Rslope = (spanX.vert2Color[0]-spanX.vert1Color[0])/(spanX.end[0]-spanX.start[0]);
									spanX.Gslope = (spanX.vert2Color[1]-spanX.vert1Color[1])/(spanX.end[0]-spanX.start[0]);
									spanX.Bslope = (spanX.vert2Color[2]-spanX.vert1Color[2])/(spanX.end[0]-spanX.start[0]);
								}
								else if(render->interp_mode == GZ_NORMALS){
									memcpy(spanX.vert1Norm,Lverts[0]->currentNorm,3*sizeof(float));
									memcpy(spanX.vert2Norm,Rverts[1]->currentNorm,3*sizeof(float));
									memcpy(spanX.currentNorm,Lverts[0]->currentNorm,3*sizeof(float));
									spanX.NXslope = (spanX.vert2Norm[0]-spanX.vert1Norm[0])/(spanX.end[0]-spanX.start[0]);
									spanX.NYslope = (spanX.vert2Norm[1]-spanX.vert1Norm[1])/(spanX.end[0]-spanX.start[0]);
									spanX.NZslope = (spanX.vert2Norm[2]-spanX.vert1Norm[2])/(spanX.end[0]-spanX.start[0]);
								}

								spanX.Ustart = Lverts[0]->Ucurrent;
								spanX.Vstart = Lverts[0]->Vcurrent;
								spanX.Uend = Rverts[1]->Ucurrent;
								spanX.Vend = Rverts[1]->Vcurrent;
								spanX.Ucurrent = Lverts[0]->Ucurrent;
								spanX.Vcurrent = Lverts[0]->Vcurrent;
								spanX.Uslope = (spanX.Uend-spanX.Ustart)/(spanX.end[0]-spanX.start[0]);
								spanX.Vslope = (spanX.Vend-spanX.Vstart)/(spanX.end[0]-spanX.start[0]);

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
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0]+1 != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								
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

								if(render->interp_mode == GZ_COLOR){
									Lverts[0]->currentColor[0] = Lverts[0]->currentColor[0]+(Lverts[0]->Rslope*dY);
									Lverts[0]->currentColor[1] = Lverts[0]->currentColor[1]+(Lverts[0]->Gslope*dY);
									Lverts[0]->currentColor[2] = Lverts[0]->currentColor[2]+(Lverts[0]->Bslope*dY);
									Rverts[0]->currentColor[0] = Rverts[0]->currentColor[0]+(Rverts[0]->Rslope*dY);
									Rverts[0]->currentColor[1] = Rverts[0]->currentColor[1]+(Rverts[0]->Gslope*dY);
									Rverts[0]->currentColor[2] = Rverts[0]->currentColor[2]+(Rverts[0]->Bslope*dY);
								}
								else if(render->interp_mode == GZ_NORMALS){
									Lverts[0]->currentNorm[0] = Lverts[0]->currentNorm[0]+(Lverts[0]->NXslope*dY);
									Lverts[0]->currentNorm[1] = Lverts[0]->currentNorm[1]+(Lverts[0]->NYslope*dY);
									Lverts[0]->currentNorm[2] = Lverts[0]->currentNorm[2]+(Lverts[0]->NZslope*dY);
									Rverts[0]->currentNorm[0] = Rverts[0]->currentNorm[0]+(Rverts[0]->NXslope*dY);
									Rverts[0]->currentNorm[1] = Rverts[0]->currentNorm[1]+(Rverts[0]->NYslope*dY);
									Rverts[0]->currentNorm[2] = Rverts[0]->currentNorm[2]+(Rverts[0]->NZslope*dY);
								}

								Lverts[0]->Ucurrent = Lverts[0]->Ucurrent + (Lverts[0]->Uslope*dY);
								Lverts[0]->Vcurrent = Lverts[0]->Vcurrent + (Lverts[0]->Vslope*dY);
								Rverts[0]->Ucurrent = Rverts[0]->Ucurrent + (Rverts[0]->Uslope*dY);
								Rverts[0]->Vcurrent = Rverts[0]->Vcurrent + (Rverts[0]->Vslope*dY);

								//set up span DDA
								minX = Lverts[0]->current[0];
								MaxX = Rverts[0]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[0]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[0]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[0]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);

								if(render->interp_mode == GZ_COLOR){
									memcpy(spanX.vert1Color,Lverts[0]->currentColor,3*sizeof(float));
									memcpy(spanX.vert2Color,Rverts[0]->currentColor,3*sizeof(float));
									memcpy(spanX.currentColor,Lverts[0]->currentColor,3*sizeof(float));
									spanX.Rslope = (spanX.vert2Color[0]-spanX.vert1Color[0])/(spanX.end[0]-spanX.start[0]);
									spanX.Gslope = (spanX.vert2Color[1]-spanX.vert1Color[1])/(spanX.end[0]-spanX.start[0]);
									spanX.Bslope = (spanX.vert2Color[2]-spanX.vert1Color[2])/(spanX.end[0]-spanX.start[0]);
								}
								else if(render->interp_mode == GZ_NORMALS){
									memcpy(spanX.vert1Norm,Lverts[0]->currentNorm,3*sizeof(float));
									memcpy(spanX.vert2Norm,Rverts[0]->currentNorm,3*sizeof(float));
									memcpy(spanX.currentNorm,Lverts[0]->currentNorm,3*sizeof(float));
									spanX.NXslope = (spanX.vert2Norm[0]-spanX.vert1Norm[0])/(spanX.end[0]-spanX.start[0]);
									spanX.NYslope = (spanX.vert2Norm[1]-spanX.vert1Norm[1])/(spanX.end[0]-spanX.start[0]);
									spanX.NZslope = (spanX.vert2Norm[2]-spanX.vert1Norm[2])/(spanX.end[0]-spanX.start[0]);
								}

								spanX.Ustart = Lverts[0]->Ucurrent;
								spanX.Vstart = Lverts[0]->Vcurrent;
								spanX.Uend = Rverts[0]->Ucurrent;
								spanX.Vend = Rverts[0]->Vcurrent;
								spanX.Ucurrent = Lverts[0]->Ucurrent;
								spanX.Vcurrent = Lverts[0]->Vcurrent;
								spanX.Uslope = (spanX.Uend-spanX.Ustart)/(spanX.end[0]-spanX.start[0]);
								spanX.Vslope = (spanX.Vend-spanX.Vstart)/(spanX.end[0]-spanX.start[0]);

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
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0]+1 != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								

							}
							else if(yLine <= midY){
								dY = 1;
								Lverts[0]->current[0] = Lverts[0]->current[0]+(Lverts[0]->Xslope*dY);
								Lverts[0]->current[1] = Lverts[0]->current[1]+dY;
								Lverts[0]->current[2] = Lverts[0]->current[2]+(Lverts[0]->Zslope*dY);
								Rverts[0]->current[0] = Rverts[0]->current[0]+(Rverts[0]->Xslope*dY);
								Rverts[0]->current[1] = Rverts[0]->current[1]+dY;
								Rverts[0]->current[2] = Rverts[0]->current[2]+(Rverts[0]->Zslope*dY);

								if(render->interp_mode == GZ_COLOR){
									Lverts[0]->currentColor[0] = Lverts[0]->currentColor[0]+(Lverts[0]->Rslope*dY);
									Lverts[0]->currentColor[1] = Lverts[0]->currentColor[1]+(Lverts[0]->Gslope*dY);
									Lverts[0]->currentColor[2] = Lverts[0]->currentColor[2]+(Lverts[0]->Bslope*dY);
									Rverts[0]->currentColor[0] = Rverts[0]->currentColor[0]+(Rverts[0]->Rslope*dY);
									Rverts[0]->currentColor[1] = Rverts[0]->currentColor[1]+(Rverts[0]->Gslope*dY);
									Rverts[0]->currentColor[2] = Rverts[0]->currentColor[2]+(Rverts[0]->Bslope*dY);
								}
								else if(render->interp_mode == GZ_NORMALS){
									Lverts[0]->currentNorm[0] = Lverts[0]->currentNorm[0]+(Lverts[0]->NXslope*dY);
									Lverts[0]->currentNorm[1] = Lverts[0]->currentNorm[1]+(Lverts[0]->NYslope*dY);
									Lverts[0]->currentNorm[2] = Lverts[0]->currentNorm[2]+(Lverts[0]->NZslope*dY);
									Rverts[0]->currentNorm[0] = Rverts[0]->currentNorm[0]+(Rverts[0]->NXslope*dY);
									Rverts[0]->currentNorm[1] = Rverts[0]->currentNorm[1]+(Rverts[0]->NYslope*dY);
									Rverts[0]->currentNorm[2] = Rverts[0]->currentNorm[2]+(Rverts[0]->NZslope*dY);
								}

								Lverts[0]->Ucurrent = Lverts[0]->Ucurrent + (Lverts[0]->Uslope*dY);
								Lverts[0]->Vcurrent = Lverts[0]->Vcurrent + (Lverts[0]->Vslope*dY);
								Rverts[0]->Ucurrent = Rverts[0]->Ucurrent + (Rverts[0]->Uslope*dY);
								Rverts[0]->Vcurrent = Rverts[0]->Vcurrent + (Rverts[0]->Vslope*dY);

								//set up span DDA
								minX = Lverts[0]->current[0];
								MaxX = Rverts[0]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[0]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[0]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[0]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);

								if(render->interp_mode == GZ_COLOR){
									memcpy(spanX.vert1Color,Lverts[0]->currentColor,3*sizeof(float));
									memcpy(spanX.vert2Color,Rverts[0]->currentColor,3*sizeof(float));
									memcpy(spanX.currentColor,Lverts[0]->currentColor,3*sizeof(float));
									spanX.Rslope = (spanX.vert2Color[0]-spanX.vert1Color[0])/(spanX.end[0]-spanX.start[0]);
									spanX.Gslope = (spanX.vert2Color[1]-spanX.vert1Color[1])/(spanX.end[0]-spanX.start[0]);
									spanX.Bslope = (spanX.vert2Color[2]-spanX.vert1Color[2])/(spanX.end[0]-spanX.start[0]);
								}
								else if(render->interp_mode == GZ_NORMALS){
									memcpy(spanX.vert1Norm,Lverts[0]->currentNorm,3*sizeof(float));
									memcpy(spanX.vert2Norm,Rverts[0]->currentNorm,3*sizeof(float));
									memcpy(spanX.currentNorm,Lverts[0]->currentNorm,3*sizeof(float));
									spanX.NXslope = (spanX.vert2Norm[0]-spanX.vert1Norm[0])/(spanX.end[0]-spanX.start[0]);
									spanX.NYslope = (spanX.vert2Norm[1]-spanX.vert1Norm[1])/(spanX.end[0]-spanX.start[0]);
									spanX.NZslope = (spanX.vert2Norm[2]-spanX.vert1Norm[2])/(spanX.end[0]-spanX.start[0]);
								}

								spanX.Ustart = Lverts[0]->Ucurrent;
								spanX.Vstart = Lverts[0]->Vcurrent;
								spanX.Uend = Rverts[0]->Ucurrent;
								spanX.Vend = Rverts[0]->Vcurrent;
								spanX.Ucurrent = Lverts[0]->Ucurrent;
								spanX.Vcurrent = Lverts[0]->Vcurrent;
								spanX.Uslope = (spanX.Uend-spanX.Ustart)/(spanX.end[0]-spanX.start[0]);
								spanX.Vslope = (spanX.Vend-spanX.Vstart)/(spanX.end[0]-spanX.start[0]);

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
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0]+1 != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								
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

								if(render->interp_mode == GZ_COLOR){
									Lverts[1]->currentColor[0] = Lverts[1]->currentColor[0]+(Lverts[1]->Rslope*dYL);
									Lverts[1]->currentColor[1] = Lverts[1]->currentColor[1]+(Lverts[1]->Gslope*dYL);
									Lverts[1]->currentColor[2] = Lverts[1]->currentColor[2]+(Lverts[1]->Bslope*dYL);
									Rverts[0]->currentColor[0] = Rverts[0]->currentColor[0]+(Rverts[0]->Rslope*dYR);
									Rverts[0]->currentColor[1] = Rverts[0]->currentColor[1]+(Rverts[0]->Gslope*dYR);
									Rverts[0]->currentColor[2] = Rverts[0]->currentColor[2]+(Rverts[0]->Bslope*dYR);
								}
								else if(render->interp_mode == GZ_NORMALS){
									Lverts[1]->currentNorm[0] = Lverts[1]->currentNorm[0]+(Lverts[1]->NXslope*dYL);
									Lverts[1]->currentNorm[1] = Lverts[1]->currentNorm[1]+(Lverts[1]->NYslope*dYL);
									Lverts[1]->currentNorm[2] = Lverts[1]->currentNorm[2]+(Lverts[1]->NZslope*dYL);
									Rverts[0]->currentNorm[0] = Rverts[0]->currentNorm[0]+(Rverts[0]->NXslope*dYR);
									Rverts[0]->currentNorm[1] = Rverts[0]->currentNorm[1]+(Rverts[0]->NYslope*dYR);
									Rverts[0]->currentNorm[2] = Rverts[0]->currentNorm[2]+(Rverts[0]->NZslope*dYR);
								}

								Lverts[1]->Ucurrent = Lverts[1]->Ucurrent + (Lverts[1]->Uslope*dYL);
								Lverts[1]->Vcurrent = Lverts[1]->Vcurrent + (Lverts[1]->Vslope*dYL);
								Rverts[0]->Ucurrent = Rverts[0]->Ucurrent + (Rverts[0]->Uslope*dYR);
								Rverts[0]->Vcurrent = Rverts[0]->Vcurrent + (Rverts[0]->Vslope*dYR);

								//set up span DDA
								minX = Lverts[1]->current[0];
								MaxX = Rverts[0]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[1]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[0]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[1]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);

								if(render->interp_mode == GZ_COLOR){
									memcpy(spanX.vert1Color,Lverts[1]->currentColor,3*sizeof(float));
									memcpy(spanX.vert2Color,Rverts[0]->currentColor,3*sizeof(float));
									memcpy(spanX.currentColor,Lverts[1]->currentColor,3*sizeof(float));
									spanX.Rslope = (spanX.vert2Color[0]-spanX.vert1Color[0])/(spanX.end[0]-spanX.start[0]);
									spanX.Gslope = (spanX.vert2Color[1]-spanX.vert1Color[1])/(spanX.end[0]-spanX.start[0]);
									spanX.Bslope = (spanX.vert2Color[2]-spanX.vert1Color[2])/(spanX.end[0]-spanX.start[0]);
								}
								else if(render->interp_mode == GZ_NORMALS){
									memcpy(spanX.vert1Norm,Lverts[1]->currentNorm,3*sizeof(float));
									memcpy(spanX.vert2Norm,Rverts[0]->currentNorm,3*sizeof(float));
									memcpy(spanX.currentNorm,Lverts[1]->currentNorm,3*sizeof(float));
									spanX.NXslope = (spanX.vert2Norm[0]-spanX.vert1Norm[0])/(spanX.end[0]-spanX.start[0]);
									spanX.NYslope = (spanX.vert2Norm[1]-spanX.vert1Norm[1])/(spanX.end[0]-spanX.start[0]);
									spanX.NZslope = (spanX.vert2Norm[2]-spanX.vert1Norm[2])/(spanX.end[0]-spanX.start[0]);
								}

								spanX.Ustart = Lverts[1]->Ucurrent;
								spanX.Vstart = Lverts[1]->Vcurrent;
								spanX.Uend = Rverts[0]->Ucurrent;
								spanX.Vend = Rverts[0]->Vcurrent;
								spanX.Ucurrent = Lverts[1]->Ucurrent;
								spanX.Vcurrent = Lverts[1]->Vcurrent;
								spanX.Uslope = (spanX.Uend-spanX.Ustart)/(spanX.end[0]-spanX.start[0]);
								spanX.Vslope = (spanX.Vend-spanX.Vstart)/(spanX.end[0]-spanX.start[0]);

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
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0] != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								
							}
							else if((yLine == midY+1 && midY != MaxY)){
								
								dY = midY - Lverts[1]->start[1]+1;
								Lverts[1]->current[0] = Lverts[1]->current[0]+(Lverts[1]->Xslope*dY);
								Lverts[1]->current[1] = Lverts[1]->current[1]+dY;
								Lverts[1]->current[2] = Lverts[1]->current[2]+(Lverts[1]->Zslope*dY);
								Rverts[0]->current[0] = Rverts[0]->current[0]+(Rverts[0]->Xslope*1);
								Rverts[0]->current[1] = Rverts[0]->current[1]+1;
								Rverts[0]->current[2] = Rverts[0]->current[2]+(Rverts[0]->Zslope*1);

								if(render->interp_mode == GZ_COLOR){
									Lverts[1]->currentColor[0] = Lverts[1]->currentColor[0]+(Lverts[1]->Rslope*dY);
									Lverts[1]->currentColor[1] = Lverts[1]->currentColor[1]+(Lverts[1]->Gslope*dY);
									Lverts[1]->currentColor[2] = Lverts[1]->currentColor[2]+(Lverts[1]->Bslope*dY);
									Rverts[0]->currentColor[0] = Rverts[0]->currentColor[0]+(Rverts[0]->Rslope*1);
									Rverts[0]->currentColor[1] = Rverts[0]->currentColor[1]+(Rverts[0]->Gslope*1);
									Rverts[0]->currentColor[2] = Rverts[0]->currentColor[2]+(Rverts[0]->Bslope*1);
								}
								else if(render->interp_mode == GZ_NORMALS){
									Lverts[1]->currentNorm[0] = Lverts[1]->currentNorm[0]+(Lverts[1]->NXslope*dY);
									Lverts[1]->currentNorm[1] = Lverts[1]->currentNorm[1]+(Lverts[1]->NYslope*dY);
									Lverts[1]->currentNorm[2] = Lverts[1]->currentNorm[2]+(Lverts[1]->NZslope*dY);
									Rverts[0]->currentNorm[0] = Rverts[0]->currentNorm[0]+(Rverts[0]->NXslope*1);
									Rverts[0]->currentNorm[1] = Rverts[0]->currentNorm[1]+(Rverts[0]->NYslope*1);
									Rverts[0]->currentNorm[2] = Rverts[0]->currentNorm[2]+(Rverts[0]->NZslope*1);
								}

								Lverts[1]->Ucurrent = Lverts[1]->Ucurrent + (Lverts[1]->Uslope*dY);
								Lverts[1]->Vcurrent = Lverts[1]->Vcurrent + (Lverts[1]->Vslope*dY);
								Rverts[0]->Ucurrent = Rverts[0]->Ucurrent + (Rverts[0]->Uslope*1);
								Rverts[0]->Vcurrent = Rverts[0]->Vcurrent + (Rverts[0]->Vslope*1);

								//set up span DDA
								minX = Lverts[1]->current[0];
								MaxX = Rverts[0]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[1]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[0]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[1]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);

								if(render->interp_mode == GZ_COLOR){
									memcpy(spanX.vert1Color,Lverts[1]->currentColor,3*sizeof(float));
									memcpy(spanX.vert2Color,Rverts[0]->currentColor,3*sizeof(float));
									memcpy(spanX.currentColor,Lverts[1]->currentColor,3*sizeof(float));
									spanX.Rslope = (spanX.vert2Color[0]-spanX.vert1Color[0])/(spanX.end[0]-spanX.start[0]);
									spanX.Gslope = (spanX.vert2Color[1]-spanX.vert1Color[1])/(spanX.end[0]-spanX.start[0]);
									spanX.Bslope = (spanX.vert2Color[2]-spanX.vert1Color[2])/(spanX.end[0]-spanX.start[0]);
								}
								else if(render->interp_mode == GZ_NORMALS){
									memcpy(spanX.vert1Norm,Lverts[1]->currentNorm,3*sizeof(float));
									memcpy(spanX.vert2Norm,Rverts[0]->currentNorm,3*sizeof(float));
									memcpy(spanX.currentNorm,Lverts[1]->currentNorm,3*sizeof(float));
									spanX.NXslope = (spanX.vert2Norm[0]-spanX.vert1Norm[0])/(spanX.end[0]-spanX.start[0]);
									spanX.NYslope = (spanX.vert2Norm[1]-spanX.vert1Norm[1])/(spanX.end[0]-spanX.start[0]);
									spanX.NZslope = (spanX.vert2Norm[2]-spanX.vert1Norm[2])/(spanX.end[0]-spanX.start[0]);
								}

								spanX.Ustart = Lverts[1]->Ucurrent;
								spanX.Vstart = Lverts[1]->Vcurrent;
								spanX.Uend = Rverts[0]->Ucurrent;
								spanX.Vend = Rverts[0]->Vcurrent;
								spanX.Ucurrent = Lverts[1]->Ucurrent;
								spanX.Vcurrent = Lverts[1]->Vcurrent;
								spanX.Uslope = (spanX.Uend-spanX.Ustart)/(spanX.end[0]-spanX.start[0]);
								spanX.Vslope = (spanX.Vend-spanX.Vstart)/(spanX.end[0]-spanX.start[0]);

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
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0] != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								
							}
							else if(midY != MaxY && yLine > midY){
								dY = 1;
								Lverts[1]->current[0] = Lverts[1]->current[0]+(Lverts[1]->Xslope*dY);
								Lverts[1]->current[1] = Lverts[1]->current[1]+dY;
								Lverts[1]->current[2] = Lverts[1]->current[2]+(Lverts[1]->Zslope*dY);
								Rverts[0]->current[0] = Rverts[0]->current[0]+(Rverts[0]->Xslope*dY);
								Rverts[0]->current[1] = Rverts[0]->current[1]+dY;
								Rverts[0]->current[2] = Rverts[0]->current[2]+(Rverts[0]->Zslope*dY);

								if(render->interp_mode == GZ_COLOR){
									Lverts[1]->currentColor[0] = Lverts[1]->currentColor[0]+(Lverts[1]->Rslope*dY);
									Lverts[1]->currentColor[1] = Lverts[1]->currentColor[1]+(Lverts[1]->Gslope*dY);
									Lverts[1]->currentColor[2] = Lverts[1]->currentColor[2]+(Lverts[1]->Bslope*dY);
									Rverts[0]->currentColor[0] = Rverts[0]->currentColor[0]+(Rverts[0]->Rslope*dY);
									Rverts[0]->currentColor[1] = Rverts[0]->currentColor[1]+(Rverts[0]->Gslope*dY);
									Rverts[0]->currentColor[2] = Rverts[0]->currentColor[2]+(Rverts[0]->Bslope*dY);
								}
								else if(render->interp_mode == GZ_NORMALS){
									Lverts[1]->currentNorm[0] = Lverts[1]->currentNorm[0]+(Lverts[1]->NXslope*dY);
									Lverts[1]->currentNorm[1] = Lverts[1]->currentNorm[1]+(Lverts[1]->NYslope*dY);
									Lverts[1]->currentNorm[2] = Lverts[1]->currentNorm[2]+(Lverts[1]->NZslope*dY);
									Rverts[0]->currentNorm[0] = Rverts[0]->currentNorm[0]+(Rverts[0]->NXslope*dY);
									Rverts[0]->currentNorm[1] = Rverts[0]->currentNorm[1]+(Rverts[0]->NYslope*dY);
									Rverts[0]->currentNorm[2] = Rverts[0]->currentNorm[2]+(Rverts[0]->NZslope*dY);
								}

								Lverts[1]->Ucurrent = Lverts[1]->Ucurrent + (Lverts[1]->Uslope*dY);
								Lverts[1]->Vcurrent = Lverts[1]->Vcurrent + (Lverts[1]->Vslope*dY);
								Rverts[0]->Ucurrent = Rverts[0]->Ucurrent + (Rverts[0]->Uslope*dY);
								Rverts[0]->Vcurrent = Rverts[0]->Vcurrent + (Rverts[0]->Vslope*dY);

								//set up span DDA
								minX = Lverts[1]->current[0];
								MaxX = Rverts[0]->current[0];
								SPANDDA spanX;
								memcpy(spanX.start,Lverts[1]->current,3*sizeof(float));
								memcpy(spanX.end,Rverts[0]->current,3*sizeof(float));
								memcpy(spanX.current,Lverts[1]->current,3*sizeof(float));
								spanX.Zslope = (spanX.end[2]-spanX.start[2])/(spanX.end[0]-spanX.start[0]);

								if(render->interp_mode == GZ_COLOR){
									memcpy(spanX.vert1Color,Lverts[1]->currentColor,3*sizeof(float));
									memcpy(spanX.vert2Color,Rverts[0]->currentColor,3*sizeof(float));
									memcpy(spanX.currentColor,Lverts[1]->currentColor,3*sizeof(float));
									spanX.Rslope = (spanX.vert2Color[0]-spanX.vert1Color[0])/(spanX.end[0]-spanX.start[0]);
									spanX.Gslope = (spanX.vert2Color[1]-spanX.vert1Color[1])/(spanX.end[0]-spanX.start[0]);
									spanX.Bslope = (spanX.vert2Color[2]-spanX.vert1Color[2])/(spanX.end[0]-spanX.start[0]);
								}
								else if(render->interp_mode == GZ_NORMALS){
									memcpy(spanX.vert1Norm,Lverts[1]->currentNorm,3*sizeof(float));
									memcpy(spanX.vert2Norm,Rverts[0]->currentNorm,3*sizeof(float));
									memcpy(spanX.currentNorm,Lverts[1]->currentNorm,3*sizeof(float));
									spanX.NXslope = (spanX.vert2Norm[0]-spanX.vert1Norm[0])/(spanX.end[0]-spanX.start[0]);
									spanX.NYslope = (spanX.vert2Norm[1]-spanX.vert1Norm[1])/(spanX.end[0]-spanX.start[0]);
									spanX.NZslope = (spanX.vert2Norm[2]-spanX.vert1Norm[2])/(spanX.end[0]-spanX.start[0]);
								}

								spanX.Ustart = Lverts[1]->Ucurrent;
								spanX.Vstart = Lverts[1]->Vcurrent;
								spanX.Uend = Rverts[0]->Ucurrent;
								spanX.Vend = Rverts[0]->Vcurrent;
								spanX.Ucurrent = Lverts[1]->Ucurrent;
								spanX.Vcurrent = Lverts[1]->Vcurrent;
								spanX.Uslope = (spanX.Uend-spanX.Ustart)/(spanX.end[0]-spanX.start[0]);
								spanX.Vslope = (spanX.Vend-spanX.Vstart)/(spanX.end[0]-spanX.start[0]);

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
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine == minX+1 && spanX.current[0] != xLine){
											dX = minX+1 -spanX.start[0];
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
										else if(xLine != minX){
											dX = 1;
											spanX.current[0] = spanX.current[0] + dX;
											spanX.current[2] = spanX.current[2] + (dX*spanX.Zslope);

											if(render->interp_mode == GZ_COLOR){
												spanX.currentColor[0] = spanX.currentColor[0] + (spanX.Rslope*dX);
												spanX.currentColor[1] = spanX.currentColor[1] + (spanX.Gslope*dX);
												spanX.currentColor[2] = spanX.currentColor[2] + (spanX.Bslope*dX);
											}
											else if(render->interp_mode == GZ_NORMALS){
												spanX.currentNorm[0] = spanX.currentNorm[0] + (spanX.NXslope*dX);
												spanX.currentNorm[1] = spanX.currentNorm[1] + (spanX.NYslope*dX);
												spanX.currentNorm[2] = spanX.currentNorm[2] + (spanX.NZslope*dX);
											}

											spanX.Ucurrent = spanX.Ucurrent + (spanX.Uslope*dX);
											spanX.Vcurrent = spanX.Vcurrent + (spanX.Vslope*dX);

											float  intZ = interpZ(xLine,yLine,A,B,C,D);
											GzIntensity r,g,b,a;
											GzDepth z;
											GzGetDisplay(render->display,xLine,yLine,&r,&g,&b,&a,&z);
											if(intZ < z){
												
												switch(render->interp_mode){
												case(GZ_NONE):
													{
														r = ctoi(render->flatcolor[0]);
														g = ctoi(render->flatcolor[1]);
														b = ctoi(render->flatcolor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_COLOR):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														r = ctoi(spanX.currentColor[0]*Kt[0]);
														g = ctoi(spanX.currentColor[1]*Kt[1]);
														b = ctoi(spanX.currentColor[2]*Kt[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												case(GZ_NORMALS):
													{
														GzColor Kt;
														//transform interpolated uv
														float tu = spanX.Ucurrent*((intZ/(MAXINT-intZ))+1);
														float tv = spanX.Vcurrent*((intZ/(MAXINT-intZ))+1);
														render->tex_fun(tu,tv,Kt);
														memcpy(render->Ka,Kt,3*sizeof(float));
														memcpy(render->Kd,Kt,3*sizeof(float));
														GzColor resColor;
														calculateColor(render,spanX.currentNorm,resColor);
														r = ctoi(resColor[0]);
														g = ctoi(resColor[1]);
														b = ctoi(resColor[2]);
														a = 4095;
														z = intZ;
														GzPutDisplay(render->display,xLine,yLine,r,g,b,a,z);
														break;
													}
												}
											}
											else{
												//skip the pixel
											}
										}
									}
								}
								
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


