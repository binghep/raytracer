#include <stdio.h>
#include <stdlib.h>
#include <iostream>   //std::cout
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>     //asin   (arcsin)     abs for float is only in cmath rather than math.h
#include <algorithm>    // std::min

using namespace std;
#include "assg.h"

void printVec3(Vec3 v){
	printf("Vec3:  x = %f, y = %f, z = %f.\n", v.x, v.y, v.z);

}

void printPoint(Point p){
	printf("Point:  x = %f, y = %f, z = %f.\n", p.x, p.y, p.z);
}
void printcolor(ColorType p){
	printf("Color:  r = %f, g = %f, b = %f.\n", p.r, p.g, p.b);
}
ColorType *c = new ColorType(0.0, 0.0, 0.0);

RayType *j = new RayType();
Vec3 *raydir = new Vec3();
Point *intersection = new Point();
RayType *reflected_ray = new RayType();
Point *from = new Point();
Point *nextfrom = new Point();


RayType *trans = new RayType();
Point *outpoint = new Point();
Vec3 *T = new Vec3();
Vec3 *N1 = new Vec3();


Vec3 *I = new Vec3();
Vec3 *N = new Vec3();
Vec3 *V = new Vec3();
Vec3 *H = new Vec3();
Vec3 *R = new Vec3();
RayType *makeray = new RayType();
ColorType *return_color = new ColorType();

Point eye;
ColorType bkgcolor;
MaterialColor materialcolor;
Vec3 viewdir, updir;
float viewdist, fovv, fovh, aspect;
int pixheight, pixwidth, numsphere, numspherecopy, numlights;
SphereType spheres[MAXSPHERE];
LightType lights[MAXLIGHT];
Vec3 L[MAXLIGHT];//opposite direction to the incomming light
Node list[MAXSPHERE];

float Fr, R0;

Point *f = new Point();

Vec3 *d = new Vec3();


/*return the index of the sphere the point belongs to
return -1 if not belong to any sphere*/
int belongto(Point *p){
	int whichsphere = -1;
	float smallest = -1;
	float small_error = 0;
	for (int i = 0; i < numsphere; ++i){
		small_error = abs((p->x - spheres[i].x)*(p->x - spheres[i].x) + (p->y - spheres[i].y)*(p->y - spheres[i].y) + (p->z - spheres[i].z)*(p->z - spheres[i].z) - spheres[i].r*spheres[i].r);
		if (small_error < 0.01){
			if ((smallest == -1) || (smallest>small_error)) {
				smallest = small_error;
				whichsphere = i;
			}
		}
	}
	//	cout << smallest<<endl;
	return whichsphere;
}

//find intersection of a ray and a sphere. construct node for each sphere intersected in "list" array. 
//In node, numintersection can be 0, 1, 2. only positive distance count.
//return the shortest positive distance--closest, -1.0 if no intersection 
float look(RayType *r, int selfexclusive_0){
	f->x = r->x;
	f->y = r->y;
	f->z = r->z;
	d->x = r->dx;
	d->y = r->dy;
	d->z = r->dz;
	for (int i = 0; i < numsphere; ++i){
		list[i].n = 0;
		list[i].p1 = 0;
		list[i].p2 = 0;
	}
	/*construct "list" array*/
	int whichsphere = belongto(f);  //selfexclusive_0  is 0 in shadow function. don't consider current sphere
	for (int i = 0; i < numsphere; ++i){
		if (selfexclusive_0 == 0) { if (whichsphere == i) continue; }
		float B = 2 * (r->dx*(r->x - spheres[i].x) + r->dy*(r->y - spheres[i].y) + r->dz*(r->z - spheres[i].z));
		float C = (r->x - spheres[i].x)*(r->x - spheres[i].x) + (r->y - spheres[i].y)*(r->y - spheres[i].y) + (r->z - spheres[i].z)*(r->z - spheres[i].z) - spheres[i].r*spheres[i].r;
		float discriminant = B*B - 4 * C;
		if (discriminant < 0) {
			list[i].n = 0;
		}
		else if (discriminant == 0)//does not mean has intersection with this ray which starts from a point.
		{
			float t = -B / 2;
			if (0.02 < t) {
				list[i].n = 1;
				list[i].p1 = t;
			}
			else {
				//cout << "has 1 intersection if trace backward" << endl;
			}
		}
		else {
			float t1 = (-B + sqrt(discriminant)) / 2;
			float t2 = (-B - sqrt(discriminant)) / 2;
			if (0.02 < t1&&t2>0.02&&t1 > t2)
			{
				list[i].n = 2;
				list[i].p1 = t2;
				list[i].p2 = t1;
			}
			else if (0.02 < t1&&t2>0.02&&t1 < t2)
			{
				list[i].n = 2;
				list[i].p1 = t1;
				list[i].p2 = t2;
			}
			else if (0.02 < t1&&t2<0.02)
			{
				list[i].n = 1;
				list[i].p1 = t1;
			}
			else if (0.02 < t1&&t2 > 0.02)
			{
				list[i].n = 1;
				list[i].p1 = t2;
			}
			else if (0.02 > t1&&t2 < 0.02){
				list[i].n = 0;
			}
			else {
				//cout << "has 2 intersection if trace backward" << endl;
			}
		}
	}
	float closest = -1.0;
	for (int i = 0; i < numsphere; ++i){
		if (list[i].n == 0) continue;
		else if (list[i].n == 1){
			if ((list[i].p1 < closest) || (closest == -1.0)) { closest = list[i].p1; }
		}
		else if (list[i].n == 2){
			if ((list[i].p1 < closest) || (closest == -1.0)) { closest = list[i].p1; }
		}
		else {
			cout << "shouldn't" << endl;
		}
	}
	if (closest != -1.0) { return closest; }// *intersection = (*f) + (*d)*closest; 
	else return -1.0;//no intersection with any sphere
}


ColorType *h = new ColorType();

ColorType Shade_Ray(Point intersec, int whichsphere)// two arguments tested correct later.
{
	Point interseccopy = intersec;
	/*------------components needed for the 2nd, 3rd part, and the first step of reflected and transmitted ray tracing part---------*/
	*N = (intersec - spheres[whichsphere].center) / (spheres[whichsphere].r);
	*V = (eye - intersec) / ((eye - intersec).GetMagnitude());
	*I = *V;
	/*reflected ray of I*/
	*R = (*N) * 2 * ((*N)&(*I)) - *I;
	*R = (*R) / (*R).GetMagnitude();
	/*------this reflected ray pointer and it's direction(Vec3 *R) are used in every call of recursive reflective ray tracing-------*/
	reflected_ray->dx = R->x;
	reflected_ray->dy = R->y;
	reflected_ray->dz = R->z;
	reflected_ray->x = intersec.x;
	reflected_ray->y = intersec.y;
	reflected_ray->z = intersec.z;
	/*----------At this intersection point, putting the NdotH of each light in an array sequantially--------------------------*/
	float NdotH[MAXLIGHT];
	float NdotH_to_nth[MAXLIGHT];

	for (int p = 0; p < MAXLIGHT; ++p){
		NdotH_to_nth[p] = 1;//later we can multiply 1 by NdotH n times
	}

	for (int i = 0; i < numlights; ++i){
		if (lights[i].w == 1){//if point light, update corresponding element in L[i]
			L[i] = (lights[i].lightlocation - intersec) / ((lights[i].lightlocation - intersec).GetMagnitude());
		}
		float NdotL = (*N)&(L[i]);//cos(theta)
		/*if (NdotL < 0), then angle between N and L exceeds 90 degrees, the light is hitting the backside of the surface*/
		NdotL = max((float)0, NdotL);

		/*------------get N&H(i) for both point and directional light-------------------------------------------------------------*/
		*H = (L[i] + *V) / (L[i] + *V).GetMagnitude();
		NdotH[i] = max((float)0, (*N)&(*H));

		//get N&H_to_the_nth(2)
		for (int j = 0; j < materialcolor.n; ++j){//add specular terms
			NdotH_to_nth[i] *= NdotH[i];
		}
	}
	/*--------------------adding the first part of the Phong Model------------------------------------------------------------*/
	return_color->r = materialcolor.ka*materialcolor.dr;
	return_color->g = materialcolor.ka*materialcolor.dg;
	return_color->b = materialcolor.ka*materialcolor.db;
	/*--------------------adding the second and the third part of the Phong Model---------------------------------------------*/
	for (int i = 0; i < numlights; ++i){
		//shooting ray from the intersection to the light
		makeray->x = intersec.x;
		makeray->y = intersec.y;
		makeray->z = intersec.z;

		makeray->dx = L[i].x;
		makeray->dy = L[i].y;
		makeray->dz = L[i].z;
		float distance = -1;/*-1 if directional light*/
		if (lights[i].w == 1){//update distance if point light
			Vec3 distance_to_light = lights[i].lightlocation - intersec;
			distance = distance_to_light.GetMagnitude();
		}
		float shadowflag = -1.0;

		shadowflag = shadow(makeray, distance, whichsphere);//shoot ray to this light, find out whether shadow is formed.
		//cout << shadowflag;
		//if (shadowflag == 0) { cout << "shadowflag:   " << shadowflag << endl; }
		if (shadowflag == -1.0){ cout << "Error: shadowflag doesn't get value/n"; }
		return_color->r += shadowflag*lights[i].r*(materialcolor.kd * materialcolor.dr*((*N)&L[i]) + materialcolor.ks*materialcolor.sr*NdotH_to_nth[i]);
		return_color->g += shadowflag*lights[i].g*(materialcolor.kd * materialcolor.dg*((*N)&L[i]) + materialcolor.ks*materialcolor.sg*NdotH_to_nth[i]);
		return_color->b += shadowflag*lights[i].b*(materialcolor.kd * materialcolor.db*((*N)&L[i]) + materialcolor.ks*materialcolor.sb*NdotH_to_nth[i]);
	}
	/*--------------------------adding recursive reflected ray tracing part-------------------------------------------------*/
	R0 = ((materialcolor.eta - 1) / (1 + materialcolor.eta))*((materialcolor.eta - 1) / (1 + materialcolor.eta));//0.04
	Fr = R0 + (1 - R0)*(((*I)&(*N))*(-1) + 1)*(((*I)&(*N))*(-1) + 1)*(((*I)&(*N))*(-1) + 1)*(((*I)&(*N))*(-1) + 1)*(((*I)&(*N))*(-1) + 1);
	// cout <<R0<<"          "<< Fr << endl;
	ColorType color;
	while ((abs((*N)&(*R) - 0) > 0.001)&&(Fr > 0.01)&&((color = reflective_trace_and_shade(reflected_ray)) != bkgcolor))
	{
		return_color->r += color.r*Fr;
		return_color->g += color.g*Fr;
		return_color->b += color.b*Fr;
		Fr = Fr*Fr;
	}
	if (color == bkgcolor){
		return_color->r += color.r*Fr;
		return_color->g += color.g*Fr;
		return_color->b += color.b*Fr;
	}
	/*-----------------------adding transmitted ray tracing part--------------------------------------------------*/
	/*-------recover the value of *N and *I we got at the beginning of this function----------------------------------------*/
	*N = (intersec - spheres[whichsphere].center) / (spheres[whichsphere].r);
	*V = (eye - intersec) / ((eye - intersec).GetMagnitude());
	*I = *V;
	//if (intersec == interseccopy) {} else {cout << " dddd" << endl;}
	*T = ((*N)*(-1))*sqrt(1 - ((1 - ((*N)&(*I))*((*N)&(*I)))*(1 / materialcolor.eta)*(1 / materialcolor.eta))) + ((*N)*((*N)&(*I)) - *I) / materialcolor.eta;
	*T = (*T) / T->GetMagnitude();
	//if (belongto(&interseccopy) != whichspherecopy) { 
	//cout << belongto(&interseccopy) << "   " << whichspherecopy << endl; //}
	trans->dx = T->x;
	trans->dy = T->y;
	trans->dz = T->z;
	trans->x = interseccopy.x;
	trans->y = interseccopy.y;
	trans->z = interseccopy.z;

	add_two_colors_due_to_this_ray_in_sphere(trans, whichsphere);//returned by ch, ch2, not recursive
	//printcolor(*h);
	//printcolor(bkgcolor);

	/*---------------------------------------------finish-----------------------------------------------------------------*/
	return_color->r = min((float)1.0, return_color->r);
	return_color->g = min((float)1.0, return_color->g);
	return_color->b = min((float)1.0, return_color->b);
	return(*return_color);
}


void add_two_colors_due_to_this_ray_in_sphere(RayType *trans, int whichsphere){
	float angle = 57.296*asin(1 / materialcolor.eta);//degrees
	h->r = 0;
	h->g = 0;
	h->b = 0;
	s g(trans);
	g.child(whichsphere);

	float u = 0;
	if ((u = look(&(g.ch), 0)) != -1.0){
		/*-------------------------------consider total reflection-------------------------------------------------------------*/
		
		if (acos(g.beforeout)*57.296 < angle){
			pointcolor(&((g.out + (g.d)*u)));
			h->r = h->r*g.transmitconstant;
			h->g = h->g*g.transmitconstant;
			h->b = h->b*g.transmitconstant;
			return_color->r += h->r;
			return_color->g += h->g;
			return_color->b += h->b;
			h->r = 0;
			h->g = 0;
			h->b = 0;
			s ggg(&g.ch);
			ggg.child(whichsphere);
			float w = 0;
			if ((w = look(&(ggg.ch), 0)) != -1.0){
				if (acos(ggg.beforeout)*57.296 < angle){
					pointcolor(&((ggg.out + (ggg.d)*w)));
					h->r = h->r*ggg.transmitconstant;
					h->g = h->g*ggg.transmitconstant;
					h->b = h->b*ggg.transmitconstant;
					return_color->r += h->r;
					return_color->g += h->g;
					return_color->b += h->b;
					ColorType a = bkgcolor*ggg.transmitconstant;
					h->r = a.r;
					h->g = a.g;
					h->b = a.b;
					return_color->r += h->r;
					return_color->g += h->g;
					return_color->b += h->b;
				}
			}
			else {
				ColorType a = bkgcolor*g.transmitconstant;
				h->r = a.r;
				h->g = a.g;
				h->b = a.b;
				return_color->r += h->r;
				return_color->g += h->g;
				return_color->b += h->b;
			}
		}
		else {
			ColorType a = bkgcolor*g.transmitconstant;
			h->r = a.r;
			h->g = a.g;
			h->b = a.b;
			return_color->r += h->r;
			return_color->g += h->g;
			return_color->b += h->b;
		}

	}
/*--------------inner reflect--------------------------------------------------*/
	h->r = 0;
	h->g = 0;
	h->b = 0;
	RayType *inner_reflect = &(g.ch2);
	s gg(inner_reflect);
	gg.child(whichsphere);
	u = 0;
	if ((u = look(&(gg.ch), 0)) != -1.0){
		if (acos(gg.beforeout)*57.296 < angle){
			pointcolor(&(gg.out + (gg.d)*u));
			h->r = h->r*gg.transmitconstant*g.Frout;
			h->g = h->g*gg.transmitconstant*g.Frout;
			h->b = h->b*gg.transmitconstant*g.Frout;
			return_color->r += h->r;
			return_color->g += h->g;
			return_color->b += h->b;
			h->r = 0;
			h->g = 0;
			h->b = 0;
			s ggg(&gg.ch2);
			ggg.child(whichsphere);
			float w = 0;
			if ((w = look(&(ggg.ch), 0)) != -1.0){
				if (acos(ggg.beforeout)*57.296 < angle){
					pointcolor(&((ggg.out + (ggg.d)*w)));
					h->r = h->r*ggg.transmitconstant;
					h->g = h->g*ggg.transmitconstant;
					h->b = h->b*ggg.transmitconstant;
					return_color->r += h->r;
					return_color->g += h->g;
					return_color->b += h->b;
					ColorType a = bkgcolor*ggg.transmitconstant;
					return_color->r += a.r;
					return_color->g += a.g;
					return_color->b += a.b;
				}
			}
			else {
				ColorType a = bkgcolor*ggg.transmitconstant;
				return_color->r += a.r;
				return_color->g += a.g;
				return_color->b += a.b;
			}
		}

	}
	else {
		ColorType a = bkgcolor*gg.transmitconstant*g.Frout;
		return_color->r += a.r;
		return_color->g += a.g;
		return_color->b += a.b;
	}
}

/*
class s{
RayType parent, ch, ch2;
float beforein, afterin, beforeout, afterout,Frin, Frout;
Point in, out;
Vec3 N, N1, T;

s(RayType* X, float Y){
parent = *X;
Frin = Y;
}
void child();
};
*/
void s::print(){
	printVec3(d);
	printVec3(d2);
	printPoint(in);
	printPoint(out);
}


bool s::child(int whichsphere){//false if the point is on a transmited ray that tangent to the sphere. don't need transmit and reflected color.
	in.x = parent.x;
	in.y = parent.y;
	in.z = parent.z;
	T.x = parent.dx;
	T.y = parent.dy;
	T.z = parent.dz;

	float B = 2 * (parent.dx*(parent.x - spheres[whichsphere].x) + parent.dy*(parent.y - spheres[whichsphere].y) + parent.dz*(parent.z - spheres[whichsphere].z));
	float C = (parent.x - spheres[whichsphere].x)*(parent.x - spheres[whichsphere].x) + (parent.y - spheres[whichsphere].y)*(parent.y - spheres[whichsphere].y) + (parent.z - spheres[whichsphere].z)*(parent.z - spheres[whichsphere].z) - spheres[whichsphere].r*spheres[whichsphere].r;
	float discriminant = B*B - 4 * C;
	//cout << discriminant << endl;
	if (discriminant <= 0){ h = c; return false; }
	//else if ((abs(discriminant))<0.01) { cout << "Er..." << endl; }//return false; }
	float t1 = (-B + sqrt(discriminant)) / 2;
	float t2 = (-B - sqrt(discriminant)) / 2;
	//if (t1>0 && t2>0) cout << 2 << endl;
	//if (t1 < 0 && t2 < 0) cout << 0 << endl;
	//if (t1>0 && t2<0) cout << 1 << endl;
	//if (t1 < 0 && t2 > 0) cout << 1 << endl;
	if (t2 > t1) 	{ out = in + T*t2; } //cout <<t2<< endl;}
	else  { out = in + T*t1; } //cout << t1 << endl;}
	//N1 is the surface normal at outpoint
	N1 = (out - spheres[whichsphere].center) / (out - spheres[whichsphere].center).GetMagnitude();
	beforeout = (N1*(-1))&(T*(-1));
	/*--------------------------------internal reflected ch2 with direction d2---------------------------*/
	d2 = (N1*(-1)) * 2 * ((N1*(-1))&(T*(-1))) - T*(-1);
	d2 = d2 / d2.GetMagnitude();
	ch2.x = out.x;
	ch2.y = out.y;
	ch2.z = out.z;
	ch2.dx = d2.x;
	ch2.dx = d2.y;
	ch2.dx = d2.z;
	/*----------------------------transmited ch with direction d----------------------------------------*/
	/*-----------------Fr, afterout--R0 is same for each ray/sphere point-------------*/
	d = N1*sqrt(1 - ((1 - beforeout*beforeout)*(materialcolor.eta / 1)*(materialcolor.eta / 1))) + (N1*(-1)*beforeout - T*(-1)) * materialcolor.eta;
	d = d / d.GetMagnitude();
	afterout = N1&d;
	Frout = R0 + (1 - R0)*(1 - afterout)*(1 - afterout)*(1 - afterout)*(1 - afterout)*(1 - afterout);//using the afterout(cos of transmitted angle) because it's bigger.
	transmitconstant = (1 - Frout)*(1 - materialcolor.alpha);
	ch.dx = d.x;
	ch.dy = d.y;
	ch.dz = d.z;
	ch.x = out.x;
	ch.y = out.y;
	ch.z = out.z;

	return true;
}
/*shade this point(verifed on sphere) using Phong model except adding recursive transmitted and reflected ray tracing part*/
void pointcolor(Point *p){ //modifies *N,*V,*I,*H, *makeray
	h->r = materialcolor.ka*materialcolor.dr;
	h->g = materialcolor.ka*materialcolor.dg;
	h->b = materialcolor.ka*materialcolor.db;
	int whichsphere = belongto(p);
	//	if (whichsphere == -1.0) { cout << "can't find whichsphere the ray start from." <<"            "<<"exit"<< endl; }
	*N = (*p - spheres[whichsphere].center) / (spheres[whichsphere].r);
	*V = (eye - *p) / ((eye - *p).GetMagnitude());
	*I = *V;
	float NdotH[MAXLIGHT];//for each light
	float NdotH_to_nth[MAXLIGHT];//for each light
	for (int a = 0; a < MAXLIGHT; ++a){
		NdotH_to_nth[a] = 1;//later we can multipy n NdotH
	}
	for (int i = 0; i < numlights; ++i){
		if (lights[i].w == 1){//point light
			L[i] = (lights[i].lightlocation - *p) / ((lights[i].lightlocation - *p).GetMagnitude());
		}
		float NdotL = (*N)&(L[i]);//cos(theta)   if (NdotL < 0) :angle between N and L exceeds 90 degrees, the light is hitting the backside of the surface
		NdotL = max((float)0, NdotL);

		//get N&H(i)
		*H = (L[i] + *V) / (L[i] + (*V)).GetMagnitude();//for either kind of light	
		NdotH[i] = max((float)0, (*N)&(*H));

		//get N&H_to_the_nth(2)
		for (int j = 0; j < materialcolor.n; ++j){//add specular terms
			NdotH_to_nth[i] *= NdotH[i];
		}
	}
	for (int i = 0; i < numlights; ++i){
		//shooting ray from the intersection to the light
		makeray->x = p->x;
		makeray->y = p->y;
		makeray->z = p->z;
		makeray->dx = L[i].x;
		makeray->dy = L[i].y;
		makeray->dz = L[i].z;
		float distance = -1;/*-1 if directional light*/
		if (lights[i].w == 1){//update distance if point light
			Vec3 distance_to_light = lights[i].lightlocation - *p;
			distance = distance_to_light.GetMagnitude();
		}
		float shadowflag = -1.0;
		shadowflag = shadow(makeray, distance, whichsphere);//shoot ray to this light, find out whether shadow is formed.
		if (shadowflag == -1.0){ cout << "Error: shadowflag doesn't get value/n"; exit(4); }
		h->r += shadowflag*lights[i].r*(materialcolor.kd * materialcolor.dr*((*N)&L[i]) + materialcolor.ks*materialcolor.sr*NdotH_to_nth[i]);
		h->g += shadowflag*lights[i].g*(materialcolor.kd * materialcolor.dg*((*N)&L[i]) + materialcolor.ks*materialcolor.sg*NdotH_to_nth[i]);
		h->b += shadowflag*lights[i].b*(materialcolor.kd * materialcolor.db*((*N)&L[i]) + materialcolor.ks*materialcolor.sb*NdotH_to_nth[i]);
	}
	//	return h;
}



ColorType reflective_trace_and_shade(RayType *ray){//ray is reflected ray
	ColorType color;
	from->x = ray->x;
	from->y = ray->y;
	from->z = ray->z;
	raydir->x = ray->dx;
	raydir->y = ray->dy;
	raydir->z = ray->dz;
	/*find out whether there is intersection on  the reflected ray*/
	float	closest = look(ray, 0);
	/*no intersection*/
	if (closest == -1.0){
		return bkgcolor;
	}
	/*intersection*/
	*nextfrom = *from + (*raydir)*closest;
	/*the intersection with which sphere?*/
	int whichsphere = belongto(nextfrom);

	*N = (*nextfrom - spheres[whichsphere].center) / (spheres[whichsphere].r);
	*V = (eye - *nextfrom) / ((eye - *nextfrom).GetMagnitude());
	*I = (*raydir)*(-1);


	float NdotH[MAXLIGHT];//for each light
	float NdotH_to_nth[MAXLIGHT];//for each light

	for (int a = 0; a < MAXLIGHT; ++a){
		NdotH_to_nth[a] = 1;//later we can multipy n NdotH
	}

	for (int i = 0; i < numlights; ++i){
		if (lights[i].w == 1){//point light
			L[i] = (lights[i].lightlocation - *nextfrom) / ((lights[i].lightlocation - *nextfrom).GetMagnitude());
		}
		float NdotL = (*N)&(L[i]);//cos(theta)
		//if (NdotL < 0){//angle between N and L exceeds 90 degrees, the light is hitting the backside of the surface
		NdotL = max((float)0, NdotL);

		//get N&H(i)
		*H = (L[i] + *V) / (L[i] + (*V)).GetMagnitude();//for either kind of light	
		NdotH[i] = max((float)0, (*N)&(*H));

		//get N&H_to_the_nth(2)
		for (int j = 0; j < materialcolor.n; ++j){//add specular terms
			NdotH_to_nth[i] *= NdotH[i];
		}
	}
	color.r = materialcolor.ka*materialcolor.dr;
	color.g = materialcolor.ka*materialcolor.dg;
	color.b = materialcolor.ka*materialcolor.db;

	for (int i = 0; i < numlights; ++i){
		//shooting ray from the intersection to the light
		makeray->x = nextfrom->x;
		makeray->y = nextfrom->y;
		makeray->z = nextfrom->z;
		makeray->dx = L[i].x;
		makeray->dy = L[i].y;
		makeray->dz = L[i].z;
		float distance = -1;/*-1 if directional light*/
		if (lights[i].w == 1){//update distance if point light
			Vec3 distance_to_light = lights[i].lightlocation - *nextfrom;
			distance = distance_to_light.GetMagnitude();
		}
		float shadowflag = -1.0;
		int whichsphere = belongto(nextfrom);
		if (whichsphere != -1.0){ shadowflag = shadow(makeray, distance, whichsphere); }//shoot ray to this light, find out whether shadow is formed.
		else { cout << "can't find a sphere that has this point" << endl; exit(3); }
		if (shadowflag == -1.0){ cout << "Error: shadowflag doesn't get value/n"; exit(4); }
		color.r += shadowflag*lights[i].r*(materialcolor.kd * materialcolor.dr*((*N)&L[i]) + materialcolor.ks*materialcolor.sr*NdotH_to_nth[i]);
		color.g += shadowflag*lights[i].g*(materialcolor.kd * materialcolor.dg*((*N)&L[i]) + materialcolor.ks*materialcolor.sg*NdotH_to_nth[i]);
		color.b += shadowflag*lights[i].b*(materialcolor.kd * materialcolor.db*((*N)&L[i]) + materialcolor.ks*materialcolor.sb*NdotH_to_nth[i]);
	}

	/*reflected ray of I*/
	*R = (*N) * 2 * ((*N)&(*I)) - *I;
	*R = (*R) / (*R).GetMagnitude();
	ray->dx = R->x;
	ray->dy = R->y;
	ray->dz = R->z;
	ray->x = nextfrom->x;
	ray->y = nextfrom->y;
	ray->z = nextfrom->z;
	return color;
}


/*shadowflag = shadow(makeray, L[i], distance, whichsphere);  shoot ray to this light, return shadow flag (0 is full shadow,1 is no shadow)*/
float shadow(RayType *ray, float distance, int whichsphere)
{
	int hitsphere = 0;//decide how dark the shadow is

	float result;

	if ((result = look(ray, 0)) != -1.0){//has intersection.
		//if directional light, then has shadow.
		for (int i = 0; i < numsphere; ++i){ if ((i != whichsphere)&&(list[i].n != 0)) ++hitsphere; }
		if (distance == -1.0){
			return max((double)0, (1.00 - hitsphere*materialcolor.alpha));//has shadow
		}
		//if point light
		if (result > distance){ return 1.0; }//no shadow
		else {
			//has shadow. return 0.0; if hit 3 spheres, then full shadow if hit 1 sphere, 30% shadow.
			return max((double)0, (1.00 - hitsphere*materialcolor.alpha));
		}
	}
	else {//no intersection on the ray direction
		return 1.0; //no  shadow
	}
}


//colormap[i][j]=Trace_Ray(ray[i][j], ray_dir[i][j]		
/*for each object in the scene, check for a
ray/object intersection; keep track of the
closest intersection point and of the identity
of the object hit at that point;
call Shade_Ray() to determine the color */
ColorType Trace_Ray(RayType ray, Vec3 ray_dir, Point from) //from is eye at first
{
	/* search for closest ray/object intersection point */
	/* call shade_ray to define the color at that point */
	ColorType return_color;
	float closest = -1.0;

	closest = look(&ray, 1);
	if (closest == -1){
		return_color = bkgcolor;
	}
	else{
		*intersection = from + ray_dir*closest;
		int whichsphere;
		if ((whichsphere = belongto(intersection)) != -1.0){ return_color = Shade_Ray(*intersection, whichsphere); }
		else { cout << "can't find whichsphere after calculating the intersection." << endl; exit(-2); }
	}
	return(return_color);
}


int main(void){
	/* read scene description from input file */
	ifstream source;                    // build a read(input)-Stream

	//source.open("input.txt", ios_base::in);  // open data
	//"C:\\temp\\in.txt"
	source.open("in.txt");
	if (!source)  {                     // if it does not work
		cerr << "Can't open Data! Abort! \n";
		return -1;
	}
	else {
		for (std::string line; std::getline(source, line);)   //read stream line by line
		{
			std::istringstream in(line);      //make a stream for the line itself

			std::string type;
			in >> type;                  //and read the first whitespace-separated token

			if (type == "eye")       //and check its value
			{
				float x, y, z;
				in >> x >> y >> z;       //now read the whitespace-separated floats
				eye.x = x;
				eye.y = y;
				eye.z = z;
			}
			else if (type == "viewdir")
			{
				float x, y, z;
				in >> x >> y >> z;
				viewdir.x = x;
				viewdir.y = y;
				viewdir.z = z;
			}
			else if (type == "updir")
			{
				float x, y, z;
				in >> x >> y >> z;
				updir.x = x;
				updir.y = y;
				updir.z = z;
			}
			else if (type == "bkgcolor")
			{
				float x, y, z;
				in >> x >> y >> z;
				bkgcolor.r = x;
				bkgcolor.g = y;
				bkgcolor.b = z;
			}
			else if (type == "materialcolor")
			{
				float dr, dg, db, sr, sg, sb, ka, kd, ks, n, alpha, eta;
				in >> dr >> dg >> db >> sr;
				in >> sg >> sb >> ka >> kd;
				in >> ks >> n >> alpha >> eta;
				materialcolor.dr = dr;
				materialcolor.dg = dg;
				materialcolor.db = db;
				materialcolor.sr = sr;
				materialcolor.sg = sg;
				materialcolor.sb = sb;
				materialcolor.ka = ka;
				materialcolor.kd = kd;
				materialcolor.ks = ks;
				materialcolor.n = n;
				materialcolor.alpha = alpha;
				materialcolor.eta = eta;
			}
			else if (type == "sphere")
			{
				float x, y, z, r;
				in >> x >> y >> z >> r;
				spheres[numsphere].x = x;
				spheres[numsphere].y = y;
				spheres[numsphere].z = z;
				spheres[numsphere].r = r;
				spheres[numsphere].center.x = x;
				spheres[numsphere].center.y = y;
				spheres[numsphere].center.z = z;
				++numsphere;//start from 0

			}
			else if (type == "pixheight")
			{
				int x;
				in >> x;
				pixheight = x;
			}
			else if (type == "pixwidth")
			{
				int x;
				in >> x;
				pixwidth = x;
			}
			else if (type == "aspect")
			{
				float x;
				in >> x;
				aspect = x;
			}
			else if (type == "fovv")
			{
				float x;
				in >> x;
				fovv = x;
			}
			else if (type == "fovh")
			{
				float x;
				in >> x;
				fovh = x;
			}
			else if (type == "viewdist")
			{
				float x;
				in >> x;
				viewdist = x;
			}
			else if (type == "light")
			{
				float x, y, z, r, g, b, w;
				in >> x >> y >> z >> w >> r >> g >> b;
				lights[numlights].x = x;
				lights[numlights].y = y;
				lights[numlights].z = z;
				lights[numlights].r = r;
				lights[numlights].g = g;
				lights[numlights].b = b;
				lights[numlights].w = w;
				if (w == 1){
					lights[numlights].lightlocation.x = x;
					lights[numlights].lightlocation.y = y;
					lights[numlights].lightlocation.z = z;
				}
				else if (w == 0){
					lights[numlights].directionallight.x = x;
					lights[numlights].directionallight.y = y;
					lights[numlights].directionallight.z = z;
				}
				++numlights;
			}
		}
	}
	source.close();
	numspherecopy = numsphere;

	for (int i = 0; i < numlights; ++i){
		if (lights[i].w == 0) {
			lights[i].directionallight = (lights[i].directionallight) / (lights[i].directionallight.GetMagnitude());
			L[i] = (lights[i].directionallight)*(-1);//store the L for directional lights
		}
	}
	//finish reading file
	//printf("%f\t%f\t%f\t%d\n",fovv,eye.y,eye.z,pixwidth);
	if (pixwidth == 0 && pixheight == 0)
	{
		cout << "need pixel height or width (any positive integer)! Abort!\n";
		return 1;
	}
	else if (pixwidth < 0 || pixheight < 0)
	{
		cout << "pixwidth or pixheight can't be negative! Abort!\n";
		return 1;
	}
	else if (pixheight != 0)
	{
		pixwidth = (int)(pixheight*aspect);
		// printf("calculate pixwidth: %d\n",pixwidth);
	}
	else //if (pixwidth!=0)
	{
		pixheight = (int)(pixwidth / aspect);
		// printf("calculate pixheight: %d\n",pixheight);
	}


	struct Point ul, ur, ll, lr;
	if (viewdir.GetMagnitude() != 1){
		viewdir = viewdir / viewdir.GetMagnitude();//nomalized viewdir 0 0 -1
	}
	struct Vec3 u, v, w;
	//printf("%f\t%f\t%f\n",w.x,w.y,w.z);//0 0 0
	w = viewdir*(-1);
	//printf("%f\t%f\t%f\n",w.x,w.y,w.z);//-0.0 -0.0 1.0

	u = viewdir%updir;//cross
	u = u / (u.GetMagnitude());//normalize u
	v = u%viewdir;
	v = v / (v.GetMagnitude());//normalize v

	float height, width;//One radian is equivalent to 180/PI degrees. height and width of view window in 3d space

	if (fovh == 0 && fovv == 0)
	{
		cout << "need fovh or fovv (positive)! Abort!\n";
		return 1;
	}
	else if (fovh > 180 || fovv > 180)
	{
		cout << "field of view can't exceed 180 degrees! Abort!\n";
		return 1;
	}
	else if (fovh == 0)
	{
		height = 2 * viewdist*tan(fovv / 2 * 3.1416 / 180.0);//10
		width = height*aspect;
		printf("height of window: %f\n width of window: %f\n", height, width);//10 10
	}
	else //if (fovv==0)
	{
		width = 2 * viewdist*tan(fovh / 2 * 3.1416 / 180.0);
		height = width / aspect;
		printf("height of window: %f\n width of window: %f\n", height, width);
	}
	/*  ul = view_origin + dn + h/2v -w/2u=view_origin+help1
	ur = view_origin + dn + h/2v + w/2u=view_origin+help2
	ll = view_origin + dn -h/2v -w/2u =..
	lr = view_origin + dn -h/2v + w/2u
	*/
	ul = eye + (viewdir*viewdist + v*(height / 2) - u*(width / 2)); //point=point+vec3
	//printPoint(ul);// -5 5 0
	ur = eye + (viewdir*viewdist + v*(height / 2) + u*(width / 2));
	//printPoint(ur);// 5 5 0
	ll = eye + (viewdir*viewdist - v*(height / 2) - u*(width / 2));
	//printPoint(ll);// -5 -5 0
	lr = eye + (viewdir*viewdist - v*(height / 2) + u*(width / 2));
	//printPoint(lr);// 5 -5 0
	Vec3 h_offset = (ur - ul) / (pixwidth - 1); //point-point=vec3. point/float=point  (10 0 0)/3=(3.333 0 0)
	Vec3 v_offset = (ll - ul) / (pixheight - 1);

	// initialize pixel array for image 
	Point coordinates[pixheight][pixwidth];//using h_offset and v_offset to calculate the coordinate of each pixel in viewing window
	//It's a vector containting 10 vectors containing 10 point. Then you can use
	ColorType colormap[pixheight][pixwidth];
	Vec3 ray_dir[pixheight][pixwidth];
	RayType ray[pixheight][pixwidth];
	
	/*	 vector<vector<Point>> coordinates(pixheight, vector<Point>(pixwidth));//using h_offset and v_offset to calculate the coordinate of each pixel in viewing window
	//It's a vector containting 10 vectors containing 10 point. Then you can use
	vector<vector<ColorType>> colormap(pixheight, vector<ColorType>(pixwidth));
	vector<vector<Vec3>> ray_dir(pixheight, vector<Vec3>(pixwidth));
	vector<vector<RayType>> ray(pixheight, vector<RayType>(pixwidth));
	*/
	//a[4][4] = p;
	for (int i = 0; i < pixheight; ++i)
	{
		for (int j = 0; j < pixwidth; ++j)
		{
			coordinates[i][j] = ul + h_offset*j + v_offset*i;//vec3+point=point
			ray_dir[i][j] = coordinates[i][j] - eye;//point-point=vec3
			ray_dir[i][j] = ray_dir[i][j] / ray_dir[i][j].GetMagnitude();//normalize all raydir
			ray[i][j].x = eye.x;//e.g. ray = (0, 0, 5) + t(0, 0, 1)
			ray[i][j].y = eye.y;
			ray[i][j].z = eye.z;
			ray[i][j].dx = ray_dir[i][j].x;
			ray[i][j].dy = ray_dir[i][j].y;
			ray[i][j].dz = ray_dir[i][j].z;
			colormap[i][j] = Trace_Ray(ray[i][j], ray_dir[i][j], eye);
		}
	}
	//ofstream out("output.ppm");
	//("G:\\Gtemp\\out.ppm")
	ofstream out("out.ppm");
	// to create a file for input: ifstream in("input.txt")    Use "ofstream out"like cout   to write to the file rather than the terminal
	out << "P3\n" << pixwidth << " " << pixheight << "\n" << 255 << "\n";
	for (int i = 0; i < pixheight; ++i)
	{
		for (int j = 0; j < pixwidth; ++j)
		{
			out << colormap[i][j].r * 255 << " " << colormap[i][j].g * 255 << " " << colormap[i][j].b * 255 << "\n";
		}
	}
	cout << "ppm file generated correctly.";
	out.close();
	return 0;
}
