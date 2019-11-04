/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

static const float Pi = 3.14159265358979323846f;
static const float InvPi = 0.31830988618379067154f;
static const float Inv2Pi = 0.15915494309189533577f;
static const float Inv4Pi = 0.07957747154594766788f;

float maxValue(float value1, float value2) {
	return (value1 > value2) ? value1 : value2; 
}

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
	
	// We get the coordinates of the point
	float Eps1 = sample[0];
	float Eps2 = sample[1];

	float X = 0; 
	float Y = 0; 
	
	// Random variables
	if (Eps1 > 0 && Eps1 < 0.5f)
		X = sqrt(2.f * Eps1) - 1.f;
	else if (Eps1 > 0.5f && Eps1 < 1)
		X = -sqrt(2.f - 2.f * Eps1) + 1.f; 

	if (Eps2 > 0 && Eps2 < 0.5f)
		Y = sqrt(2.f * Eps2) - 1.f;
	else if (Eps2 > 0.5f && Eps2 < 1.f)
		Y = -sqrt(2.f - 2.f * Eps2) + 1.f;
	
	return Point2f(X, Y);
}

float Warp::squareToTentPdf(const Point2f &p) {
	
	// We get the coordinates of the point 
	float x = p[0]; 
	float y = p[1]; 

	// Set at 0 per default
	float px = 0; 
	float py = 0; 

	// We check the sign of the coordinate because of the absolute value in the PDF
	if(x > -1 && x < 1)
		px = (x >= 0) ? 1 - x : 1 + x; 

	if (y > -1 && y < 1)
		py = (y >= 0) ? 1 - y : 1 + y;

	// The probabality density function is separable: p(x,y) = px(x)py(y)
	return px * py;

    //throw NoriException("Warp::squareToTentPdf() is not yet implemented!");
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {

	// Random variables uniformely distributed on the square
	float x = sample[0];
	float y = sample[1];

	// Transformation in the new coordinate system
	float r = sqrt(x);
	float theta = 2*Pi*y;

	// Return random variables uniformely distributed on the disk
	return Point2f(r*cos(theta), r*sin(theta));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {

	// Random variables uniformely distributed on the square
	float x = p[0];
	float y = p[1];

	// We compute the radius
	float r = sqrt(pow(x, 2) + pow(y, 2)); 

	// We initialize the probability density 
	float pxy = 0;

	// We check if we are inside the disk
	if (r < 1)
		pxy = InvPi; 

	// Return the probability density function
	return pxy;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {

	// Random variables uniformely distributed on the square
	float Eps1 = sample[0]; 
	float Eps2 = sample[1]; 

	// Transformation in the new coordinate system
	float wz = 2 * Eps1 - 1.f;
	float r = sqrt(1.f - wz * wz); 
	float phi = 2 * M_PI*Eps2; 
	float wx = r * cos(phi); 
	float wy = r * sin(phi); 

	// Return random variables uniformely distributed on the sphere
	return Vector3f(wx, wy, wz);
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {

	return Inv4Pi;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
	
	// Random variables uniformely distributed on the square
	float Eps1 = sample[0];
	float Eps2 = sample[1];

	// Transformation in the new coordinate system
	float r = sqrt(std::max(0.f, 1.f - Eps1 * Eps1));
	float phi = 2 * M_PI * Eps2; 

	// Return random variables uniformely distributed on the hemisphere
	return Vector3f(r * std::cos(phi), r * std::sin(phi), Eps1);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {

	// Get the coordinate indicating the hemisphere
	float z = v[2];

	// Here we are interested in the upper hemisphere (z>0)
	return (z > 0 ) ? Inv2Pi : 0;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {

	// Transformation to get random variables uniformely distributed on the disk
	Point2f d = squareToUniformDisk(sample); 

	// We compute the z coordinate by projecting on the hemisphere
	float z = std::sqrt(maxValue(0.f, 1.f - pow(d[0],2) - pow(d[1], 2))); 

	// Return random variables distributed on the hemisphere with a cosine density function
	return Vector3f(d[0], d[1], z);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {

	// We get the coordinates information
	float x = v[0]; 
	float y = v[1]; 
	float z = v[2];

	// We compute the radius
	float r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)); 

	// We compute the angle
	float cosTheta = z/r; 

	// Here we are interested in the upper hemisphere (z>0)
	return (z >= 0) ? cosTheta * InvPi : 0; 
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
	
	// Random variables uniformely distributed on the square
	float Eps1 = sample[0];
	float Eps2 = sample[1];

	// Transformation to get random variables uniformely distributed according to the Beckman distribution
	float phi = 2 * M_PI*Eps2; 
	float cosTheta = 1 / sqrt((1 - pow(alpha, 2) * log(1 - Eps1))); 
	float sinTheta = sqrt(1 - pow(cosTheta, 2)); 
    
	// Random variables distributed according to the Beckman distribution
	return Vector3f(sinTheta*cos(phi), sinTheta*sin(phi), cosTheta); 

	//throw NoriException("Warp::squareToBeckmann() is not yet implemented!");
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
	 
	// We get the coordinates information
	float x = m[0];
	float y = m[1];
	float z = m[2];

	// We compute the radius
	float r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

	// We compute the angle
	float cosTheta = z/r; 
	float tanTheta = sqrt((1 - pow(cosTheta, 2))) / cosTheta; 

	float D = Inv2Pi * 2 * exp(-pow(tanTheta, 2) / pow(alpha, 2)) / (pow(alpha, 2)*pow(cosTheta, 3)); 

	// Upper hemisphere
	return (cosTheta > 0) ? D : 0;
}

NORI_NAMESPACE_END
