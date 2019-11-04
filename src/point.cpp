#pragma once

#include <nori/emitter.h>
#include <algorithm> 
#include <stdio.h>
#include <nori/transform.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class PointLight : public Emitter {
private:
	Color3f intensity;
	Point3f position;
	Transform lightToWorld;
	Point3f pLight; 

public:

	PointLight(const PropertyList &propList) {
		intensity = propList.getColor("intensity");
		position = propList.getPoint("position");
		lightToWorld.translate(position); 
		/*without transform*/
		//pLight = position; 
		/*with transform*/
		lightToWorld * Point3f(0.f, 0.f, 0.f);
	}

	Color3f eval(const EmitterQueryRecord &eRec, const Point3f renderedPoint) const {
		return 4*M_PI*intensity;
	}

	float pdf(const EmitterQueryRecord &eRec) const {
		return 1.f; 
	}

	Color3f sample(EmitterQueryRecord &eRec, const Point3f renderedPoint, const Point2f &sample) const {
	
		eRec.wl = (pLight - renderedPoint).normalized(); 
		eRec.emSample.p = pLight; 

		return intensity / (pLight - renderedPoint).squaredNorm(); 
	}

	std::string toString() const {
		return "PointLightIntegrator[]";
	}

};

NORI_REGISTER_CLASS(PointLight, "point");
NORI_NAMESPACE_END
