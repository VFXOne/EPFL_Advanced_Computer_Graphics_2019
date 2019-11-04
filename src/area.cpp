#pragma once

#include <nori/emitter.h>
#include <algorithm> 
#include <stdio.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
private: 
	Color3f radiance; 

public:

	AreaLight(const PropertyList &propList) {
		radiance = propList.getColor("radiance"); 
	}

	Color3f eval(const EmitterQueryRecord &eRec, const Point3f renderedPoint) const {
		return radiance;
	}

	float pdf(const EmitterQueryRecord &eRec) const {
		return eRec.m->pdf();
	}

	Color3f sample(EmitterQueryRecord &eRec, const Point3f renderedPoint, const Point2f &sample) const {
		eRec.emSample = eRec.m->sampleSurface(sample);
		eRec.wl = (renderedPoint - eRec.emSample.p).normalized();

		//float invPdf = std::max(abs(eRec.emSample.n.dot(eRec.wl)),0.f) / (eRec.emSample.pdf*(renderedPoint - eRec.emSample.p).squaredNorm());
		
		return (radiance * std::max(eRec.emSample.n.dot(eRec.wl), 0.f)) / ((eRec.emSample.pdf) * (renderedPoint - eRec.emSample.p).squaredNorm());
	}

	std::string toString() const {
		return "AreaLightIntegrator[]";
	}

};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END
