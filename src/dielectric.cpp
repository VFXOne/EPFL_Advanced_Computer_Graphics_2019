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

#include <nori/bsdf.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
		
		// Incident angle
		float cosThetaI = Frame::cosTheta(bRec.wi);
		float Eps1 = sample[0];

		// Compute the reflected index
		float fr = fresnel(cosThetaI, m_extIOR, m_intIOR);

		// If the reflected index is bigger than Eps1 than there is reflection
		if (fr > Eps1) { 

			bRec.measure = EDiscrete;

			Normal3f normal(0.f);

			if (cosThetaI > 0) {
				bRec.eta = 1; // Light goes from int to int (reflection in the interior)
				normal.z() = 1; // The normal points upwards
			}
			else {
				bRec.eta = 1; // Light goes from ext to ext (reflection in the exterior)
				normal.z() = -1; // The normal points downwards
			}

			// Recalculate cosTheta with the right direction (cosTheta > 0)
			float cosThetaN = bRec.wi.dot(normal);

			// The output is the reflected direction
			bRec.wo = 2 * cosThetaN*normal - bRec.wi; 

			return Color3f(1.f);
		}
		else { // If the reflected index is smaller than Eps1 there is refraction

			bRec.measure = EDiscrete;

			Normal3f normal(0.f);

			if (cosThetaI > 0) {
				bRec.eta = m_extIOR / m_intIOR; // Light goes from int to ext
				normal.z() = 1; // The normal points upwards
			}
			else {
				bRec.eta = m_intIOR / m_extIOR; // Light goes from ext to int
				normal.z() = -1; // The normal points downwards
			}

			// Recalculate cosTheta with the right direction (cosTheta > 0)
			float cosThetaN = bRec.wi.dot(normal);

			// The output direction is actually the transmitted direction (Refracted here)
			bRec.wo = - bRec.eta * (bRec.wi - cosThetaN * normal) - normal * sqrt(1 - bRec.eta * bRec.eta * (1 - cosThetaN * cosThetaN));

			// Transmitted index
			float ft = bRec.eta * bRec.eta * (1 - fr); 

			return Color3f(1.f);
		}
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
