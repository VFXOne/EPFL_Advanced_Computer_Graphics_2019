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
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

float compute_G1(Vector3f wv, Vector3f wh, float alpha) {
	if (Frame::cosTheta(wv) == 0) return 0.f;
	float c = wv.dot(wh) / Frame::cosTheta(wv); // cosThetaV = Z coordinate of wv in local coordinate system
	float G1 = 0.f;
	if (c > 0) {
		float b = 1 / (alpha * Frame::tanTheta(wv));
		G1 = (b < 1.6f) ? ((3.555f*b + 2.181f*b*b) / (1 + 2.276f*b + 2.577f*b*b)) : 1.f;
	}
	return G1;
}

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
		
		if(Frame::cosTheta(bRec.wi) <= 0 || Frame::cosTheta(bRec.wo) <= 0)
			return Color3f(0.0f);

		// Evaluate the half way direction
		Vector3f wh = (bRec.wi + bRec.wo).normalized();

		// Evaluate the fresnel reflection coefficient
		float F = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR); 

		// Compute the PDF for the Beckmann distribution
		float D = Warp::squareToBeckmannPdf(wh, m_alpha);

		// Compute the shadowing term
		float G = compute_G1(bRec.wi,wh, m_alpha) * compute_G1(bRec.wo,wh, m_alpha); 

		// Compute cosine factors
		float cosThetaI = Frame::cosTheta(bRec.wi);
		float cosThetaO = Frame::cosTheta(bRec.wo);
		float cosThetaH = Frame::cosTheta(wh);

		Color3f color = (m_kd * INV_PI) + (m_ks*D*F*G / (4.f*cosThetaI*cosThetaO*cosThetaH));
		return color;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {

		if (Frame::cosTheta(bRec.wi) <= 0 || Frame::cosTheta(bRec.wo) <= 0)
			return 0.0f;

		Vector3f wh = (bRec.wi + bRec.wo).normalized();
		float D = Warp::squareToBeckmannPdf(wh, m_alpha);
		float Jh = 1 / (4 * wh.dot(bRec.wo));
		float pdf = m_ks * D * Jh + (1 - m_ks) * Warp::squareToCosineHemispherePdf(bRec.wo);

		return pdf;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
   
		float Eps1 = _sample[0];
		float Eps2 = _sample[1];

		if (Frame::cosTheta(bRec.wi) <= 0)
			return Color3f(0.0f);

		if (m_ks > Eps1) { // Specular
			Eps1 /= m_ks;
			Normal3f n = Warp::squareToBeckmann(Point2f(Eps1, Eps2), m_alpha);
			bRec.wo = 2.f * n * n.dot(bRec.wi) - bRec.wi; // ThetaI = ThetaO for specular

		}
		else { // Diffuse
			Eps1 = (Eps1 - m_ks) / (1.f - m_ks);
			bRec.wo = Warp::squareToCosineHemisphere(Point2f(Eps1, Eps2));
		}

		bRec.measure = ESolidAngle;
		bRec.eta = m_extIOR / m_intIOR;

        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        
		return (pdf(bRec) > 0)? eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec) : Color3f(0.f);

    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
