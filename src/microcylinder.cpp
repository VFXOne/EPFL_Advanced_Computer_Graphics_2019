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
#include <nori/texture.h>

NORI_NAMESPACE_BEGIN

int Loop = 2;

float g(float x, float sigma, float mu) {
	return ((1/(sigma*std::sqrt(2*M_PI))) * std::exp(-(x - mu) / (2 * sigma*sigma)));
}

/**
 * \brief Microcylinder / BRDF model for cloth rendering
 */
	class Microcylinder : public BSDF {
	public:
		Microcylinder(const PropertyList &propList) {

			cloth = propList.getString("cloth_type", ""); 

			std::cout << "Cloth: " << cloth << std::endl; 

			if (cloth == "Linen Plain") {
				std::cout << "Inside the if" << std::endl; 
				N = Point2i(2, 2); 

				TO.resize(2, N[0]);
				TL.resize(2, N[0]-1); 
				
				TL << 1.f, 
					  1.f; 

				TO << -25.f, 25.f,
					  -25.f, 25.f;

				std::cout << TL << std::endl; 
				std::cout << TO << std::endl; 

				intIOR = 1.46f; 
				extIOR = 1.000277f; 
				kd = 0.3f; 
				gs = 12.f;
				gv = 24.f;
				a = 0.33f;
				A = Color3f(0.2f, 0.8f, 1)*0.3f;
			}

			std::cout << this->toString() << std::endl; 
		}

		/* Added helper functions */

		Vector3f projection(const Vector3f v, const Vector3f x, const Vector3f y, const Vector3f z) const {
			return Vector3f(v.dot(x), v.dot(y), v.dot(z)); 
		}

		float surfaceReflection(const float Fr, const float PhiD, const float ThetaH) const {
			float frs = Fr * std::cos(PhiD / 2) * g(ThetaH, gs, 0);
			return frs;
		}

		Color3f volumeScattering(const float F, const float ThetaI, const float ThetaO, const float ThetaH) const {
			Color3f frv = A * F * ((1.f - kd) * g(ThetaH, gv, 0) + kd) / (std::cos(ThetaI) + std::cos(ThetaO));
			return frv; 
		}

		Color3f threadScattering(const Vector3f wi, const Vector3f wo, const Vector3f wh, const Point2f sci, const Point2f sco) const {

			//std::cout << "threadScattering" << std::endl; 
			float ThetaI = sci.x(), PhiI = sci.y();
			float ThetaO = sco.x(), PhiO = sco.y();

			float ThetaH = (ThetaI + ThetaO) / 2;
			float ThetaD = (ThetaI - ThetaO) / 2;
			float PhiD = (PhiI + PhiO) / 2; 

			float eta = extIOR / intIOR;
			float etaInv = 1 / eta; 
		
			float FrI = fresnel(wh.dot(wi), extIOR, intIOR);
			float FtI = eta * eta * (1 - FrI);

			float FrO = fresnel(wh.dot(wo), intIOR, extIOR);
			float FtO = etaInv * etaInv * (1 - FrO);

			float F = FtI * FtO;

			float frs = surfaceReflection(FrI, PhiD, ThetaH); 
			Color3f frv = volumeScattering(F, ThetaI, ThetaO, ThetaH); 

			Color3f fs = (frs + frv) / std::cos(ThetaD) * std::cos(ThetaD);

			return fs; 
		}

		float masking(const Point2f sci, const Point2f sco) const {
			float PhiI = sci.y(), PhiO = sco.y();
			float PhiD = (PhiI - PhiO) / 2;

			float MI = std::max(std::cos(PhiI), 0.f);
			float MO = std::max(std::cos(PhiO), 0.f);

			float M = (1 - g(PhiD, 1, 0)) * MI * MO + g(PhiD, 1, 0) * std::min(MI, MO);
			return M; 
		}

		float reweighting(const Point2f sci, const Point2f sco) const {

			float PsiI = sci.x(), PsiO = sco.y();
			float PsiD = (PsiI - PsiO) / 2; // Is it /2 ??

			float PI = std::max(std::cos(PsiI), 0.f);
			float PO = std::max(std::cos(PsiO), 0.f);

			float P = (1 - g(PsiD, 1, 0)) * PI * PO + g(PsiD, 1, 0) * std::min(PI, PO);

			return P; 
		}

		Color3f outgoingRadiance(const Vector3f wi, const Vector3f wo, const Normal3f n, const Vector3f t, const Vector3f nxt) const {

			// Evaluate the half way direction
			Vector3f wh = (wi + wo).normalized();

			/* Compute the longitudinal angles */

			// Projection on the local coordinates system (spherical coordinates)

			Vector3f wipB1 = projection(wi, t, nxt, n);
			Vector3f wopB1 = projection(wo, t, nxt, n);

			Vector3f wipB2 = projection(wi, n, t, nxt);
			Vector3f wopB2 = projection(wo, n, t, nxt);

			Point2f sciB1 = sphericalCoordinates(wipB1);
			Point2f scoB1 = sphericalCoordinates(wopB1);

			Point2f sciB2 = sphericalCoordinates(wipB2);
			Point2f scoB2 = sphericalCoordinates(wopB2);

			/*if (Loop > 0)
			{
				std::cout << "n: " << n.toString() << " - t:" << t.toString() << " - nxt: " << nxt.toString() << std::endl; 
				std::cout << "wi: " << wi.toString() << "wo: " << wo.toString() << "wh: " << wh.toString() << std::endl;
				std::cout << "wipB1: " << wipB1.toString() << " - wopB1: " << wopB1.toString() << std::endl;
				std::cout << "sciB1: " << sciB2.toString() << " - scoB1: " << scoB1.toString() << std::endl;
				std::cout << "wipB2: " << wipB2.toString() << " - wopB2: " << wopB2.toString() << std::endl;
				std::cout << "sciB2: " << sciB2.toString() << " - scoB2: " << scoB2.toString() << std::endl;

				Loop--; 
			}*/

			Color3f fs = threadScattering(wi, wo, wh, sciB1, scoB1);

			/* Compute Masking */
			//float M = masking(sciB1, scoB1);

			/* Compute Reweighting */
			//float P = reweighting(sciB2, sciB2);

			/* Compute Normalization */
			/*float a1 = a.x(), a2 = a.y();
			int N1 = N.x(), N2 = N.y();

			float Q = (a1 / N1) * P + (a2 / N2) * (P + (1 - a1 - a2) * wo.dot(n));

			fs *= (1 / Q) * (1 / N1) * M * P;*/

			if (fs.maxCoeff() < 0) fs = abs(fs);

			return fs; 
		}

		// direction: x = Rx, y=Ry, z=Rz
		Vector3f rotateVector(const Vector3f v, const Vector3f axis, const float offset, const std::string direction) const {

			//std::cout << "rotateVector" << std::endl; 
			float alpha = std::acos(v.dot(axis)) + offset;
			MatrixXf RxV;
			RxV.resize(3, 3);

			if (direction == "x") {
				RxV << 1.f,       0.f,         0.f,
					   0.f, cos(alpha), -sin(alpha),
					   0.f, sin(alpha),  cos(alpha);
			}
			else if (direction == "y") {
				RxV << cos(alpha), 0.f, sin(alpha),
					          0.f, 1.f,        0.f,
					  -sin(alpha), 0.f, cos(alpha);
			}
			else if (direction == "z") 
				std::cout << "Error: the rotation shouldn't be in this directions." << std::endl; 
			else 
				std::cout << "Error: the indicated direction isn't valid." << std::endl;

			return RxV * v; 
		}

		/// Evaluate the BRDF model
		Color3f eval(const BSDFQueryRecord &bRec) const {

			//std::cout << "eval" << std::endl;

			// Get the tangent and bitangent defining the direction of the threads
			/*Vector3f tu = bRec.si.shFrame.toLocal(bRec.si.mesh->getVextexTangent().col(bRec.si.fidx)).normalized();
			Vector3f tv = bRec.si.shFrame.toLocal(bRec.si.mesh->getVextexBiTangent().col(bRec.si.fidx)).normalized();*/
			/*Vector3f tu = bRec.si.mesh->getVextexTangent().col(bRec.si.fidx).normalized();
			Vector3f tv = bRec.si.mesh->getVextexBiTangent().col(bRec.si.fidx).normalized();*/

			// Get the 3D coordinate sytem n, tu, nxtu
			/*Normal3f n = bRec.si.geoFrame.n.normalized();
			Vector3f nxtu = tu.cross(n).normalized();
			Vector3f nxtv = tv.cross(n).normalized();

			Color3f fsU = Color3f(0.f, 0.f, 0.f);
			Color3f fsV = Color3f(0.f, 0.f, 0.f);*/

			/*fsV = outgoingRadiance(bRec.wi, bRec.wo, n, tv, nxtv);
			fsU = outgoingRadiance(bRec.wi, bRec.wo, n, tu, nxtu);

			for (int j = 0; j < N[0]; j++) {

				float  offset = (TO(0, j) * M_PI ) / 180;

				Vector3f new_tu = rotateVector(tu, tv, offset, "x");
				Normal3f new_n = rotateVector(n, tv, offset, "x");
				Vector3f nxtu = new_n.cross(new_tu).normalized();

				fsU += outgoingRadiance(bRec.wi, bRec.wo, new_n, new_tu, nxtu);
			}

			for (int j = 0; j < N[1]; j++) {

				float offset = (TO(1, j) * M_PI) / 180;

				Vector3f new_tv = rotateVector(tv, tu, offset, "y");
				Normal3f new_n = rotateVector(n, tu, offset, "y");
				Vector3f nxtv = new_n.cross(new_tv).normalized();

				fsV += outgoingRadiance(bRec.wi, bRec.wo, new_n, new_tv, nxtv);
			}*/


			/*fsV = outgoingRadiance(bRec.wi, bRec.wo, n, tv, nxtv);
			fsU = outgoingRadiance(bRec.wi, bRec.wo, n, tu, nxtu);

			Color3f fs = a[0] * fsU + a[1] * fsV;

			return fs; */
			return Color3f(0, 0, 0); 
		}

		/// Compute the density of \ref sample() wrt. solid angles
		float pdf(const BSDFQueryRecord &bRec) const {
			float pdf = Warp::squareToUniformHemispherePdf(bRec.wo);
			return pdf; 
		}

		/// Draw a a sample from the BRDF model
		Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {

			// Use samples between 0 and 1
			float Eps1 = sample[0];
			float Eps2 = sample[1];

			if (Frame::cosTheta(bRec.wi) <= 0)
				return Color3f(0.0f);

			// Sample an outgoing direction (reflected)
			bRec.wo = Warp::squareToUniformHemisphere(Point2f(Eps1, Eps2));
			bRec.measure = ESolidAngle;
			bRec.eta = extIOR / intIOR;

			Color3f fs = (pdf(bRec) > 0) ? eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec) : Color3f(0.f);

			return fs;
		}

		bool isDiffuse() const {
			return true;
		}

		/// Return a human-readable summary
		std::string toString() const {
			return tfm::format("Microcylinder[intIOR = %f, extIOR = %f, gs = %f, gv = %f kd = %f]",
				                intIOR, extIOR, gs, gv, kd);
		}

		EClassType getClassType() const { return EBSDF; }

	private:

		float intIOR; 
		float extIOR;
		float gs;
		float gv; 
		float kd;

		Point2i N;
		MatrixXf TO;
		MatrixXf TL; // Tangent lengths

		Point2f a; 
		Color3f A; 
		std::string cloth; 
};

NORI_REGISTER_CLASS(Microcylinder, "microcylinder");
NORI_NAMESPACE_END
