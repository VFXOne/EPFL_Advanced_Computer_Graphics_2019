#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/emitter.h>
#include <math.h>

NORI_NAMESPACE_BEGIN

class PathMisIntegrator : public Integrator {

public:
	PathMisIntegrator(const PropertyList &props) {
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		Color3f L(0.f), beta(1.f), Ld(0.f);
		Ray3f new_ray(ray);
		int bounces = 0;
		float eta = 1.f;
		int maxDepth = 1000;
		bool specularBounce = false;

		float prevbsdfpdf = 0.f;
		float prevcos = 0.f;

		for (bounces = 0;; ++bounces) {

			/* Intersection */
			Intersection its;
			bool foundIntersection = scene->rayIntersect(new_ray, its);
			if (!foundIntersection || bounces >= maxDepth) break;

			/* Sample new direction */
			Vector3f wi = new_ray.d.normalized();
			const BSDF *bsdf = its.mesh->getBSDF();
			BSDFQueryRecord bRec(its.shFrame.toLocal(-wi));

			Color3f fr = bsdf->sample(bRec, sampler->next2D()); // fr = eval * cosine / pdf

			prevbsdfpdf = bsdf->pdf(bRec);
			prevcos = its.shFrame.cosTheta(its.toLocal(-wi));
			
			if (its.mesh->isEmitter() && its.shFrame.n.dot(-new_ray.d) > 0) {

				/* Get emitter*/
				const Emitter *e = its.mesh->getEmitter();
				EmitterQueryRecord eRecIts(its.mesh);

				/* MIS */
				float Cosine = prevcos;
				float pLight = Cosine > 0.f ? e->pdf(eRecIts) * its.t * its.t / Cosine : 0.f;
				float pBRDF = prevbsdfpdf > 0.f ? prevbsdfpdf : 0.f;
				float wBRDF = (specularBounce || (pLight + pBRDF) == 0) ? 1.f : pBRDF / (pLight + pBRDF);

				/*std::cout << "pLight: " << pLight << std::endl; 
				std::cout << "pBRDF: " << pBRDF << std::endl; 
				std::cout << "wBRDF: " << wBRDF << std::endl;
				std::cout << "specularBounce: " << specularBounce << std::endl; */

				L += beta * e->eval(eRecIts, its.p) *wBRDF;
			}

			if (!bsdf->isDiffuse()) {
				specularBounce = true;
			}
			else
			{
				/* Emitter sampling */
				const std::vector<Emitter *> &emitterList = scene->getEmitters();
				int randomEmitterIdx = static_cast<int>(floor(sampler->next1D() * emitterList.size()));
				Emitter *emitter = emitterList[randomEmitterIdx];
				EmitterQueryRecord eRec(emitter->getMesh());
				Ld = emitter->sample(eRec, its.p, sampler->next2D());

				/* Shadow ray */
				Vector3f wo = (eRec.emSample.p - its.p).normalized();
				float mint = Epsilon;
				float maxt = (1 - Epsilon)*(eRec.emSample.p - its.p).norm();
				Ray3f current_ray(its.p, wo, mint, maxt);
				int visibility = (!scene->rayIntersect(current_ray)) ? true : false;
				float G = visibility * std::max(its.shFrame.n.dot(wo), 0.f);

				/* Contribution from emitter */
				BSDFQueryRecord bRecDirect(its.shFrame.toLocal(-wi), its.shFrame.toLocal(wo), ESolidAngle);
				Color3f frDirect = bsdf->eval(bRecDirect);

				/* MIS*/

				float Cosine = std::max(eRec.wl.dot(eRec.emSample.n), 0.f);
				float pLight = Cosine > 0.f ? eRec.emSample.pdf * (its.p - eRec.emSample.p).squaredNorm() / Cosine : 0.f;
				float pBRDF = bsdf->pdf(bRecDirect);
				float wLight = (pLight + pBRDF) != 0 ? pLight / (pLight + pBRDF) : 1.f;

				L += beta * G * Ld * frDirect * wLight * emitterList.size();
			}

			/* Update throughput */
			eta *= bRec.eta * bRec.eta;
			beta *= fr;

			/* Russian Roulette termination */
			float q = std::min(beta.maxCoeff()*eta, 0.99f);

			if (bounces > 3)
			{
				beta /= q;
				if (sampler->next1D() >= q) break;
			}

			// We update the ray for the next iteration
			new_ray = Ray3f(its.p, its.shFrame.toWorld(bRec.wo));
		}
		return L;
	}

	std::string toString() const {
		return "PathMisIntegrator[]";
	}
};

NORI_REGISTER_CLASS(PathMisIntegrator, "path_mis");

NORI_NAMESPACE_END