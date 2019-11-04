#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/emitter.h>
#include <math.h>

NORI_NAMESPACE_BEGIN

class PathMatsIntegrator : public Integrator {

public:
	PathMatsIntegrator(const PropertyList &props) {
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		
		Color3f L(0.f), beta(1.f);
		Ray3f new_ray(ray);
		int bounces = 0;
		float eta = 1.f; 
		int maxDepth = 1000; // Very high value as we want to use Russian Roulette

		for (bounces = 0;; ++bounces) {

			/* Intersection */
			Intersection its;
			bool foundIntersection = scene->rayIntersect(new_ray, its);
			if (!foundIntersection || bounces >= maxDepth) break;

			if (its.mesh->isEmitter() && its.shFrame.n.dot(-new_ray.d) > 0) {
				EmitterQueryRecord eRecIts(its.mesh);
				L += beta * its.mesh->getEmitter()->eval(eRecIts, its.p);
			}

			/* Sample new direction */
			const BSDF *bsdf = its.mesh->getBSDF();
			Vector3f wi = new_ray.d.normalized();
			BSDFQueryRecord bRec(its.shFrame.toLocal(-wi));
			Color3f fr = bsdf->sample(bRec, sampler->next2D()); // fr = eval * cosine / pdf
			
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
		return "PathMatsIntegratorIntegrator[]";
	}
};

NORI_REGISTER_CLASS(PathMatsIntegrator, "path_mats");

NORI_NAMESPACE_END