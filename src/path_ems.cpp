#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/emitter.h>
#include <nori/texture.h>
#include <math.h>

NORI_NAMESPACE_BEGIN

class PathEmsIntegrator : public Integrator {

public:
	PathEmsIntegrator(const PropertyList &props) {
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {

		Color3f L(0.f), beta(1.f), Ld(0.f);
		Ray3f new_ray(ray);
		int bounces = 0;
		float eta = 1.f;
		int maxDepth = 1000; // Very high value as we want to use Russian Roulette
		bool specularBounce = false; 

		for (bounces = 0;; ++bounces) {

			/* Intersection */
			Intersection its;
			bool foundIntersection = scene->rayIntersect(new_ray, its);
			if (!foundIntersection || bounces >= maxDepth) break;

			/* Sample new direction */
			const BSDF *bsdf = its.mesh->getBSDF();
			Vector3f wi = new_ray.d.normalized();
			BSDFQueryRecord bRec(its.shFrame.toLocal(-wi), its.uv, its);

			
			Color3f fr = bsdf->sample(bRec, sampler->next2D()); // fr = eval * cosine / pdf

			//std::cout << "before: " << its.shFrame.n.toString() << std::endl; 

			its.shFrame = bRec.si.shFrame;

			//std::cout << "after: " << its.shFrame.n.toString() << std::endl;



			if(bounces == 0 || specularBounce) {
				if (its.mesh->isEmitter() && its.shFrame.n.dot(-new_ray.d) > 0) {
					EmitterQueryRecord eRecIts(its.mesh);
					L += beta * its.mesh->getEmitter()->eval(eRecIts, its.p);
				}
			}
			
			if (!bsdf->isDiffuse()) specularBounce = true; 
			else
			{
				specularBounce = false; 
				/* Emitter sampling */
				const std::vector<Emitter *> &emitterList = scene->getEmitters();
				int randomEmitterIdx = static_cast<int>(floor(sampler->next1D() * emitterList.size()));

				while (randomEmitterIdx > emitterList.size() - 1) {
					randomEmitterIdx --;
				} 

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
				BSDFQueryRecord bRecDirect(its.shFrame.toLocal(-wi), its.shFrame.toLocal(wo), its.uv, its, ESolidAngle);
				Color3f frDirect = bsdf->eval(bRecDirect); // no sampling = no need for pdf

				L += beta * G * Ld * frDirect * emitterList.size(); 

			}

			/* Update throughput */
			eta *= bRec.eta * bRec.eta;
			beta *= fr;

			// Russian Roulette probability
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
		return "PathEmsIntegratorIntegrator[]";
	}
};

NORI_REGISTER_CLASS(PathEmsIntegrator, "path_ems");

NORI_NAMESPACE_END