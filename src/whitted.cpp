#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/emitter.h>
#include <math.h>

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {

public:
	WhittedIntegrator(const PropertyList &props) {
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		// We get the BSDF from the intersection found
		const BSDF *bsdf = its.mesh->getBSDF();

		Color3f C(0.0f);

		if (!bsdf->isDiffuse()) {
			Vector3f wi = ray.d.normalized();
			BSDFQueryRecord bRec(its.shFrame.toLocal(-wi));
			Color3f fr = bsdf->sample(bRec, sampler->next2D());

			if (sampler->next1D() < 0.95f) {
				// Reflected or refracted ray: origine is intersected point and direction is given by the BSDF sampling
				const Ray3f new_ray(its.p, its.shFrame.toWorld(bRec.wo)); 
				return (1 / 0.95f) * Li(scene, sampler, new_ray); // Recursive call
			}
		}

		/* Choosing an emitter and computing Le */

		// List of all the emitters in the scene
		const std::vector<Emitter *> &emitterList = scene->getEmitters();

		// Convert random emitter sample [0,1[ to [|0,L|] (L number of emitters in the scene)
		int randomEmitterIdx = static_cast<int>(floor(sampler->next1D() * emitterList.size()));

		// Select a random emitter
		Emitter *emitter = emitterList[randomEmitterIdx];

		// We create an EmmiterQueryRecord from the intersected mesh
		EmitterQueryRecord eRec(emitter->getMesh());

		// We select an emitter randomely from the mesh by using the sample function of Emitter
		// We compute Le from the Emmiter's sample function which will be adapted to the emitter type (polymorph function)
		Color3f Le = emitter->sample(eRec, its.p, sampler->next2D());

		// Check if the point being rendered belongs to an emitter
		if (its.mesh->isEmitter()) {
			if (its.shFrame.n.dot(-ray.d) > 0) {
				EmitterQueryRecord eRecIts(its.mesh);
				C += its.mesh->getEmitter()->eval(eRecIts, its.p);
			}
		}

		/* Computing fr */

		// Incident direction (emitter's position - point being rendered)
		Vector3f wo = (eRec.emSample.p - its.p).normalized();
		Vector3f wi = ray.d.normalized();

		// The direction should always be pointing outside
		BSDFQueryRecord bRec(its.shFrame.toLocal(-wi), its.shFrame.toLocal(wo), ESolidAngle);
		Color3f fr = bsdf->eval(bRec);

		/* Checking visibility */

		// Starting from the origin to almost the point 
		float mint = Epsilon;
		float maxt = (1 - Epsilon)*(eRec.emSample.p - its.p).norm(); // Insures that we are not hitting the point

		// Creating a new ray that goes from the emitter to the current point position
		Ray3f current_ray(its.p, wo, mint, maxt); // New ray with origin at the point and in the direction of the position of the light source

		// Check the visibility between currentPoint and emitter
		if (!scene->rayIntersect(current_ray)) {

			/* Computing G */
			float G = std::max(its.shFrame.n.dot(wo), 0.f);

			/* Out color */
			C += fr*G*Le;
		}

		return C;
	}

	std::string toString() const {
		return "WhittedIntegrator[]";
	}
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");

NORI_NAMESPACE_END