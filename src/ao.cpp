#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

class aoIntegrator : public Integrator {

public:
	aoIntegrator(const PropertyList &props) {
		/* Nothing to do here */
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		Color3f ambOcc(0.f); // Initialize the output color to 0 (black)
		Normal3f normal = its.shFrame.n; // Get the normal direction
		Point3f point = its.p; // Get the point to be rendered
		Point3f normalized = its.p / its.p.norm(); 

		// Sampling the point on the hemisphere
		Vector3f wHem = Warp::squareToCosineHemisphere(sampler->next2D());
		Frame f(normal); 
		Vector3f wSurf = f.toWorld(wHem); 

		// Starting from the origin to almost the point 
		float mint = Epsilon;
		float maxt = std::numeric_limits<float>::infinity();

		// Creating a new ray that goes from the light source to the current point position
		Ray3f current_ray(point, wSurf, mint, maxt); // New ray with origin at the point and in the direction of the position of the light source

		// Check the visibility of both x and p using a shadow ray query
		if (visibility(scene, current_ray) == true) {
			ambOcc = 1.f; 
		}
		return ambOcc; 
	}

	bool visibility(const Scene *scene, const Ray3f &ray) const {
		// We get the acceleration structure to go through the octree and make intersection queries
		const Accel* accel = scene->getAccel();

		Intersection its;

		// ShadowRay query to check if there's something blocking the visibility between the point and the light source
		bool pointIntersect = accel->rayIntersect(ray, its, true);

		// Return true if there's nothing in between
		return !pointIntersect;
	}

	std::string toString() const {
		return "aoIntegrator[]";
	}
};

NORI_REGISTER_CLASS(aoIntegrator, "ao");
NORI_NAMESPACE_END