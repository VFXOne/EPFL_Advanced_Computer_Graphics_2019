#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
private:
	Point3f position;
	Color3f energy;

public:
	SimpleIntegrator(const PropertyList &props) {
		position = props.getPoint("position");
		energy = props.getColor("energy"); 

		std::cout << "e :" << energy << "- pos: " << position;
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		Color3f L(0.f); // Initialize the output color to 0 (black)
		Normal3f normal = its.shFrame.n; // Get the normal direction
		Point3f point = its.p; // Get the point to be rendered

		Vector3f v = position - point; // direction from the point to the position of the point light source
		
		// Starting from the origin to almost the point 
		float mint = Epsilon; 
		float maxt = 0.999f * v.norm(); // Insures that we are not hitting the point

		// Creating a new ray that goes from the light source to the current point position
		Ray3f current_ray(point, v, mint, maxt); // New ray with origin at the point and in the direction of the position of the light source

		// Check the visibility of both x and p using a shadow ray query
		if (visibility(scene, current_ray) == true) { 
			float cosTheta = normal.dot(v) / v.norm()*normal.norm(); // angle between the direction normal and the shading surface normal
			//std::cout << cosTheta << std::endl; 
			L = (energy * std::max(0.f, cosTheta)) / (4 * pow(M_PI, 2) * pow(v.norm(), 2)); 
		}

		return L; 
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
		return "SimpleIntegrator[]";
	}
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END