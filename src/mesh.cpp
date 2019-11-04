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

#include <nori/mesh.h>
#include <nori/bbox.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/warp.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

Mesh::Mesh() { }

Mesh::~Mesh() {
    delete m_bsdf;
    delete m_emitter;
}

void Mesh::activate() {
    if (!m_bsdf) {
        /* If no material was assigned, instantiate a diffuse BRDF */
        m_bsdf = static_cast<BSDF *>(
            NoriObjectFactory::createInstance("diffuse", PropertyList()));
    }

	uint32_t nbTriangles = getTriangleCount();
	m_dpdf.reserve(nbTriangles);
	m_dpdf.clear();

	float surfArea = 0;
	for (uint32_t idx = 0; idx < nbTriangles; idx++) {
		surfArea = surfaceArea(idx); 
		m_dpdf.append(surfArea);
	}
	surface = m_dpdf.normalize();

	computeTangents(); 
}

float Mesh::surfaceArea(uint32_t index) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    return 0.5f * Vector3f((p1 - p0).cross(p2 - p0)).norm();
}

bool Mesh::rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    /* Find vectors for two edges sharing v[0] */
    Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

    /* Begin calculating determinant - also used to calculate U parameter */
    Vector3f pvec = ray.d.cross(edge2);

    /* If determinant is near zero, ray lies in plane of triangle */
    float det = edge1.dot(pvec);

    if (det > -1e-8f && det < 1e-8f)
        return false;
    float inv_det = 1.0f / det;

    /* Calculate distance from v[0] to ray origin */
    Vector3f tvec = ray.o - p0;

    /* Calculate U parameter and test bounds */
    u = tvec.dot(pvec) * inv_det;
    if (u < 0.0 || u > 1.0)
        return false;

    /* Prepare to test V parameter */
    Vector3f qvec = tvec.cross(edge1);

    /* Calculate V parameter and test bounds */
    v = ray.d.dot(qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return false;

    /* Ray intersects triangle -> compute t */
    t = edge2.dot(qvec) * inv_det;

    return t >= ray.mint && t <= ray.maxt;
}

BoundingBox3f Mesh::getBoundingBox(uint32_t index) const {
    BoundingBox3f result(m_V.col(m_F(0, index)));
    result.expandBy(m_V.col(m_F(1, index)));
    result.expandBy(m_V.col(m_F(2, index)));
    return result;
}

Point3f Mesh::getCentroid(uint32_t index) const {
    return (1.0f / 3.0f) *
        (m_V.col(m_F(0, index)) +
         m_V.col(m_F(1, index)) +
         m_V.col(m_F(2, index)));
}

void Mesh::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
        case EBSDF:
            if (m_bsdf)
                throw NoriException(
                    "Mesh: tried to register multiple BSDF instances!");
            m_bsdf = static_cast<BSDF *>(obj);
            break;

        case EEmitter: {
                Emitter *emitter = static_cast<Emitter *>(obj);
                if (m_emitter)
                    throw NoriException(
                        "Mesh: tried to register multiple Emitter instances!");
                m_emitter = emitter;
            }
            break;

        default:
            throw NoriException("Mesh::addChild(<%s>) is not supported!",
                                classTypeName(obj->getClassType()));
    }
}

std::string Mesh::toString() const {
    return tfm::format(
        "Mesh[\n"
        "  name = \"%s\",\n"
        "  vertexCount = %i,\n"
        "  triangleCount = %i,\n"
        "  bsdf = %s,\n"
        "  emitter = %s\n"
        "]",
        m_name,
        m_V.cols(),
        m_F.cols(),
        m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
        m_emitter ? indent(m_emitter->toString()) : std::string("null")
    );
}

std::string Intersection::toString() const {
    if (!mesh)
        return "Intersection[invalid]";

    return tfm::format(
        "Intersection[\n"
        "  p = %s,\n"
        "  t = %f,\n"
        "  uv = %s,\n"
        "  shFrame = %s,\n"
        "  geoFrame = %s,\n"
        "  mesh = %s\n"
        "]",
        p.toString(),
        t,
        uv.toString(),
        indent(shFrame.toString()),
        indent(geoFrame.toString()),
        mesh ? mesh->toString() : std::string("null")
    );
}

std::string EmitterSample::toString() const {

	return tfm::format(
		"EmitterSample["
		" p = [%f,%f,%f],"
		" n = [%f,%f,%f],"
		" pdf = %f, "
		" isVN = %d]",
		p[0], p[1], p[2],
		n[0], n[1], n[2],
		pdf, isVertexNormal
	);
}

float Mesh::pdf() const {
	return 1/surface; 
}

EmitterSample Mesh::sampleSurface(const Point2f &sample) const {
	float Eps1 = sample[0];
	float Eps2 = sample[1];

	// Choosing randomely a triangle
	size_t tIdx = m_dpdf.sampleReuse(Eps1);

	// Vertex indices of the triangle
	uint32_t idx0 = m_F(0, tIdx), idx1 = m_F(1, tIdx), idx2 = m_F(2, tIdx);

	// Vertexes of the triangle
	Point3f x0 = m_V.col(idx0), x1 = m_V.col(idx1), x2 = m_V.col(idx2);

	float alpha = 1 - sqrt(1 - Eps1); 
	float beta = Eps2 * sqrt(1 - Eps1);

	// Sampled point as a barycenter of the triangle weighted by alpha and beta
	Point3f p = alpha * x0 + beta * x1 + (1 - alpha - beta)*x2;

	// Normal computed as the cross product of the vertexes
	Normal3f n = (x1 - x0).cross(x2 - x0).normalized();
	bool isVN = false;

	if (m_N.size() > 0) {

		// Normals of the triangle
		Normal3f nx0 = m_N.col(idx0), nx1 = m_N.col(idx1), nx2 = m_N.col(idx2);

		// Normal computed as a barycenter of the triangle normal's weighted by alpha and beta
		n = (alpha * nx0 + beta * nx1 + (1 - alpha - beta)*nx2).normalized();
	}

	// Probability density of the sample (always the same 1/whole surface mesh)
	float pdf = 1/m_dpdf.getSum();

	EmitterSample ES(p, n, isVN, pdf);

	return ES;
}

/* Added functions - Project */
void Mesh::computeTangents() {
	if (m_UV.isZero()) {
		std::cout << "Mesh::computeTangents \n In order to generate tangent values texture coordinates are required." << std::endl;
		return;
	}

	uint32_t nbTriangles = getTriangleCount();
	m_Tu.resize(3, nbTriangles);
	m_Tv.resize(3, nbTriangles);

	for (uint32_t idx = 0; idx < nbTriangles; idx++) {

		// Vertex indices of the triangle
		uint32_t idx0 = m_F(0, idx), idx1 = m_F(1, idx), idx2 = m_F(2, idx);

		// Vertexes of the triangle
		const Point3f v0 = m_V.col(idx0), v1 = m_V.col(idx1), v2 = m_V.col(idx2);

		// UV coordinates of each vertex of the triangle
		const Point2f uv0 = m_UV.col(idx0), uv1 = m_UV.col(idx1), uv2 = m_UV.col(idx2); 
		
		Vector3f dP1 = v1 - v0, dP2 = v2 - v0; // Corresponds to (dP1/dx, dP1/dy, dP1/dz) 
		Vector2f dUV1 = uv1 - uv0, dUV2 = uv2 - uv0; // Corresponds to (du,dv) 
		Normal3f n = Normal3f(dP1.cross(dP2)); 

		float determinant = dUV1.x() * dUV2.y() - dUV1.y() * dUV2.x(); 

		if (determinant != 0) {
			float invDet = 1.f / determinant;
			Vector3f dpdu = (dUV2.y() * dP1 - dUV1.y() * dP2) * invDet; // Cramer
			Vector3f dpdv = (-dUV2.x() * dP1 + dUV1.x() * dP2) * invDet; 

			m_Tu.col(idx) = dpdu; 
			m_Tv.col(idx) = dpdv; 
			//std::cout << "dpdu: " << dpdu.toString() << " - dpdv: " << dpdv.toString() << " - dot: " << dpdu.dot(dpdv) << std::endl;
		}
		else {
			//std::cout << "Generating arbitrary tangents as the determinant is 0." << std::endl;
			m_Tu.col(idx) = Vector3f(-1.f, -1.f, -1.f);
			m_Tv.col(idx) = Vector3f(-1.f, -1.f, -1.f);
		}
	}
}



NORI_NAMESPACE_END
