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

#pragma once

#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

struct EmitterQueryRecord {
	const Mesh *m; // Mesh associated with the emitter
	EmitterSample emSample; // Information about the sampled point on the emitter
	Vector3f wl; // Direction from the emmiter to the point

	EmitterQueryRecord()
		: m(nullptr), emSample(), wl()
	{}

	EmitterQueryRecord(const Mesh* mesh)
		: m(mesh), emSample(), wl()
	{}
};

/**
 * \brief Superclass of all emitters
 */
class Emitter : public NoriObject {
public:

	virtual Color3f sample(EmitterQueryRecord &eRec, const Point3f renderedPoint, const Point2f &sample) const = 0;

	virtual Color3f eval(const EmitterQueryRecord &eRec, const Point3f renderedPoint) const = 0;

	virtual float pdf(const EmitterQueryRecord &eRec) const = 0;

	void setMesh( Mesh *mesh) {
		m = mesh;
	}

	const Mesh* getMesh() const {
		return m;
	}

    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.) 
     * provided by this instance
     * */
    EClassType getClassType() const { return EEmitter; }
private:
	Mesh *m; 
};

NORI_NAMESPACE_END
