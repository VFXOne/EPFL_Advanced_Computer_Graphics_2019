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

#include <nori/bitmap.h>
#include <stb_image.h>
#include <filesystem/resolver.h>
#include <filesystem/path.h>
#include <nori/frame.h>
#include "mesh.h"
#include "nori/common.h"
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class NormalMap {

private:
	uint8_t* image;
	std::string normalmap_filename;
	int width = 0;
	int height = 0;
	int comp = 0;

public:
	NormalMap();
	NormalMap(std::string filename, int channels = 3);
	Color3f getNorxel(Point2f uv) const;
	Frame getShFrame(const BSDFQueryRecord &bRec) const;
	Vector2i getSize() const;
	bool isEmpty() const;
	std::string toString() const;
};

NORI_NAMESPACE_END