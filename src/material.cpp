#include "nori/material.h"

NORI_NAMESPACE_BEGIN

float mapInterval(float X, float A, float B, float C, float D) {
	float Y = (X - A) / (B - A) * (D - C) + C;
	return Y;
}

NormalMap::NormalMap()
	: image(nullptr), width(0), height(0), comp(0)
{}

NormalMap::NormalMap(std::string filename, int channels) {
	normalmap_filename = filename;
	filesystem::path path = getFileResolver()->resolve(filesystem::path(filename));
	//stbi_set_flip_vertically_on_load(true);
	image = stbi_load(path.str().c_str(), &width, &height, &comp, channels);
}

Color3f NormalMap::getNorxel(Point2f uv) const {

	// Convert UV coordinates into texel coordinates (multiply by the width and the height of the image
	int nx = static_cast<int>(floor(uv[0] * (width - 1)));
	int ny = static_cast<int>(floor(uv[1] * (height - 1)));

	// Access the corresponding texel location to get the rgb values
	uint8_t *norxel = image + (this->comp * (ny * height + nx));
	uint8_t r = norxel[0];
	uint8_t g = norxel[1];
	uint8_t b = norxel[2];

	// Returns the corresponding texture color extracted from the texture image
	return Color3f(r, g, b) / 255.f;
}

Frame NormalMap::getShFrame(const BSDFQueryRecord &bRec) const {
	Color3f rgb = this->getNorxel(bRec.si.uv);
	Normal3f normal(mapInterval(rgb.x(), 0, 1, -1, 1), mapInterval(rgb.y(), 0, 1, -1, 1), mapInterval(rgb.z(), 0.5019607f, 1, -1, 0));

	Vector3f u = bRec.si.mesh->getVextexTangent().col(bRec.si.fidx);
	Vector3f v = bRec.si.mesh->getVextexBiTangent().col(bRec.si.fidx);
	Normal3f n = bRec.si.shFrame.n * (1 - normal.z()) + u * normal.x() + v * normal.y();

	return Frame(u, v, n);
}

Vector2i NormalMap::getSize() const {
	return Vector2i(width, height);
}

bool NormalMap::isEmpty() const {
	if (this->image == nullptr)
		return true;
	return false;
}

std::string NormalMap::toString() const {
	return tfm::format("NormalMap[filename = %s, width = %d, height = %d, comp = %d]", normalmap_filename, width, height, comp);
}

NORI_NAMESPACE_END


