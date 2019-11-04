#include "nori/texture.h"
#include "nori/common.h"
#include <filesystem/resolver.h>
#include <filesystem/path.h>

NORI_NAMESPACE_BEGIN

Texture::Texture()
	: image(nullptr), width(0), height(0), comp(0)
{}

Texture::Texture(std::string filename, int channels) {
	texture_filename = filename;
	filesystem::path path = getFileResolver()->resolve(filesystem::path(filename));
	//stbi_set_flip_vertically_on_load(true);
	image = stbi_load(path.str().c_str(), &width, &height, &comp, channels);
}

Color3f Texture::getTexel(Point2f uv) const {
	
	// Convert UV coordinates into texel coordinates (multiply by the width and the height of the image
	int tx = static_cast<int>(floor(uv[0] * (width-1)));
	int ty = static_cast<int>(floor(uv[1] * (height-1)));

	// Access the corresponding texel location to get the rgb values
	uint8_t *texel = image + (this->comp * (ty * height + tx));
	uint8_t r = texel[0];
	uint8_t g = texel[1];
	uint8_t b = texel[2];

	// Returns the corresponding texture color extracted from the texture image
	return Color3f(r, g, b) / 255.f;
}

Vector2i Texture::getSize() const {
	return Vector2i(width, height);
}

bool Texture::isEmpty() const {
	if (this->image == nullptr)
		return true; 
	return false; 	
}

std::string Texture::toString() const {
	return tfm::format("Texture[filename = %s, width = %d, height = %d, comp = %d]", texture_filename, width, height, comp);
}



NORI_NAMESPACE_END


