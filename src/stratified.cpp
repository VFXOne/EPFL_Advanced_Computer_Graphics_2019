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

#include <nori/sampler.h>
#include <nori/block.h>
#include <pcg32.h>

NORI_NAMESPACE_BEGIN

/**
 * Stratified sampling .
 *
 */
	class Stratified : public Sampler {
	public:
		Stratified(const PropertyList &propList) {

			size_t desiredSampleCount = (size_t)propList.getInteger("sampleCount", 4);
			m_maxDimension = propList.getInteger("dimension", 4);

			size_t findResolution = 1;
			while (findResolution * findResolution < desiredSampleCount) ++findResolution;

			m_sampleCount = findResolution * findResolution;

			if (m_sampleCount != desiredSampleCount)
				std::cout<<"Sample count should be a perfect square -- rounding to " << findResolution * findResolution << std::endl;
			
			m_resolution = static_cast<int>(findResolution);
			m_sampleCount = m_resolution * m_resolution;
			m_invResolution = 1.f / static_cast<float>(m_resolution);
			m_invResolutionSquare = 1.f / static_cast<float>(m_sampleCount);

			std::cout << "res: " << m_resolution << "invRes: " << m_invResolution << " - invResSquare: " << m_invResolutionSquare << std::endl;

			m_permutations1D = new uint32_t*[m_maxDimension];
			m_permutations2D = new uint32_t*[m_maxDimension];

			for (int i = 0; i < m_maxDimension; i++) {
				m_permutations1D[i] = new uint32_t[m_sampleCount];
				m_permutations2D[i] = new uint32_t[m_sampleCount];
			}

			std::cout << this->toString() << std::endl; 
		}

		virtual ~Stratified() {
			for (int i = 0; i < m_maxDimension; i++) {
				delete[] m_permutations1D[i];
				delete[] m_permutations2D[i];
			}

			delete[] m_permutations1D;
			delete[] m_permutations2D;
		}

		std::unique_ptr<Sampler> clone() const {

			std::unique_ptr<Stratified> cloned(new Stratified());

			cloned->m_sampleCount = m_sampleCount;
			cloned->m_random = m_random;
			cloned->m_maxDimension = m_maxDimension; 
			cloned->m_resolution = m_resolution;
			cloned->m_invResolution = m_invResolution; 
			cloned->m_invResolutionSquare = m_invResolutionSquare;

			cloned->m_permutations1D = new uint32_t*[m_maxDimension];
			cloned->m_permutations2D = new uint32_t*[m_maxDimension];

			for (int i = 0; i < m_maxDimension; i++) {
				cloned->m_permutations1D[i] = new uint32_t[m_sampleCount]; 
				cloned->m_permutations2D[i] = new uint32_t[m_sampleCount];
			}

			return std::move(cloned);
		}

		void prepare(const ImageBlock &block) {
			m_random.seed(block.getOffset().x(), block.getOffset().y());
		}

		void generate() {
			for (int i = 0; i < m_maxDimension; i++) {

				for (size_t j = 0; j < m_sampleCount; j++)
					m_permutations1D[i][j] = static_cast<uint32_t>(j);
				m_random.shuffle(&m_permutations1D[i][0], &m_permutations1D[i][m_sampleCount]);

				for (size_t j = 0; j < m_sampleCount; j++)
					m_permutations2D[i][j] = static_cast<uint32_t>(j);
				m_random.shuffle(&m_permutations2D[i][0], &m_permutations2D[i][m_sampleCount]);
			}

			m_sampleIndex = 0;
			m_dimension1D = m_dimension2D = 0;
		}

		void advance() { 
			m_sampleIndex++;
			m_dimension1D = m_dimension2D = 0;
		}

		float next1D() {
			assert(m_sampleIndex < m_sampleCount); 
			if (m_dimension1D < m_maxDimension) {
				uint32_t k = m_permutations1D[m_dimension1D++][m_sampleIndex];
				float sample = (k + m_random.nextFloat()) * m_invResolutionSquare;
				return sample;
			}
			else {
				return m_random.nextFloat();
			}
		}

		Point2f next2D() {
			assert(m_sampleIndex < m_sampleCount);
			if (m_dimension2D < m_maxDimension) {
				uint32_t k = m_permutations2D[m_dimension2D++][m_sampleIndex];
				int x = k % m_resolution; 
				int y = k / m_resolution; 

				Point2f sample((x + m_random.nextFloat())*m_invResolution,
							   (y + m_random.nextFloat())*m_invResolution);
				return sample; 
			}
			else {
				return Point2f(m_random.nextFloat(), m_random.nextFloat());
			}
		}

		std::string toString() const {
			return tfm::format("Stratified[sampleCount=%i, resolution=%i, maxDimension=%i]", m_sampleCount, m_resolution, m_maxDimension);
		}
	protected:
		Stratified() { }

	private:
		pcg32 m_random;
		int m_resolution; 
		int m_maxDimension; 
		float m_invResolution, m_invResolutionSquare; 
		uint32_t **m_permutations1D, **m_permutations2D; 
		int m_dimension1D, m_dimension2D; 
		
		size_t m_sampleIndex;
};

NORI_REGISTER_CLASS(Stratified, "stratified");
NORI_NAMESPACE_END
