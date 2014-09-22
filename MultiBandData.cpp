#include "MultiBandData.h"

namespace Lensing2
{

MultiBandData MultiBandData::instance;

MultiBandData::MultiBandData()
{

}

void MultiBandData::load(const char* metadata_file, const char* image_file,
				const char* sigma_file, const char* psf_file)
{
	Data data;
	data.load(metadata_file, image_file, sigma_file, psf_file);
	images.push_back(data);
}

} // namespace Lensing2

