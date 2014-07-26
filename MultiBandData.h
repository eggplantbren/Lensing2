#ifndef _MultiBandData_
#define _MultiBandData_

#include "Data.h"
#include <vector>

namespace Lensing2
{

class MultiBandData
{
	private:
		std::vector<Data> images;

	public:
		MultiBandData();

		// Same load function as regular data. Appends to the 'images'
		// vector
		void load(const char* metadata_file, const char* image_file,
				const char* sigma_file, const char* psf_file);


		static MultiBandData instance;
};

} // namespace Lensing2

#endif

