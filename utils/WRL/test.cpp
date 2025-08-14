#include "weak_rate_library.h"
#include <iostream>
#include <array>
#include <cstddef>

int main() {
  std::size_t totalBytes = sizeof(gridfire::rates::weak::UNIFIED_WEAK_DATA);
  std::cout << "Total size in bytes of packed weak reaction data is : " << totalBytes/1.0e6 << " MB" << std::endl;
}
