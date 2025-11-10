/* ***********************************************************************
//
//   Copyright (C) 2025 -- The 4D-STAR Collaboration
//   File Authors: Aaron Dotter, Emily Boudreaux
//   Last Modified: March 21, 2025
//
//   4DSSE is free software; you can use it and/or modify
//   it under the terms and restrictions the GNU General Library Public
//   License version 3 (GPLv3) as published by the Free Software Foundation.
//
//   4DSSE is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//   See the GNU Library General Public License for more details.
//
//   You should have received a copy of the GNU Library General Public License
//   along with this software; if not, write to the Free Software
//   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
// *********************************************************************** */
#include "gridfire/network.h"
#include "gridfire/reaction/reaction.h"

#include <ranges>

namespace gridfire {
    // Trim whitespace from both ends of a string
    std::string trim_whitespace(const std::string& str) {
        auto startIt = str.begin();
        const auto endIt   = str.end();

        while (startIt != endIt && std::isspace(static_cast<unsigned char>(*startIt))) {
            ++startIt;
        }
        if (startIt == endIt) {
            return "";
        }
        const auto ritr = std::find_if(str.rbegin(), std::string::const_reverse_iterator(startIt),
                                 [](const unsigned char ch){ return !std::isspace(ch); });
        return {startIt, ritr.base()};
    }

}
