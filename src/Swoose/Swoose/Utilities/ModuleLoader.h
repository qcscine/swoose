/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSEUTILITIES_MODULELOADER_H
#define SWOOSEUTILITIES_MODULELOADER_H

#include <Core/ModuleManager.h>
#include <string>
#include <vector>

namespace Scine {
namespace SwooseUtilities {

/**
 * @brief Function to load modules via the module manager (used in app and python bindings).
 * @param manager A reference to an instance of the module manager.
 * @param modules A list of module names to load.
 * @param throwError Whether to throw an error if the module could not be loaded.
 */
void loadModules(Core::ModuleManager& manager, std::vector<std::string>&& modules, bool throwError) {
  for (auto module : modules) {
    if (manager.moduleLoaded(module))
      continue;
    module.at(0) = std::tolower(module.at(0)); // convert to lowercase
    try {
      manager.load(module);
    }
    catch (const std::runtime_error& e) {
      if (throwError)
        throw std::runtime_error("Installation did not work properly. Essential modules could not be loaded.");
    }
  }
}

} // namespace SwooseUtilities
} // namespace Scine

#endif // SWOOSEUTILITIES_MODULELOADER_H
