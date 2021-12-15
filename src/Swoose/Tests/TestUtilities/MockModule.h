/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SWOOSE_MOCK_MODULE_H
#define SWOOSE_MOCK_MODULE_H

#include <Core/Module.h>
#include <boost/dll/alias.hpp>
#include <memory>

namespace Scine {
namespace Swoose {

/**
 * @brief Mock module for QM/MM tests.
 */
class MockModule : public Scine::Core::Module {
 public:
  std::string name() const noexcept final;

  boost::any get(const std::string& interface, const std::string& model) const final;

  bool has(const std::string& interface, const std::string& model) const noexcept final;

  std::vector<std::string> announceInterfaces() const noexcept final;

  std::vector<std::string> announceModels(const std::string& interface) const noexcept final;

  static std::shared_ptr<Scine::Core::Module> make();
};

// Shared library entry point creating pointers to all contained modules
std::vector<std::shared_ptr<Scine::Core::Module>> moduleFactory();

} // namespace Swoose
} // namespace Scine

#ifdef __MINGW32__
/* MinGW builds are problematic. We build with default visibility, and adding
 * an attribute __dllexport__ specifically for this singular symbol leads to the
 * loss of all other weak symbols. Essentially, here we have just expanded the
 * BOOST_DLL_ALIAS macro in order to declare the type-erased const void*
 * 'moduleFactory' without any symbol visibility attribute additions that could
 * confuse the MinGW linker, which per Boost DLL documentation is unable to mix
 * weak attributes and __dllexport__ correctly.
 *
 * If ever the default visibility for this translation unit is changed, we
 * will have to revisit this bit of code for the MinGW platform again.
 *
 * Additionally, more recent Boost releases may have fixed this problem.
 * See the macro BOOST_DLL_FORCE_ALIAS_INSTANTIATIONS as used in the library's
 * example files.
 */
extern "C" {
const void* moduleFactory = reinterpret_cast<const void*>(reinterpret_cast<intptr_t>(&sample_namespace::moduleFactory));
}
#else
BOOST_DLL_ALIAS(Scine::Swoose::moduleFactory, moduleFactory)
#endif

#endif // SWOOSE_MOCK_MODULE_H
