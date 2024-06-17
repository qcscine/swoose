from dev.conan.base import ScineConan


class ScineSwooseConan(ScineConan):
    name = "scine_swoose"
    version = "2.0.0"
    url = "https://github.com/qcscine/swoose"
    description = """
Treat large molecular systems with self-parametrizing system-focused atomistic
models"""
    options = {
        "python": [True, False],
        "tests": [True, False],
        "coverage": [True, False],
        "microarch": ["detect", "none"],
        "database": [True, False]
    }
    default_options = {"python": False, "tests": False,
                       "coverage": False, "microarch": "none",
                       "database": False}
    exports = "dev/conan/*.py"
    exports_sources = ["dev/cmake/*", "src/*", "CMakeLists.txt", "README.rst",
                       "LICENSE.txt", "dev/conan/hook.cmake", "dev/conan/glue/*"]
    requires = ["yaml-cpp/0.6.3",
                "scine_utilities/[==9.0.0]",
                "scine_molassembler/[==2.0.1]"]
    cmake_name = "Swoose"

    def requirements(self):
        self.cmake_definitions = {
            "SWOOSE_COMPILE_DATABASE": self.options.get_safe("database")
        }

        if self.options.get_safe("database"):
            self.requires("scine_database/[==1.3.0]")

        if hasattr(super(), "requirements"):
            super().requirements()
