from conan import ConanFile
from conan.tools.files import copy, get, replace_in_file
from conan.tools.layout import basic_layout
import os


class cgrangesRecipe(ConanFile):
    name = "xoos-cgranges"
    homepage = "https://github.com/lh3/cgranges"
    description = "A C/C++ library for fast interval overlap queries"
    license = "MIT"
    package_type = "header-library"
    generators = "MakeDeps"

    def requirements(self):
        self.test_requires(
            "zlib-ng/2.3.3", options={"zlib_compat": True, "shared": False}
        )

    def source(self):
        get(self, **self.conan_data["sources"][self.version])

    def build(self):
        if not self.conf.get("tools.build:skip_test", default=False):
            replace_in_file(
                self,
                "test/Makefile",
                "-I../cpp",
                "-I../cpp -I${CONAN_INCLUDE_DIRS} -L${CONAN_LIB_DIRS}",
            )
            self.run(
                "cd test && make -f ../conandeps.mk -f Makefile bedcov-iitree && ./bedcov-iitree test1.bed test2.bed"
            )

    def package(self):
        copy(
            self,
            "LICENSE.txt",
            self.source_folder,
            os.path.join(self.package_folder, "licenses"),
        )
        copy(
            self,
            "*.h",
            os.path.join(self.source_folder, "cpp"),
            os.path.join(self.package_folder, "include"),
        )

    def package_info(self):
        self.cpp_info.bindirs = []
        self.cpp_info.libdirs = []

    def package_id(self):
        self.info.clear()
