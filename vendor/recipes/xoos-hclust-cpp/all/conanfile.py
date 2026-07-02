from conan import ConanFile
from conan.tools.cmake import CMake, cmake_layout, CMakeToolchain
from conan.tools.files import copy, get
import os, shutil


class HclustCppRecipe(ConanFile):
    name = "xoos-hclust-cpp"
    homepage = "https://github.com/cdalitz/hclust-cpp"
    settings = "os", "compiler", "build_type", "arch"

    def export_sources(self):
        copy(self, "CMakeLists.txt", self.recipe_folder, self.export_sources_folder)

    def layout(self):
        cmake_layout(self, src_folder="src")

    def source(self):
        get(self, **self.conan_data["sources"][self.version])
        shutil.copy(os.path.join(self.export_sources_folder, "CMakeLists.txt"), ".")

    def generate(self):
        tc = CMakeToolchain(self)
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        copy(
            self,
            "LICENSE",
            self.source_folder,
            os.path.join(self.package_folder, "licenses"),
        )

        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["hclust-cpp"]
