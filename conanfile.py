from conans import ConanFile, CMake, tools
import os

class CHMConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
   

    name = "mesher"
    version = "1.0"
    license = "https://github.com/Chrismarsh/CHM/blob/master/LICENSE"
    author = "Chris Marsh"
    url = "https://github.com/Chrismarsh/CHM"
    description = "Canadian hydrological model"
    generators = "cmake_find_package"

    options = {"verbose_cmake":[True,False], "build_tests":[True,False] }

    default_options = {"gperftools:heapprof":True,
                       "verbose_cmake":False,
                       "build_tests":True}

    def source(self):
        git = tools.Git()
        git.clone("https://github.com/Chrismarsh/mesher.git",branch=branch)


    def requirements(self):
        self.requires( "cgal/[>=5.2]@CHM/stable" )
        self.requires( "boost/[>=1.75.0]@CHM/stable" )
        self.requires( "gdal/[>=3]@CHM/stable" )
        self.requires("proj/[>=7]@CHM/stable")

    def _configure_cmake(self):
        cmake = CMake(self)

        cmake.configure(source_folder=self.source_folder)

        return cmake

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()
        cmake.test(target="check")

    def package(self):
        cmake = self._configure_cmake()
        cmake.install()


    def imports(self):
        self.copy("*.so*", dst="lib", src="lib")  # From bin to bin
        self.copy("*.dylib*", dst="lib", src="lib")  # From lib to bin
