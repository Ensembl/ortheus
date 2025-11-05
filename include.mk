binPath=${rootPath}/src/ortheus/bin
libPath=${rootPath}lib
#Modify this variable to set the location of sonLib
sonLibRootPath=${rootPath}/submodules/sonLib
sonLibPath=${sonLibRootPath}/lib

#Modify this variable to set the location of test data
export SON_TRACE_DATASETS := ${PWD}/cactusTestData

include  ${sonLibRootPath}/include.mk

cflags += -I ${sonLibPath}
basicLibs = ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a
basicLibsDependencies = ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a 
